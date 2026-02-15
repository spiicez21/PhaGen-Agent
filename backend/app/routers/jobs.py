from __future__ import annotations

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, Query, Response

from agents.exceptions import LLMError, PhaGenError, RetrievalError, WorkerError
from agents.logging_config import get_logger

from ..chemistry import build_structure_payload
from ..config import get_settings
from ..database import SessionLocal
from ..jobs import InMemoryJobStore, PostgresJobStore
from ..reporting import generate_report_pdf
from ..storage import store_report_pdf
from ..schemas import JobCreateRequest, JobResponse, JobStatus
from ..security.pii_redactor import get_dlp_policy

logger = get_logger(__name__)

try:
    from ..ml import generate_repurposing_suggestions
except (ImportError, OSError):
    generate_repurposing_suggestions = None

from agents.master import MasterAgent

settings = get_settings()
router = APIRouter(prefix=f"{settings.api_prefix}/jobs", tags=["jobs"])


# -- Dependency injection ----------------------------------------------------

def _get_job_store():
    """Provide the job store, falling back to in-memory when DB is unavailable."""
    try:
        return PostgresJobStore(SessionLocal, ttl_minutes=settings.job_ttl_minutes)
    except Exception:
        logger.warning("postgres_unavailable_falling_back_to_memory")
        return InMemoryJobStore(ttl_minutes=settings.job_ttl_minutes)


def _get_master_agent():
    return MasterAgent(top_k=settings.rag_top_k)


# Module-level defaults (used by background tasks which can't use Depends)
_job_store = _get_job_store()
_master_agent = _get_master_agent()


def _ensure_structure_metadata(job_id: str, payload: dict) -> None:
    smiles = payload.get("smiles") or payload.get("molecule_smiles")
    if not smiles:
        return

    existing = payload.get("structure") or {}
    if existing.get("svg") and not existing.get("error"):
        return

    molecule_label = (payload.get("molecule") or "molecule").strip() or "molecule"
    payload["structure"] = build_structure_payload(
        smiles=smiles,
        job_id=job_id,
        molecule_label=molecule_label,
    )


def _run_job(job_id: str, payload: JobCreateRequest) -> None:
    try:
        logger.info("job_started", job_id=job_id, molecule=payload.molecule)

        _job_store.update_job(job_id, status=JobStatus.running)
        result = _master_agent.run(
            molecule=payload.molecule,
            synonyms=payload.synonyms,
            smiles=payload.smiles,
        )
        if not result.success:
            logger.error("job_worker_failures", job_id=job_id, count=len(result.failures))
            _job_store.update_job(
                job_id,
                status=JobStatus.failed,
                payload={"failures": [failure.model_dump() for failure in result.failures]},
            )
            return
        serialized = _master_agent.serialize(result.output)
        serialized.setdefault("molecule", payload.molecule)
        if payload.smiles:
            serialized.setdefault("smiles", payload.smiles)
        version = _job_store.assign_report_version(job_id, serialized.get("molecule"))
        serialized["report_version"] = version
        _ensure_structure_metadata(job_id, serialized)

        if generate_repurposing_suggestions:
            try:
                suggestions = generate_repurposing_suggestions(job_id, serialized)
                if suggestions:
                    serialized["repurposing_suggestions"] = suggestions
            except Exception as exc:
                logger.warning("ml_suggestions_failed", job_id=job_id, error=str(exc))

        _job_store.update_job(
            job_id,
            status=JobStatus.completed,
            payload=serialized,
            recommendation=result.output.recommendation,
            report_version=version,
        )
        _job_store.persist_artifacts(job_id, serialized, version)

        logger.info(
            "job_completed",
            job_id=job_id,
            recommendation=result.output.recommendation,
            report_version=version,
        )

    except (LLMError, RetrievalError, WorkerError) as exc:
        logger.error("job_agent_error", job_id=job_id, error_code=exc.code, error=str(exc))
        _job_store.update_job(
            job_id,
            status=JobStatus.failed,
            payload={"error": str(exc), "code": exc.code},
        )
    except PhaGenError as exc:
        logger.error("job_phagen_error", job_id=job_id, error_code=exc.code, error=str(exc))
        _job_store.update_job(
            job_id,
            status=JobStatus.failed,
            payload={"error": str(exc), "code": exc.code},
        )
    except Exception as exc:
        logger.error("job_unexpected_error", job_id=job_id, error=str(exc), exc_info=True)
        _job_store.update_job(
            job_id,
            status=JobStatus.failed,
            payload={"error": str(exc)},
        )


@router.post("", response_model=JobResponse)
def create_job(
    request: JobCreateRequest, 
    background_tasks: BackgroundTasks,
    use_celery: bool = Query(default=False)
) -> JobResponse:
    job = _job_store.create_job(request)
    
    if use_celery:
        # Use Celery for distributed execution
        try:
            from ..celery_tasks import run_analysis_job
            run_analysis_job.delay(job.job_id, request.model_dump())
        except (ImportError, Exception):
            # Fall back to background task if Celery not available
            background_tasks.add_task(_run_job, job.job_id, request)
    else:
        # Use FastAPI BackgroundTasks for single-instance deployment
        background_tasks.add_task(_run_job, job.job_id, request)
    
    return job


@router.get("/{job_id}", response_model=JobResponse)
def get_job(job_id: str) -> JobResponse:
    try:
        job = _job_store.get_job(job_id)

        if job.payload:
            dlp = get_dlp_policy()
            payload_str = str(job.payload)
            allowed, reason = dlp.enforce_policy(payload_str, operation=f"get_job:{job_id}")

            if not allowed:
                logger.warning("dlp_blocked_job_retrieval", job_id=job_id)
                raise HTTPException(
                    status_code=403,
                    detail="Job contains sensitive data that cannot be exported"
                )
        
        return job
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Job not found") from exc


@router.get("/compare", response_model=list[JobResponse])
def compare_jobs(job_ids: list[str] = Query(..., min_items=2, description="Provide at least two job IDs to compare")) -> list[JobResponse]:
    jobs: list[JobResponse] = []
    missing: list[str] = []
    for job_id in job_ids:
        try:
            jobs.append(_job_store.get_job(job_id))
        except KeyError:
            missing.append(job_id)

    if missing:
        missing_str = ", ".join(missing)
        raise HTTPException(status_code=404, detail=f"Job(s) not found: {missing_str}")

    return jobs


@router.get("/{job_id}/report.pdf")
def download_report(job_id: str) -> Response:
    try:
        job = _job_store.get_job(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Job not found") from exc

    if job.status != JobStatus.completed or not job.payload:
        raise HTTPException(status_code=400, detail="Job is not ready for export")

    try:
        pdf_bytes = generate_report_pdf(job)
    except Exception as exc:  # pragma: no cover
        raise HTTPException(status_code=500, detail=f"Failed to render PDF: {exc}") from exc

    artifact_uri = store_report_pdf(job.job_id, job.report_version or 1, pdf_bytes)
    if artifact_uri:
        _job_store.record_report_artifact(job.job_id, job.report_version or 1, artifact_uri)

    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={
            "Content-Disposition": f'attachment; filename="phagen-report-{job_id}.pdf"'
        },
    )
