from __future__ import annotations

from fastapi import APIRouter, BackgroundTasks, HTTPException, Query, Response

from ..config import get_settings
from ..jobs import InMemoryJobStore
from ..reporting import generate_report_pdf
from ..schemas import JobCreateRequest, JobResponse, JobStatus

try:
    from agents.master import MasterAgent
except ModuleNotFoundError as exc:  # pragma: no cover
    raise RuntimeError(
        "Agents package not found. Ensure repository root is on PYTHONPATH."
    ) from exc

settings = get_settings()
router = APIRouter(prefix=f"{settings.api_prefix}/jobs", tags=["jobs"])
job_store = InMemoryJobStore(ttl_minutes=settings.job_ttl_minutes)
master_agent = MasterAgent(top_k=settings.rag_top_k)


def _run_job(job_id: str, payload: JobCreateRequest) -> None:
    try:
        job_store.update_job(job_id, status=JobStatus.running)
        result = master_agent.run(
            molecule=payload.molecule,
            synonyms=payload.synonyms,
            smiles=payload.smiles,
        )
        if not result.success:
            job_store.update_job(
                job_id,
                status=JobStatus.failed,
                payload={"failures": [failure.__dict__ for failure in result.failures]},
            )
            return
        serialized = master_agent.serialize(result.output)
        serialized.setdefault("molecule", payload.molecule)
        job_store.update_job(
            job_id,
            status=JobStatus.completed,
            payload=serialized,
            recommendation=result.output.recommendation,
        )
    except Exception as exc:  # pragma: no cover - logging stub
        job_store.update_job(
            job_id,
            status=JobStatus.failed,
            payload={"error": str(exc)},
        )


@router.post("", response_model=JobResponse)
def create_job(request: JobCreateRequest, background_tasks: BackgroundTasks) -> JobResponse:
    job = job_store.create_job(request)
    background_tasks.add_task(_run_job, job.job_id, request)
    return job


@router.get("/{job_id}", response_model=JobResponse)
def get_job(job_id: str) -> JobResponse:
    try:
        return job_store.get_job(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Job not found") from exc


@router.get("/compare", response_model=list[JobResponse])
def compare_jobs(job_ids: list[str] = Query(..., min_items=2, description="Provide at least two job IDs to compare")) -> list[JobResponse]:
    jobs: list[JobResponse] = []
    missing: list[str] = []
    for job_id in job_ids:
        try:
            jobs.append(job_store.get_job(job_id))
        except KeyError:
            missing.append(job_id)

    if missing:
        missing_str = ", ".join(missing)
        raise HTTPException(status_code=404, detail=f"Job(s) not found: {missing_str}")

    return jobs


@router.get("/{job_id}/report.pdf")
def download_report(job_id: str) -> Response:
    try:
        job = job_store.get_job(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Job not found") from exc

    if job.status != JobStatus.completed or not job.payload:
        raise HTTPException(status_code=400, detail="Job is not ready for export")

    try:
        pdf_bytes = generate_report_pdf(job)
    except Exception as exc:  # pragma: no cover - rendering depends on runtime libs
        raise HTTPException(status_code=500, detail=f"Failed to render PDF: {exc}") from exc

    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={
            "Content-Disposition": f'attachment; filename="phagen-report-{job_id}.pdf"'
        },
    )
