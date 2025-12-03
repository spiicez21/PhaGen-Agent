from __future__ import annotations

from fastapi import APIRouter, BackgroundTasks, HTTPException

from ..config import get_settings
from ..jobs import InMemoryJobStore
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
        )
        if not result.success:
            job_store.update_job(
                job_id,
                status=JobStatus.failed,
                payload={"failures": [failure.__dict__ for failure in result.failures]},
            )
            return
        job_store.update_job(
            job_id,
            status=JobStatus.completed,
            payload=master_agent.serialize(result.output),
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
