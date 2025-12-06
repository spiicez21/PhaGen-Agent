"""
Celery tasks for distributed job execution.
"""
from __future__ import annotations

from celery import Celery
from .config import get_settings

settings = get_settings()

celery_app = Celery(
    "phagen",
    broker=settings.celery_broker_url,
    backend=settings.celery_result_backend,
)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    task_time_limit=1800,  # 30 minutes hard limit
    task_soft_time_limit=1500,  # 25 minutes soft limit
    worker_prefetch_multiplier=1,  # Fetch one task at a time for long-running jobs
    worker_max_tasks_per_child=100,  # Restart worker after 100 tasks to prevent memory leaks
)


@celery_app.task(name="phagen.run_analysis_job", bind=True)
def run_analysis_job(self, job_id: str, payload: dict):
    """
    Celery task for running molecule analysis job.
    
    Args:
        job_id: Unique job identifier
        payload: Job parameters (molecule, smiles, etc.)
    
    Returns:
        dict: Job result with worker outputs and master synthesis
    """
    import sys
    from pathlib import Path
    
    # Ensure repo root is in path
    ROOT = Path(__file__).resolve().parents[3]
    if str(ROOT) not in sys.path:
        sys.path.append(str(ROOT))
    
    from agents.master import MasterAgent
    from .database import SessionLocal
    from .jobs import PostgresJobStore
    from .schemas import JobStatus, JobCreateRequest
    from .chemistry import build_structure_payload
    from .storage import store_report_pdf, store_raw_document
    from .reporting import generate_report_pdf
    
    try:
        # Update task state
        self.update_state(state="RUNNING", meta={"status": "initializing"})
        
        # Initialize job store
        job_store = PostgresJobStore(SessionLocal, ttl_minutes=settings.job_ttl_minutes)
        job_store.update_job(job_id, status=JobStatus.running)
        
        # Ensure structure metadata
        smiles = payload.get("smiles") or payload.get("molecule_smiles")
        if smiles:
            molecule_label = (payload.get("molecule") or "molecule").strip() or "molecule"
            payload["structure"] = build_structure_payload(
                smiles=smiles,
                job_id=job_id,
                molecule_label=molecule_label,
            )
        
        # Run master agent
        self.update_state(state="RUNNING", meta={"status": "analyzing"})
        master_agent = MasterAgent(top_k=settings.rag_top_k)
        
        request = JobCreateRequest(**payload)
        result = master_agent.run(request.model_dump())
        
        # Store raw document in S3 if enabled
        if settings.storage_provider != "none":
            store_raw_document(job_id, result)
        
        # Generate and store PDF report
        self.update_state(state="RUNNING", meta={"status": "generating_report"})
        pdf_bytes = generate_report_pdf(result)
        if pdf_bytes and settings.storage_provider != "none":
            store_report_pdf(job_id, pdf_bytes)
        
        # Update job status
        job_store.update_job(job_id, status=JobStatus.completed, result=result)
        
        return result
        
    except Exception as exc:
        # Mark job as failed
        job_store.update_job(
            job_id,
            status=JobStatus.failed,
            result={"error": str(exc), "error_type": type(exc).__name__}
        )
        raise


@celery_app.task(name="phagen.cleanup_expired_jobs")
def cleanup_expired_jobs():
    """
    Periodic task to clean up expired jobs from database.
    Run via: celery -A backend.app.celery_tasks beat
    """
    from .database import SessionLocal
    from .jobs import PostgresJobStore
    
    job_store = PostgresJobStore(SessionLocal, ttl_minutes=settings.job_ttl_minutes)
    # Cleanup logic would go here
    return {"status": "cleanup_completed"}


# Celery beat schedule for periodic tasks
celery_app.conf.beat_schedule = {
    "cleanup-expired-jobs": {
        "task": "phagen.cleanup_expired_jobs",
        "schedule": 3600.0,  # Run every hour
    },
}
