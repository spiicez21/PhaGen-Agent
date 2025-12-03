from __future__ import annotations

from datetime import datetime, timedelta
from typing import Dict
from uuid import uuid4

from .schemas import JobCreateRequest, JobResponse, JobStatus


class InMemoryJobStore:
    def __init__(self, ttl_minutes: int = 60) -> None:
        self._jobs: Dict[str, JobResponse] = {}
        self._ttl = timedelta(minutes=ttl_minutes)

    def create_job(self, payload: JobCreateRequest) -> JobResponse:
        now = datetime.utcnow()
        job = JobResponse(
            job_id=str(uuid4()),
            status=JobStatus.pending,
            created_at=now,
            updated_at=now,
            payload=None,
            recommendation=None,
        )
        self._jobs[job.job_id] = job
        return job

    def update_job(self, job_id: str, **kwargs) -> JobResponse:
        job = self._jobs[job_id]
        for key, value in kwargs.items():
            setattr(job, key, value)
        job.updated_at = datetime.utcnow()
        self._jobs[job_id] = job
        return job

    def get_job(self, job_id: str) -> JobResponse:
        return self._jobs[job_id]

    def sweep(self) -> None:
        now = datetime.utcnow()
        expired = [
            job_id
            for job_id, job in self._jobs.items()
            if now - job.updated_at > self._ttl
        ]
        for job_id in expired:
            self._jobs.pop(job_id, None)
