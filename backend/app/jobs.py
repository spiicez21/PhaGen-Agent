from __future__ import annotations

from datetime import datetime, timedelta, timezone
from typing import Callable, Dict, Iterable, Optional
from uuid import uuid4

from sqlalchemy import delete, func, select
from sqlalchemy.orm import Session

from .models import Document, Job as JobModel, Molecule, Passage, Report
from .schemas import JobCreateRequest, JobResponse, JobStatus
from .storage import store_raw_document


class InMemoryJobStore:
    def __init__(self, ttl_minutes: int = 60) -> None:
        self._jobs: Dict[str, JobResponse] = {}
        self._ttl = timedelta(minutes=ttl_minutes)
        self._version_index: Dict[str, int] = {}

    def create_job(self, payload: JobCreateRequest) -> JobResponse:
        now = datetime.now(timezone.utc)
        job = JobResponse(
            job_id=str(uuid4()),
            status=JobStatus.pending,
            created_at=now,
            updated_at=now,
            payload=None,
            recommendation=None,
            report_version=None,
        )
        self._jobs[job.job_id] = job
        return job

    def update_job(self, job_id: str, **kwargs) -> JobResponse:
        job = self._jobs[job_id]
        for key, value in kwargs.items():
            setattr(job, key, value)
        job.updated_at = datetime.now(timezone.utc)
        self._jobs[job_id] = job
        return job

    def assign_report_version(self, job_id: str, molecule: str | None) -> int:
        key = (molecule or "").strip().lower() or "__unlabeled__"
        next_version = self._version_index.get(key, 0) + 1
        self._version_index[key] = next_version
        job = self._jobs[job_id]
        job.report_version = next_version
        job.updated_at = datetime.now(timezone.utc)
        self._jobs[job_id] = job
        return next_version

    def get_job(self, job_id: str) -> JobResponse:
        return self._jobs[job_id]

    def sweep(self) -> None:
        now = datetime.now(timezone.utc)
        expired = [
            job_id
            for job_id, job in self._jobs.items()
            if now - job.updated_at > self._ttl
        ]
        for job_id in expired:
            self._jobs.pop(job_id, None)

    def persist_artifacts(self, job_id: str, payload: dict, report_version: int | None) -> None:
        raw_uri = store_raw_document(job_id, payload)
        job = self._jobs.get(job_id)
        if not job:
            return
        storage_meta = payload.setdefault("storage", {})
        if raw_uri:
            storage_meta["raw_document_uri"] = raw_uri
        job.payload = payload
        job.report_version = report_version or job.report_version
        job.updated_at = datetime.now(timezone.utc)
        self._jobs[job_id] = job

    def record_report_artifact(self, job_id: str, report_version: int, artifact_path: str) -> None:
        job = self._jobs.get(job_id)
        if not job:
            return
        payload = job.payload or {}
        storage_meta = payload.setdefault("storage", {})
        storage_meta["report_pdf_uri"] = artifact_path
        job.payload = payload
        job.updated_at = datetime.now(timezone.utc)
        self._jobs[job_id] = job


class PostgresJobStore:
    def __init__(self, session_factory: Callable[[], Session], ttl_minutes: int = 60) -> None:
        self._session_factory = session_factory
        self._ttl = timedelta(minutes=ttl_minutes)

    def _session(self) -> Session:
        return self._session_factory()

    def create_job(self, payload: JobCreateRequest) -> JobResponse:
        with self._session() as session:
            molecule = self._get_or_create_molecule(session, payload)
            now = datetime.now(timezone.utc)
            job = JobModel(
                id=str(uuid4()),
                molecule_id=molecule.id,
                status=JobStatus.pending.value,
                created_at=now,
                updated_at=now,
            )
            session.add(job)
            session.commit()
            session.refresh(job)
            return self._to_response(job)

    def update_job(self, job_id: str, **kwargs) -> JobResponse:
        with self._session() as session:
            job = self._get_job_or_error(session, job_id)
            for key, value in kwargs.items():
                if not hasattr(job, key):
                    continue
                if key == "status" and isinstance(value, JobStatus):
                    setattr(job, key, value.value)
                    continue
                setattr(job, key, value)
            job.updated_at = datetime.now(timezone.utc)
            session.add(job)
            session.commit()
            session.refresh(job)
            return self._to_response(job)

    def assign_report_version(self, job_id: str, molecule: str | None) -> int:
        with self._session() as session:
            job = self._get_job_or_error(session, job_id)
            max_version: Optional[int] = session.scalar(
                select(func.max(JobModel.report_version)).where(
                    JobModel.molecule_id == job.molecule_id
                )
            )
            job.report_version = (max_version or 0) + 1
            job.updated_at = datetime.now(timezone.utc)
            session.commit()
            return job.report_version

    def get_job(self, job_id: str) -> JobResponse:
        with self._session() as session:
            job = self._get_job_or_error(session, job_id)
            return self._to_response(job)

    def sweep(self) -> None:
        threshold = datetime.now(timezone.utc) - self._ttl
        with self._session() as session:
            session.execute(
                delete(JobModel).where(
                    JobModel.updated_at < threshold,
                    JobModel.status != JobStatus.running.value,
                )
            )
            session.commit()

    def persist_artifacts(self, job_id: str, payload: dict, report_version: int | None) -> None:
        workers = payload.get("workers") or {}
        with self._session() as session:
            job = self._get_job_or_error(session, job_id)

            subquery = select(Document.id).where(Document.job_id == job_id).scalar_subquery()
            session.execute(delete(Passage).where(Passage.document_id.in_(subquery)))
            session.execute(delete(Document).where(Document.job_id == job_id))
            session.execute(delete(Report).where(Report.job_id == job_id))

            for worker_name, worker_payload in workers.items():
                evidence_list: Iterable[dict] = worker_payload.get("evidence", []) or []
                for record in evidence_list:
                    document = Document(
                        job_id=job_id,
                        source=worker_name,
                        document_type=str(record.get("type") or "evidence"),
                        url=record.get("url"),
                        meta={
                            "evidence_id": record.get("evidence_id"),
                            "confidence": record.get("confidence"),
                        },
                    )
                    session.add(document)
                    session.flush()
                    passage = Passage(
                        document_id=document.id,
                        text=str(record.get("text") or ""),
                        confidence=record.get("confidence"),
                        meta={"worker": worker_name},
                    )
                    session.add(passage)

            report = Report(
                job_id=job_id,
                molecule_id=job.molecule_id,
                version=report_version or job.report_version or 1,
                payload=payload,
            )
            session.add(report)

            raw_uri = store_raw_document(job_id, payload)
            storage_meta = payload.setdefault("storage", {})
            if raw_uri:
                storage_meta["raw_document_uri"] = raw_uri
            report.payload = payload
            job.payload = payload
            session.add(job)
            session.commit()

    def record_report_artifact(self, job_id: str, report_version: int, artifact_path: str) -> None:
        with self._session() as session:
            job = self._get_job_or_error(session, job_id)
            statement = (
                select(Report)
                .where(Report.job_id == job_id, Report.version == report_version)
                .order_by(Report.created_at.desc())
                .limit(1)
            )
            report = session.scalar(statement)
            if report:
                report.artifact_path = artifact_path
                session.add(report)

            payload = job.payload or {}
            storage_meta = payload.setdefault("storage", {})
            storage_meta["report_pdf_uri"] = artifact_path
            job.payload = payload
            session.add(job)
            session.commit()

    def _get_or_create_molecule(self, session: Session, payload: JobCreateRequest) -> Molecule:
        name = payload.molecule.strip()
        stmt = select(Molecule).where(func.lower(Molecule.name) == name.lower())
        molecule = session.scalar(stmt)
        if molecule:
            return molecule

        molecule = Molecule(
            name=name,
            smiles=payload.smiles,
            synonyms=payload.synonyms or [],
        )
        session.add(molecule)
        session.flush()
        return molecule

    def _get_job_or_error(self, session: Session, job_id: str) -> JobModel:
        job = session.get(JobModel, job_id)
        if not job:
            raise KeyError(job_id)
        return job

    def _to_response(self, job: JobModel) -> JobResponse:
        return JobResponse(
            job_id=job.id,
            status=JobStatus(job.status),
            created_at=job.created_at,
            updated_at=job.updated_at,
            payload=job.payload,
            recommendation=job.recommendation,
            report_version=job.report_version,
        )
