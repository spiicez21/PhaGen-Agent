from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from typing import Any, Optional

from .config import get_settings

try:  # pragma: no cover - optional heavy dependency
    import boto3
    from botocore.client import Config as BotoConfig
    from botocore.exceptions import BotoCoreError, ClientError
except Exception as exc:  # noqa: BLE001
    boto3 = None  # type: ignore[assignment]
    BotoConfig = None  # type: ignore[assignment]
    BotoError = exc
else:
    BotoError = None

_logger = logging.getLogger(__name__)
_settings = get_settings()
_storage_client = None


def _is_enabled() -> bool:
    return (
        _settings.storage_provider.lower() == "s3"
        and boto3 is not None
        and BotoError is None
    )


def _client():  # type: ignore[return-type]
    global _storage_client
    if not _is_enabled():
        return None
    if _storage_client is None:
        _storage_client = boto3.client(  # type: ignore[call-arg]
            "s3",
            region_name=_settings.s3_region,
            endpoint_url=_settings.s3_endpoint_url,
            aws_access_key_id=_settings.s3_access_key,
            aws_secret_access_key=_settings.s3_secret_key,
            use_ssl=_settings.s3_use_ssl,
            config=BotoConfig(signature_version="s3v4") if BotoConfig else None,
        )
    return _storage_client


def _ensure_bucket(client, bucket: str) -> None:
    if not bucket:
        raise ValueError("Bucket name must be provided")
    try:
        client.head_bucket(Bucket=bucket)
        return
    except ClientError as exc:  # pragma: no cover - network side effects
        error_code = exc.response.get("Error", {}).get("Code", "")
        if error_code not in {"404", "NoSuchBucket"}:
            raise
    create_params: dict[str, Any] = {"Bucket": bucket}
    if _settings.s3_region != "us-east-1":
        create_params["CreateBucketConfiguration"] = {
            "LocationConstraint": _settings.s3_region,
        }
    client.create_bucket(**create_params)


def _put_object(*, bucket: str, key: str, body: bytes, content_type: str) -> str | None:
    client = _client()
    if not client:
        return None
    try:
        _ensure_bucket(client, bucket)
        client.put_object(
            Bucket=bucket,
            Key=key,
            Body=body,
            ContentType=content_type,
        )
        return f"s3://{bucket}/{key}"
    except (BotoCoreError, ClientError, ValueError) as exc:  # pragma: no cover
        _logger.error("Failed to upload %s: %s", key, exc)
        return None


def store_raw_document(job_id: str, payload: dict) -> Optional[str]:
    if not _is_enabled():
        return None
    timestamp = datetime.utcnow().replace(tzinfo=timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    key = f"jobs/{job_id}/raw/{timestamp}.json"
    body = json.dumps(payload, ensure_ascii=False, default=str).encode("utf-8")
    return _put_object(
        bucket=_settings.s3_raw_documents_bucket,
        key=key,
        body=body,
        content_type="application/json",
    )


def store_report_pdf(job_id: str, report_version: int, pdf_bytes: bytes) -> Optional[str]:
    if not _is_enabled():
        return None
    timestamp = datetime.utcnow().replace(tzinfo=timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    key = f"jobs/{job_id}/reports/v{report_version}-{timestamp}.pdf"
    return _put_object(
        bucket=_settings.s3_reports_bucket,
        key=key,
        body=pdf_bytes,
        content_type="application/pdf",
    )


__all__ = ["store_raw_document", "store_report_pdf"]
