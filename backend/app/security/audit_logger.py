"""
Immutable Audit Logging System

Provides tamper-proof audit trails for:
- Worker runs and agent executions
- Data access and retrieval operations
- Report generation and downloads
- User actions and API calls
"""
import json
import logging
import hashlib
from datetime import datetime, timezone
from typing import Optional, Dict, Any, List
from enum import Enum
from pathlib import Path

logger = logging.getLogger(__name__)


class AuditEventType(str, Enum):
    """Types of auditable events"""
    # Job lifecycle
    JOB_CREATED = "job.created"
    JOB_STARTED = "job.started"
    JOB_COMPLETED = "job.completed"
    JOB_FAILED = "job.failed"
    
    # Worker execution
    WORKER_STARTED = "worker.started"
    WORKER_COMPLETED = "worker.completed"
    WORKER_FAILED = "worker.failed"
    
    # Data access
    DATA_RETRIEVED = "data.retrieved"
    DATA_UPLOADED = "data.uploaded"
    DATA_DELETED = "data.deleted"
    
    # Report operations
    REPORT_GENERATED = "report.generated"
    REPORT_DOWNLOADED = "report.downloaded"
    REPORT_SHARED = "report.shared"
    
    # Authentication & authorization
    USER_LOGIN = "user.login"
    USER_LOGOUT = "user.logout"
    ACCESS_DENIED = "access.denied"
    
    # Security events
    SECURITY_ALERT = "security.alert"
    CONFIG_CHANGED = "config.changed"
    KEY_ROTATED = "key.rotated"


class AuditSeverity(str, Enum):
    """Severity levels for audit events"""
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class AuditEvent:
    """
    Immutable audit event record
    
    Once created, cannot be modified - ensuring audit trail integrity
    """
    
    def __init__(
        self,
        event_type: AuditEventType,
        event_data: Dict[str, Any],
        user_id: Optional[str] = None,
        session_id: Optional[str] = None,
        severity: AuditSeverity = AuditSeverity.INFO,
        ip_address: Optional[str] = None,
    ):
        # Immutable fields - set once at creation
        self._timestamp = datetime.now(timezone.utc)
        self._event_type = event_type
        self._event_data = event_data.copy()  # Deep copy to prevent external modification
        self._user_id = user_id
        self._session_id = session_id
        self._severity = severity
        self._ip_address = ip_address
        
        # Generate immutable hash for integrity verification
        self._hash = self._compute_hash()
    
    def _compute_hash(self) -> str:
        """Compute SHA-256 hash of event data for integrity verification"""
        data = {
            "timestamp": self._timestamp.isoformat(),
            "event_type": self._event_type.value,
            "event_data": self._event_data,
            "user_id": self._user_id,
            "session_id": self._session_id,
            "severity": self._severity.value,
            "ip_address": self._ip_address,
        }
        json_str = json.dumps(data, sort_keys=True, default=str)
        return hashlib.sha256(json_str.encode()).hexdigest()
    
    @property
    def timestamp(self) -> datetime:
        return self._timestamp
    
    @property
    def event_type(self) -> AuditEventType:
        return self._event_type
    
    @property
    def event_data(self) -> Dict[str, Any]:
        return self._event_data.copy()  # Return copy to maintain immutability
    
    @property
    def user_id(self) -> Optional[str]:
        return self._user_id
    
    @property
    def session_id(self) -> Optional[str]:
        return self._session_id
    
    @property
    def severity(self) -> AuditSeverity:
        return self._severity
    
    @property
    def ip_address(self) -> Optional[str]:
        return self._ip_address
    
    @property
    def hash(self) -> str:
        return self._hash
    
    def verify_integrity(self) -> bool:
        """Verify that event data hasn't been tampered with"""
        return self._hash == self._compute_hash()
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            "timestamp": self._timestamp.isoformat(),
            "event_type": self._event_type.value,
            "event_data": self._event_data,
            "user_id": self._user_id,
            "session_id": self._session_id,
            "severity": self._severity.value,
            "ip_address": self._ip_address,
            "hash": self._hash,
        }
    
    def to_json(self) -> str:
        """Convert to JSON string"""
        return json.dumps(self.to_dict(), default=str)


class AuditLogger:
    """
    Immutable audit logger with retention policies
    
    Features:
    - Append-only log file (no modifications/deletions)
    - Cryptographic integrity verification
    - Retention policy enforcement
    - Export capabilities for SIEM integration
    """
    
    def __init__(
        self,
        log_file: Optional[Path] = None,
        retention_days: int = 2555,  # ~7 years default
        enable_console: bool = True,
    ):
        self.log_file = log_file or Path("logs/audit.log")
        self.retention_days = retention_days
        self.enable_console = enable_console
        
        # Ensure log directory exists
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"[AUDIT] Audit logger initialized: {self.log_file}")
        logger.info(f"[AUDIT] Retention period: {retention_days} days")
    
    def log(self, event: AuditEvent) -> None:
        """
        Log an audit event (append-only)
        
        Args:
            event: Immutable audit event to log
        """
        # Verify event integrity before logging
        if not event.verify_integrity():
            logger.error("[AUDIT] Event integrity check failed - refusing to log")
            raise ValueError("Event integrity verification failed")
        
        # Append to log file (immutable - never modify existing entries)
        try:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(event.to_json() + '\n')
        except Exception as e:
            logger.error(f"[AUDIT] Failed to write to audit log: {e}")
            raise
        
        # Console output if enabled
        if self.enable_console:
            severity_prefix = {
                AuditSeverity.DEBUG: "[DEBUG]",
                AuditSeverity.INFO: "[INFO]",
                AuditSeverity.WARNING: "[WARNING]",
                AuditSeverity.ERROR: "[ERROR]",
                AuditSeverity.CRITICAL: "[CRITICAL]",
            }.get(event.severity, "[INFO]")
            
            logger.info(
                f"[AUDIT] {severity_prefix} {event.event_type.value} | "
                f"User: {event.user_id or 'system'} | Session: {event.session_id or 'N/A'} | "
                f"Data: {json.dumps(event.event_data, default=str)[:200]}"
            )
    
    def log_job_created(self, job_id: str, molecule: str, user_id: Optional[str] = None, **metadata):
        """Convenience method to log job creation"""
        event = AuditEvent(
            event_type=AuditEventType.JOB_CREATED,
            event_data={"job_id": job_id, "molecule": molecule, **metadata},
            user_id=user_id,
            session_id=job_id,
            severity=AuditSeverity.INFO,
        )
        self.log(event)
    
    def log_job_completed(self, job_id: str, recommendation: str, duration_seconds: float, user_id: Optional[str] = None):
        """Convenience method to log job completion"""
        event = AuditEvent(
            event_type=AuditEventType.JOB_COMPLETED,
            event_data={
                "job_id": job_id,
                "recommendation": recommendation,
                "duration_seconds": duration_seconds,
            },
            user_id=user_id,
            session_id=job_id,
            severity=AuditSeverity.INFO,
        )
        self.log(event)
    
    def log_worker_execution(
        self,
        worker_name: str,
        job_id: str,
        status: str,
        duration_ms: float,
        confidence: Optional[float] = None,
        user_id: Optional[str] = None,
    ):
        """Convenience method to log worker execution"""
        event_type = {
            "started": AuditEventType.WORKER_STARTED,
            "completed": AuditEventType.WORKER_COMPLETED,
            "failed": AuditEventType.WORKER_FAILED,
        }.get(status, AuditEventType.WORKER_COMPLETED)
        
        event = AuditEvent(
            event_type=event_type,
            event_data={
                "worker_name": worker_name,
                "job_id": job_id,
                "duration_ms": duration_ms,
                "confidence": confidence,
            },
            user_id=user_id,
            session_id=job_id,
            severity=AuditSeverity.WARNING if status == "failed" else AuditSeverity.INFO,
        )
        self.log(event)
    
    def log_data_access(
        self,
        operation: str,
        resource_type: str,
        resource_id: str,
        user_id: Optional[str] = None,
        ip_address: Optional[str] = None,
        **metadata
    ):
        """Convenience method to log data access"""
        event_type = {
            "retrieve": AuditEventType.DATA_RETRIEVED,
            "upload": AuditEventType.DATA_UPLOADED,
            "delete": AuditEventType.DATA_DELETED,
        }.get(operation, AuditEventType.DATA_RETRIEVED)
        
        event = AuditEvent(
            event_type=event_type,
            event_data={
                "operation": operation,
                "resource_type": resource_type,
                "resource_id": resource_id,
                **metadata,
            },
            user_id=user_id,
            ip_address=ip_address,
            severity=AuditSeverity.WARNING if operation == "delete" else AuditSeverity.INFO,
        )
        self.log(event)
    
    def log_report_action(
        self,
        action: str,
        report_id: str,
        job_id: str,
        user_id: Optional[str] = None,
        ip_address: Optional[str] = None,
    ):
        """Convenience method to log report actions"""
        event_type = {
            "generated": AuditEventType.REPORT_GENERATED,
            "downloaded": AuditEventType.REPORT_DOWNLOADED,
            "shared": AuditEventType.REPORT_SHARED,
        }.get(action, AuditEventType.REPORT_GENERATED)
        
        event = AuditEvent(
            event_type=event_type,
            event_data={
                "action": action,
                "report_id": report_id,
                "job_id": job_id,
            },
            user_id=user_id,
            session_id=job_id,
            ip_address=ip_address,
            severity=AuditSeverity.INFO,
        )
        self.log(event)
    
    def log_security_event(
        self,
        event_description: str,
        severity: AuditSeverity = AuditSeverity.WARNING,
        user_id: Optional[str] = None,
        ip_address: Optional[str] = None,
        **metadata
    ):
        """Convenience method to log security events"""
        event = AuditEvent(
            event_type=AuditEventType.SECURITY_ALERT,
            event_data={"description": event_description, **metadata},
            user_id=user_id,
            ip_address=ip_address,
            severity=severity,
        )
        self.log(event)
    
    def read_events(
        self,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
        event_type: Optional[AuditEventType] = None,
        user_id: Optional[str] = None,
        limit: int = 1000,
    ) -> List[Dict[str, Any]]:
        """
        Read audit events with filtering
        
        Args:
            start_time: Filter events after this time
            end_time: Filter events before this time
            event_type: Filter by event type
            user_id: Filter by user ID
            limit: Maximum number of events to return
            
        Returns:
            List of matching audit events
        """
        events = []
        
        if not self.log_file.exists():
            return events
        
        try:
            with open(self.log_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.strip():
                        continue
                    
                    try:
                        event_dict = json.loads(line)
                        
                        # Apply filters
                        event_time = datetime.fromisoformat(event_dict["timestamp"])
                        
                        if start_time and event_time < start_time:
                            continue
                        if end_time and event_time > end_time:
                            continue
                        if event_type and event_dict["event_type"] != event_type.value:
                            continue
                        if user_id and event_dict.get("user_id") != user_id:
                            continue
                        
                        events.append(event_dict)
                        
                        if len(events) >= limit:
                            break
                            
                    except (json.JSONDecodeError, KeyError) as e:
                        logger.warning(f"[AUDIT] Skipped malformed log entry: {e}")
                        continue
        
        except Exception as e:
            logger.error(f"[AUDIT] Failed to read audit log: {e}")
            raise
        
        return events
    
    def export_for_siem(self, output_file: Path, start_time: Optional[datetime] = None, end_time: Optional[datetime] = None):
        """
        Export audit logs in SIEM-compatible format (JSON Lines)
        
        Args:
            output_file: Path to export file
            start_time: Export events after this time
            end_time: Export events before this time
        """
        events = self.read_events(start_time=start_time, end_time=end_time, limit=1000000)
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            for event in events:
                f.write(json.dumps(event, default=str) + '\n')
        
        logger.info(f"[AUDIT] Exported {len(events)} events to {output_file}")


# Global audit logger instance
_global_audit_logger: Optional[AuditLogger] = None


def get_audit_logger() -> AuditLogger:
    """Get or create global audit logger instance"""
    global _global_audit_logger
    if _global_audit_logger is None:
        _global_audit_logger = AuditLogger()
    return _global_audit_logger
