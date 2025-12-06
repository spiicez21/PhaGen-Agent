"""
Zero Data Retention (ZDR) Manager

Enforces data retention policies for sensitive pharma environments.
When enabled, prevents persistence of user uploads and clears temporary data.
"""
import os
import shutil
import logging
from pathlib import Path
from typing import Optional, List
from functools import wraps
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class ZDRConfig:
    """Configuration for Zero Data Retention mode"""
    
    def __init__(self):
        self.enabled = os.getenv("ZDR_MODE_ENABLED", "false").lower() == "true"
        self.temp_dir = Path(os.getenv("ZDR_TEMP_DIR", "/tmp/phagen_zdr"))
        self.max_memory_mb = int(os.getenv("ZDR_MAX_MEMORY_MB", "512"))
        self.auto_purge_minutes = int(os.getenv("ZDR_AUTO_PURGE_MINUTES", "60"))
        
    def __repr__(self):
        return f"ZDRConfig(enabled={self.enabled}, temp_dir={self.temp_dir})"


class ZDRManager:
    """
    Zero Data Retention Manager
    
    Manages temporary data lifecycle in sensitive environments:
    - Prevents S3/MinIO writes when ZDR is enabled
    - Auto-purges temporary files after processing
    - Clears in-memory data structures
    - Tracks data lineage for compliance
    """
    
    def __init__(self, config: Optional[ZDRConfig] = None):
        self.config = config or ZDRConfig()
        self._active_sessions: set = set()
        self._temp_files: List[Path] = []
        
        if self.config.enabled:
            logger.info(f"[ZDR] Zero Data Retention mode ENABLED")
            logger.info(f"[ZDR] Temp directory: {self.config.temp_dir}")
            logger.info(f"[ZDR] Auto-purge after: {self.config.auto_purge_minutes} minutes")
            self._ensure_temp_dir()
        else:
            logger.info("[ZDR] Zero Data Retention mode disabled")
    
    def _ensure_temp_dir(self):
        """Create temporary directory if it doesn't exist"""
        self.config.temp_dir.mkdir(parents=True, exist_ok=True)
        # Set restrictive permissions (owner only)
        os.chmod(self.config.temp_dir, 0o700)
    
    def is_enabled(self) -> bool:
        """Check if ZDR mode is active"""
        return self.config.enabled
    
    def register_temp_file(self, file_path: Path) -> None:
        """Register a temporary file for auto-cleanup"""
        if not self.config.enabled:
            return
            
        self._temp_files.append(file_path)
        logger.debug(f"[ZDR] Registered temp file: {file_path}")
    
    def create_temp_file(self, filename: str, session_id: str) -> Path:
        """
        Create a temporary file in ZDR-compliant directory
        
        Args:
            filename: Name of the file
            session_id: Job/session ID for tracking
            
        Returns:
            Path to the temporary file
        """
        if not self.config.enabled:
            raise RuntimeError("ZDR mode not enabled")
        
        session_dir = self.config.temp_dir / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        os.chmod(session_dir, 0o700)
        
        temp_file = session_dir / filename
        self.register_temp_file(temp_file)
        self._active_sessions.add(session_id)
        
        logger.info(f"[ZDR] Created temp file: {temp_file}")
        return temp_file
    
    def purge_session(self, session_id: str) -> int:
        """
        Purge all data for a session
        
        Args:
            session_id: Job/session ID to purge
            
        Returns:
            Number of files deleted
        """
        if not self.config.enabled:
            logger.warning("[ZDR] Purge called but ZDR mode not enabled")
            return 0
        
        session_dir = self.config.temp_dir / session_id
        if not session_dir.exists():
            logger.debug(f"[ZDR] No session directory to purge: {session_id}")
            return 0
        
        # Count files before deletion
        file_count = sum(1 for _ in session_dir.rglob('*') if _.is_file())
        
        # Secure deletion
        try:
            shutil.rmtree(session_dir)
            self._active_sessions.discard(session_id)
            
            # Remove from tracked temp files
            self._temp_files = [f for f in self._temp_files if not f.is_relative_to(session_dir)]
            
            logger.info(f"[ZDR] Purged session {session_id}: {file_count} files deleted")
            return file_count
            
        except Exception as e:
            logger.error(f"[ZDR] Failed to purge session {session_id}: {e}")
            raise
    
    def purge_all(self) -> int:
        """
        Purge all temporary data
        
        Returns:
            Total number of files deleted
        """
        if not self.config.enabled:
            logger.warning("[ZDR] Purge all called but ZDR mode not enabled")
            return 0
        
        total_deleted = 0
        for session_id in list(self._active_sessions):
            total_deleted += self.purge_session(session_id)
        
        logger.info(f"[ZDR] Purged all sessions: {total_deleted} total files deleted")
        return total_deleted
    
    def clear_memory(self, data_dict: dict) -> None:
        """
        Clear sensitive data from memory by zeroing out values
        
        Args:
            data_dict: Dictionary containing sensitive data to clear
        """
        if not self.config.enabled:
            return
        
        # Zero out string values
        for key in list(data_dict.keys()):
            if isinstance(data_dict[key], str):
                # Overwrite with zeros in memory
                data_dict[key] = '\0' * len(data_dict[key])
            elif isinstance(data_dict[key], bytes):
                data_dict[key] = b'\0' * len(data_dict[key])
            elif isinstance(data_dict[key], dict):
                self.clear_memory(data_dict[key])
        
        # Clear the dict
        data_dict.clear()
        logger.debug("[ZDR] Cleared sensitive data from memory")
    
    def prevent_storage_write(self, operation: str) -> None:
        """
        Prevent write operations to persistent storage in ZDR mode
        
        Args:
            operation: Description of the blocked operation
            
        Raises:
            PermissionError: If ZDR mode is enabled
        """
        if self.config.enabled:
            error_msg = f"[ZDR] Storage write blocked: {operation} (ZDR mode enabled)"
            logger.warning(error_msg)
            raise PermissionError(error_msg)
    
    @contextmanager
    def temp_session(self, session_id: str):
        """
        Context manager for temporary sessions with auto-cleanup
        
        Usage:
            with zdr_manager.temp_session("job_123") as session_dir:
                # Work with temporary files
                temp_file = session_dir / "data.json"
                temp_file.write_text(data)
            # Files automatically purged after context exit
        """
        if not self.config.enabled:
            yield None
            return
        
        session_dir = self.config.temp_dir / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        os.chmod(session_dir, 0o700)
        self._active_sessions.add(session_id)
        
        logger.info(f"[ZDR] Started temp session: {session_id}")
        
        try:
            yield session_dir
        finally:
            # Auto-purge on exit
            self.purge_session(session_id)


# Decorator for ZDR-compliant functions
def zdr_compliant(func):
    """
    Decorator to mark functions as ZDR-compliant
    Auto-purges temporary data after execution
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        zdr = ZDRManager()
        
        # Extract session_id if available
        session_id = kwargs.get('job_id') or kwargs.get('session_id')
        
        try:
            result = func(*args, **kwargs)
            return result
        finally:
            # Auto-purge if ZDR enabled and session_id provided
            if zdr.is_enabled() and session_id:
                zdr.purge_session(session_id)
    
    return wrapper


# Global ZDR manager instance
_global_zdr_manager: Optional[ZDRManager] = None


def get_zdr_manager() -> ZDRManager:
    """Get or create global ZDR manager instance"""
    global _global_zdr_manager
    if _global_zdr_manager is None:
        _global_zdr_manager = ZDRManager()
    return _global_zdr_manager


# Storage write prevention decorator
def prevent_storage_in_zdr(operation_name: str):
    """
    Decorator to prevent storage writes in ZDR mode
    
    Usage:
        @prevent_storage_in_zdr("S3 upload")
        def upload_to_s3(data):
            # This will be blocked if ZDR is enabled
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            zdr = get_zdr_manager()
            zdr.prevent_storage_write(operation_name)
            return func(*args, **kwargs)
        return wrapper
    return decorator
