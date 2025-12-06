"""
Security module for PhaGen

Provides enterprise-grade security features:
- HSM integration for hardware-backed key storage
- Zero Data Retention (ZDR) mode for sensitive environments
- PII/PHI redaction and Data Loss Prevention (DLP)
"""
from .hsm_manager import HSMManager, HSMConfig
from .zdr_manager import ZDRManager, ZDRConfig, get_zdr_manager, zdr_compliant, prevent_storage_in_zdr
from .pii_redactor import (
    PIIRedactor,
    PIIType,
    DLPPolicy,
    get_pii_redactor,
    get_dlp_policy,
    redact_pii
)

__all__ = [
    "HSMManager",
    "HSMConfig",
    "ZDRManager",
    "ZDRConfig",
    "get_zdr_manager",
    "zdr_compliant",
    "prevent_storage_in_zdr",
    "PIIRedactor",
    "PIIType",
    "DLPPolicy",
    "get_pii_redactor",
    "get_dlp_policy",
    "redact_pii",
]
