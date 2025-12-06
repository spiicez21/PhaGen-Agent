"""
Security module for PhaGen

Provides enterprise-grade security features:
- HSM integration for hardware-backed key storage
- Zero Data Retention (ZDR) mode for sensitive environments
- PII/PHI redaction and Data Loss Prevention (DLP)
- Tenant isolation and Role-Based Access Control (RBAC)
- End-to-end encryption (TLS 1.3 + AES-256-GCM)
- Centralized secrets management (Vault/KMS)
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
from .tenant_rbac import (
    Role,
    Permission,
    Tenant,
    TenantUser,
    TenantContext,
    TenantManager,
    require_permission,
    require_role,
    get_tenant_context,
    apply_tenant_filter,
)
from .encryption import (
    EncryptionConfig,
    DataEncryptionKey,
    EncryptionManager,
    get_tls_config,
    generate_self_signed_cert,
)
from .secrets_manager import (
    SecretsManager,
    SecretProvider,
    LocalSecretsProvider,
    VaultSecretsProvider,
    AWSSecretsProvider,
    AzureSecretsProvider,
    get_secrets_manager,
    get_secret,
    set_secret,
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
    "Role",
    "Permission",
    "Tenant",
    "TenantUser",
    "TenantContext",
    "TenantManager",
    "require_permission",
    "require_role",
    "get_tenant_context",
    "apply_tenant_filter",
    "EncryptionConfig",
    "DataEncryptionKey",
    "EncryptionManager",
    "get_tls_config",
    "generate_self_signed_cert",
    "SecretsManager",
    "SecretProvider",
    "LocalSecretsProvider",
    "VaultSecretsProvider",
    "AWSSecretsProvider",
    "AzureSecretsProvider",
    "get_secrets_manager",
    "get_secret",
    "set_secret",
]
