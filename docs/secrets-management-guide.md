# Centralized Secrets Management Guide

## Overview

PhaGen implements centralized secrets management to securely store and rotate sensitive configuration data, API keys, credentials, and certificates. The system supports multiple backend providers with automatic rotation and versioning.

## Features

### 🔐 Multi-Provider Support
- **HashiCorp Vault**: Industry-standard secrets management
- **AWS Secrets Manager**: Native AWS integration
- **Azure Key Vault**: Azure-native secrets storage
- **Local Provider**: Development/testing only

### 🔄 Automatic Rotation
- Configurable rotation schedules
- Zero-downtime rotation with dual keys
- Version history for rollback
- Audit logging for all operations

### ⚡ Performance
- In-memory caching with TTL
- Connection pooling
- Batch operations
- Cache invalidation on updates

### 🎯 Use Cases
- Database credentials
- API keys (NCBI, OpenFDA, PubChem, LLM providers)
- OAuth tokens
- TLS certificates
- Encryption keys
- Service account credentials

## Architecture

### Secrets Hierarchy

```
┌─────────────────────────────────────────────────────────┐
│                                                         │
│  Application                                            │
│       ↓                                                 │
│  SecretsManager (unified interface)                     │
│       ↓                                                 │
│  SecretProvider (Vault/AWS/Azure/Local)                 │
│       ↓                                                 │
│  Backend Storage (encrypted at rest)                    │
│       ↓                                                 │
│  Audit Logs                                             │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

### Secret Versioning

```
Secret: db_password
├── Version 1: "password123" (archived)
├── Version 2: "password456" (archived)
└── Version 3: "password789" (current)
```

## Configuration

### Environment Variables

```bash
# Secrets provider: local, vault, aws, azure
SECRETS_PROVIDER=vault

# Cache TTL (seconds)
SECRETS_CACHE_TTL=300  # 5 minutes
```

### HashiCorp Vault

```bash
# Vault configuration
VAULT_ADDR=https://vault.example.com:8200
VAULT_TOKEN=s.xxxxxxxxxxxxxxxxxxxxxxxxx

# Or use Kubernetes auth
VAULT_ROLE=phagen-backend
VAULT_MOUNT_POINT=secret
```

### AWS Secrets Manager

```bash
# AWS credentials
AWS_REGION=us-east-1
AWS_ACCESS_KEY_ID=<your-key>
AWS_SECRET_ACCESS_KEY=<your-secret>

# Or use IAM role (EC2, ECS, Lambda)
```

### Azure Key Vault

```bash
# Azure credentials
AZURE_TENANT_ID=<tenant-id>
AZURE_CLIENT_ID=<client-id>
AZURE_CLIENT_SECRET=<client-secret>

# Key vault URL
AZURE_KEY_VAULT_URL=https://phagen-vault.vault.azure.net/
```

## Usage

### Basic Operations

```python
from app.security import get_secrets_manager

manager = get_secrets_manager()

# Set secret
manager.set("db_password", "secure_password_123")

# Get secret
password = manager.get("db_password")

# List secrets
secrets = manager.list()
# Returns: ["db_password", "api_key_ncbi", "tls_cert_path", ...]

# Delete secret
manager.delete("old_api_key")
```

### Secret Rotation

```python
# Rotate single secret
new_value = manager.rotate("api_key")

# Rotate all API keys
manager.rotate_all(pattern="api_key")

# Rotate database credentials
new_password = manager.rotate("db_password")
# Update database connection with new password
```

### Convenience Functions

```python
from app.security import get_secret, set_secret

# Quick get
api_key = get_secret("ncbi_api_key")

# Quick set
set_secret("new_token", "sk-...")
```

### Database Connection Example

```python
from app.security import get_secrets_manager

manager = get_secrets_manager()

# Retrieve credentials
DB_HOST = manager.get("db_host")
DB_PORT = manager.get("db_port")
DB_USER = manager.get("db_user")
DB_PASSWORD = manager.get("db_password")
DB_NAME = manager.get("db_name")

# Build connection string
DATABASE_URL = f"postgresql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"
```

### API Key Management

```python
# Store API keys
manager.set("ncbi_api_key", "your-ncbi-key")
manager.set("openfda_api_key", "your-fda-key")
manager.set("openai_api_key", "sk-...")

# Retrieve in application
ncbi_key = manager.get("ncbi_api_key")
response = requests.get(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    params={"api_key": ncbi_key, ...}
)
```

## HashiCorp Vault Integration

### Setup

1. **Install Vault CLI**:
   ```bash
   brew install vault  # macOS
   # or
   apt-get install vault  # Ubuntu
   ```

2. **Start Vault Server** (dev mode):
   ```bash
   vault server -dev
   ```

3. **Enable KV v2 Secrets Engine**:
   ```bash
   vault secrets enable -path=secret kv-v2
   ```

4. **Set Secrets**:
   ```bash
   vault kv put secret/db_password value="secure_password"
   vault kv put secret/api_keys/ncbi value="ncbi-key-123"
   ```

### Production Setup

```bash
# Enable Vault with TLS
vault server -config=/etc/vault/config.hcl

# Initialize and unseal
vault operator init
vault operator unseal <key1>
vault operator unseal <key2>
vault operator unseal <key3>

# Create policy for PhaGen
vault policy write phagen-policy phagen-policy.hcl

# Enable AppRole auth
vault auth enable approle
vault write auth/approle/role/phagen \
    secret_id_ttl=24h \
    token_ttl=1h \
    token_max_ttl=4h \
    policies="phagen-policy"
```

**Policy File** (`phagen-policy.hcl`):
```hcl
# Allow read/write to phagen secrets
path "secret/data/phagen/*" {
  capabilities = ["create", "read", "update", "delete", "list"]
}

# Allow secret rotation
path "secret/metadata/phagen/*" {
  capabilities = ["list"]
}
```

## AWS Secrets Manager Integration

### Setup

1. **Create Secret**:
   ```bash
   aws secretsmanager create-secret \
       --name phagen/db_password \
       --secret-string "secure_password_123"
   ```

2. **Enable Rotation**:
   ```bash
   aws secretsmanager rotate-secret \
       --secret-id phagen/db_password \
       --rotation-lambda-arn arn:aws:lambda:us-east-1:123456789012:function:phagen-rotation
   ```

3. **IAM Policy**:
   ```json
   {
     "Version": "2012-10-17",
     "Statement": [
       {
         "Effect": "Allow",
         "Action": [
           "secretsmanager:GetSecretValue",
           "secretsmanager:DescribeSecret",
           "secretsmanager:ListSecrets",
           "secretsmanager:PutSecretValue",
           "secretsmanager:UpdateSecret"
         ],
         "Resource": "arn:aws:secretsmanager:us-east-1:123456789012:secret:phagen/*"
       }
     ]
   }
   ```

## Azure Key Vault Integration

### Setup

1. **Create Key Vault**:
   ```bash
   az keyvault create \
       --name phagen-vault \
       --resource-group phagen-rg \
       --location eastus
   ```

2. **Set Secret**:
   ```bash
   az keyvault secret set \
       --vault-name phagen-vault \
       --name db-password \
       --value "secure_password_123"
   ```

3. **Grant Access**:
   ```bash
   az keyvault set-policy \
       --name phagen-vault \
       --object-id <app-object-id> \
       --secret-permissions get list set delete
   ```

## Secret Rotation

### Automatic Rotation Schedule

```python
from celery import Celery
from app.security import get_secrets_manager

app = Celery('phagen')

@app.task
def rotate_api_keys():
    """Rotate API keys monthly"""
    manager = get_secrets_manager()
    manager.rotate_all(pattern="api_key")

@app.task
def rotate_db_credentials():
    """Rotate database credentials quarterly"""
    manager = get_secrets_manager()
    
    # Rotate password
    new_password = manager.rotate("db_password")
    
    # Update database user
    update_database_password(new_password)

# Schedule tasks
app.conf.beat_schedule = {
    'rotate-api-keys': {
        'task': 'app.tasks.rotate_api_keys',
        'schedule': 30 * 86400  # 30 days
    },
    'rotate-db-credentials': {
        'task': 'app.tasks.rotate_db_credentials',
        'schedule': 90 * 86400  # 90 days
    }
}
```

### Zero-Downtime Rotation

```python
def rotate_database_credentials():
    """Rotate DB password with zero downtime"""
    manager = get_secrets_manager()
    
    # 1. Generate new password
    new_password = manager.rotate("db_password")
    
    # 2. Create new database user with new password
    db.execute(f"CREATE USER phagen_new WITH PASSWORD '{new_password}'")
    db.execute("GRANT ALL PRIVILEGES ON DATABASE phagen_db TO phagen_new")
    
    # 3. Update application config (rolling restart)
    update_app_config("db_user", "phagen_new")
    update_app_config("db_password", new_password)
    
    # 4. Wait for all instances to restart
    time.sleep(60)
    
    # 5. Drop old user
    db.execute("DROP USER phagen_old")
```

## Security Best Practices

### Production Checklist

- [ ] Use managed secrets service (Vault/AWS/Azure)
- [ ] Never commit secrets to Git
- [ ] Rotate secrets regularly (30-90 days)
- [ ] Use least-privilege IAM policies
- [ ] Enable audit logging for all secret access
- [ ] Encrypt secrets at rest in backend
- [ ] Use TLS for secrets transit
- [ ] Implement secret versioning
- [ ] Document recovery procedures
- [ ] Test rotation procedures in staging
- [ ] Monitor for unauthorized access
- [ ] Use separate secrets per environment

### Access Control

```python
# Role-based access to secrets
from app.security import get_secrets_manager, require_role, Role

@require_role(Role.SUPER_ADMIN)
def rotate_all_secrets():
    """Only super admins can rotate secrets"""
    manager = get_secrets_manager()
    manager.rotate_all()

@require_role(Role.ANALYST)
def get_api_keys():
    """Analysts can read API keys"""
    manager = get_secrets_manager()
    return {
        "ncbi": manager.get("api_key_ncbi"),
        "openfda": manager.get("api_key_openfda")
    }
```

## Compliance Mapping

| Requirement | Implementation | Evidence |
|------------|----------------|----------|
| **SOC2 CC6.1** | Secure secrets storage | `secrets_manager.py` L1-455 |
| **ISO 27001 A.9.4.1** | Access control to secrets | RBAC integration |
| **ISO 27001 A.10.1.2** | Key management | Version control + rotation |
| **PCI DSS 3.5** | Protect stored keys | Encryption at rest |
| **HIPAA §164.312(a)(2)(iv)** | Key management procedures | This document |

## Troubleshooting

### Common Issues

#### 1. "Secret not found"

**Cause**: Secret never created or deleted  
**Solution**: Create secret or check spelling

```python
# List all secrets
secrets = manager.list()
print(secrets)

# Create missing secret
manager.set("missing_secret", "value")
```

#### 2. "Vault authentication failed"

**Cause**: Invalid token or expired credentials  
**Solution**: Refresh token or re-authenticate

```bash
# Check Vault status
vault status

# Login again
vault login <token>
```

#### 3. "AWS credentials not found"

**Cause**: Missing ~/.aws/credentials or IAM role  
**Solution**: Configure AWS credentials

```bash
aws configure
# or
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
```

#### 4. "Cache serving stale secret"

**Cause**: Cache not invalidated after rotation  
**Solution**: Disable cache or reduce TTL

```python
# Get fresh secret (bypass cache)
secret = manager.get("api_key", use_cache=False)

# Or reduce cache TTL
manager._cache_ttl = 60  # 1 minute
```

## Testing

Run secrets management tests:

```bash
pytest backend/tests/test_secrets_manager.py -v
```

**Coverage**:
- ✅ Local provider CRUD operations
- ✅ Secret versioning
- ✅ Automatic rotation
- ✅ Secrets manager caching
- ✅ Cache invalidation
- ✅ Batch operations
- ✅ Real-world scenarios (DB creds, API keys, TLS certs)

## Migration Guide

### From Environment Variables

```python
# Before: Reading from env vars
DB_PASSWORD = os.getenv("DB_PASSWORD")
API_KEY = os.getenv("NCBI_API_KEY")

# After: Reading from secrets manager
from app.security import get_secret

DB_PASSWORD = get_secret("db_password")
API_KEY = get_secret("api_key_ncbi")
```

### Migration Script

```python
import os
from app.security import get_secrets_manager

manager = get_secrets_manager()

# Migrate environment variables to secrets vault
env_to_secrets = {
    "DB_PASSWORD": "db_password",
    "DB_USER": "db_user",
    "NCBI_API_KEY": "api_key_ncbi",
    "OPENFDA_API_KEY": "api_key_openfda",
    "OPENAI_API_KEY": "api_key_openai",
    "SSL_CERT_FILE": "tls_cert_path",
    "SSL_KEY_FILE": "tls_key_path"
}

for env_var, secret_name in env_to_secrets.items():
    value = os.getenv(env_var)
    if value:
        manager.set(secret_name, value)
        print(f"Migrated {env_var} -> {secret_name}")

print("Migration complete! Remove env vars from .env file")
```

## References

- [HashiCorp Vault Documentation](https://www.vaultproject.io/docs)
- [AWS Secrets Manager Best Practices](https://docs.aws.amazon.com/secretsmanager/latest/userguide/best-practices.html)
- [Azure Key Vault Overview](https://docs.microsoft.com/en-us/azure/key-vault/)
- [OWASP Secrets Management Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Secrets_Management_Cheat_Sheet.html)
- [12-Factor App: Config](https://12factor.net/config)
