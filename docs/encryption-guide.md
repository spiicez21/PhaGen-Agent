# End-to-End Encryption Guide

## Overview

PhaGen implements enterprise-grade end-to-end encryption to protect sensitive pharmaceutical data at rest and in transit. This system ensures compliance with HIPAA, SOC2, ISO 27001, and GxP regulations.

## Features

### 🔐 Encryption at Rest
- **Algorithm**: AES-256-GCM (Galois/Counter Mode)
- **Authentication**: AEAD (Authenticated Encryption with Associated Data)
- **Key Management**: Envelope encryption with Data Encryption Keys (DEKs) wrapped by Key Encryption Keys (KEKs)
- **Scope**: Database fields, file storage, backups

### 🚀 Encryption in Transit
- **Protocol**: TLS 1.3 (latest version)
- **Ciphersuites**: Modern AEAD ciphers only
  - `TLS_AES_256_GCM_SHA384`
  - `TLS_CHACHA20_POLY1305_SHA256`
- **Certificate Management**: Automated generation and rotation

### 🔑 Key Management
- **Cloud KMS Support**:
  - AWS KMS
  - Azure Key Vault
  - Google Cloud KMS
- **Local Development**: File-based keys (dev/testing only)
- **Key Rotation**: Automatic rotation based on configurable intervals (default: 90 days)

## Architecture

### Envelope Encryption Pattern

```
┌─────────────────────────────────────────────────────────┐
│                                                         │
│  Plaintext Data                                         │
│       ↓                                                 │
│  Encrypt with DEK (AES-256-GCM)                         │
│       ↓                                                 │
│  Ciphertext + DEK                                       │
│       ↓                                                 │
│  Encrypt DEK with Master Key (KMS)                      │
│       ↓                                                 │
│  Ciphertext + Encrypted DEK                             │
│       ↓                                                 │
│  Store in Database/S3                                   │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

**Benefits**:
- Performance: Data encrypted with fast symmetric DEK
- Security: DEK protected by hardware-backed master key
- Rotation: Easy key rotation without re-encrypting all data
- Multi-tenant: Per-tenant DEKs for isolation

## Configuration

### Environment Variables

```bash
# Enable encryption
ENCRYPTION_ENABLED=true

# KMS provider: local, aws, azure, gcp
KMS_PROVIDER=aws

# Key rotation interval (days)
KEY_ROTATION_DAYS=90

# Encryption algorithm (default: AES-256-GCM)
ENCRYPTION_ALGORITHM=AES-256-GCM
```

### AWS KMS Configuration

```bash
# AWS credentials
AWS_REGION=us-east-1
AWS_ACCESS_KEY_ID=<your-key>
AWS_SECRET_ACCESS_KEY=<your-secret>

# Master key ID
AWS_KMS_KEY_ID=arn:aws:kms:us-east-1:123456789012:key/abcd1234-...
```

### Azure Key Vault Configuration

```bash
# Azure credentials
AZURE_TENANT_ID=<tenant-id>
AZURE_CLIENT_ID=<client-id>
AZURE_CLIENT_SECRET=<client-secret>

# Key vault URL
AZURE_KEY_VAULT_URL=https://phagen-vault.vault.azure.net/

# Key name
AZURE_KEY_NAME=phagen-master-key
```

### GCP Cloud KMS Configuration

```bash
# GCP credentials
GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json

# Project and key ring
GCP_PROJECT_ID=phagen-project
GCP_LOCATION=us-east1
GCP_KEY_RING=phagen-keyring
GCP_KEY_NAME=phagen-master-key
```

## Usage

### Database Field Encryption

```python
from app.security import EncryptionManager

manager = EncryptionManager()

# Encrypt sensitive field
plaintext = "Patient ID: 12345"
encrypted = manager.encrypt_field(plaintext, field_name="patient_id")

# Store encrypted in database
# {
#   "ciphertext": "...",
#   "nonce": "...",
#   "key_id": "dek-...",
#   "algorithm": "AES-256-GCM"
# }

# Decrypt when reading
decrypted = manager.decrypt_field(encrypted, field_name="patient_id")
# Returns: "Patient ID: 12345"
```

### File Encryption

```python
# Encrypt file
manager.encrypt_file(
    input_file="clinical_data.pdf",
    output_file="clinical_data.pdf.enc"
)

# Decrypt file
manager.decrypt_file(
    input_file="clinical_data.pdf.enc",
    output_file="clinical_data_decrypted.pdf"
)
```

### SQLAlchemy Model Integration

```python
from sqlalchemy import Column, String, JSON
from app.security import EncryptionManager

encryption = EncryptionManager()

class Patient(Base):
    __tablename__ = "patients"
    
    id = Column(String, primary_key=True)
    name_encrypted = Column(JSON)  # Store encrypted data
    
    @property
    def name(self):
        """Decrypt name on read"""
        return encryption.decrypt_field(self.name_encrypted, "name")
    
    @name.setter
    def name(self, value):
        """Encrypt name on write"""
        self.name_encrypted = encryption.encrypt_field(value, "name")
```

## TLS 1.3 Configuration

### Uvicorn Server

```python
from app.security import get_tls_config

tls_config = get_tls_config()

uvicorn.run(
    app,
    host="0.0.0.0",
    port=8000,
    ssl_version=tls_config["ssl_version"],
    ssl_cert_reqs=ssl.CERT_REQUIRED,
    ssl_ca_certs=tls_config["ssl_ca_file"],
    ssl_certfile=tls_config["ssl_cert_file"],
    ssl_keyfile=tls_config["ssl_key_file"],
    ssl_ciphers=":".join(tls_config["ssl_ciphers"])
)
```

### Generate Self-Signed Certificate (Development)

```python
from app.security import generate_self_signed_cert

generate_self_signed_cert(
    cert_file="certs/phagen.pem",
    key_file="certs/phagen.key",
    days_valid=365
)
```

**Note**: Use proper CA-signed certificates in production!

## Key Rotation

### Automatic Rotation

Keys are automatically rotated when they exceed the configured age:

```python
# Check if key needs rotation
if dek.needs_rotation(days=90):
    new_dek = manager.rotate_keys()
```

### Manual Rotation

```python
# Rotate all keys immediately
manager.rotate_keys()
```

### Zero-Downtime Rotation

The system maintains old keys temporarily to decrypt existing data while encrypting new data with rotated keys:

1. Generate new DEK
2. Encrypt new data with new DEK
3. Keep old DEK in cache for decryption
4. Gradually re-encrypt existing data in background
5. Retire old DEK after grace period

## Security Best Practices

### Production Checklist

- [ ] Use cloud KMS (AWS/Azure/GCP), never local keys
- [ ] Enable automatic key rotation (90 days or less)
- [ ] Use CA-signed TLS certificates
- [ ] Restrict KMS key access with IAM policies
- [ ] Enable CloudTrail/Azure Monitor for key usage auditing
- [ ] Store certificates in secure vault (HashiCorp Vault, AWS Secrets Manager)
- [ ] Use separate keys per environment (dev, staging, prod)
- [ ] Implement key backup and disaster recovery procedures
- [ ] Document key recovery procedures for auditors
- [ ] Test encryption/decryption performance under load

### Access Control

```json
// AWS KMS key policy example
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "AllowPhaGenBackend",
      "Effect": "Allow",
      "Principal": {
        "AWS": "arn:aws:iam::123456789012:role/PhaGenBackendRole"
      },
      "Action": [
        "kms:Encrypt",
        "kms:Decrypt",
        "kms:GenerateDataKey"
      ],
      "Resource": "*"
    }
  ]
}
```

## Compliance Mapping

| Requirement | Implementation | Evidence |
|------------|----------------|----------|
| **HIPAA §164.312(a)(2)(iv)** | AES-256-GCM encryption | `encryption.py` L40-70 |
| **HIPAA §164.312(e)(1)** | TLS 1.3 transport security | `encryption.py` L390-420 |
| **SOC2 CC6.7** | Cryptographic key management | `kms_providers.py` L1-250 |
| **ISO 27001 A.10.1.1** | Encryption policy | This document |
| **ISO 27001 A.10.1.2** | Key management | `encryption.py` L170-200 |
| **21 CFR Part 11** | Data integrity controls | Envelope encryption pattern |

## Troubleshooting

### Common Issues

#### 1. "Encryption key not found"

**Cause**: DEK expired or cache cleared  
**Solution**: Ensure key rotation grace period allows old keys to remain available

```python
# Extend key cache TTL
manager._dek_cache_ttl = 86400  # 24 hours
```

#### 2. "KMS authentication failed"

**Cause**: Invalid credentials or IAM permissions  
**Solution**: Verify AWS/Azure/GCP credentials and key access policies

```bash
# Test AWS KMS access
aws kms describe-key --key-id $AWS_KMS_KEY_ID

# Test Azure Key Vault access
az keyvault key show --vault-name phagen-vault --name phagen-master-key
```

#### 3. "Decryption failed: invalid tag"

**Cause**: Ciphertext tampered or wrong associated data  
**Solution**: Verify data integrity and field name matches

```python
# Ensure field_name matches encryption
encrypted = manager.encrypt_field(data, field_name="patient_id")
decrypted = manager.decrypt_field(encrypted, field_name="patient_id")
```

#### 4. "TLS handshake failed"

**Cause**: Certificate expired or cipher mismatch  
**Solution**: Regenerate certificate or update ciphersuites

```python
# Generate new certificate
generate_self_signed_cert("cert.pem", "key.pem", days_valid=365)
```

## Performance

### Benchmarks

- **Field Encryption**: ~0.5ms per operation (1KB data)
- **File Encryption**: ~50MB/s throughput
- **Key Generation**: ~10ms per DEK
- **KMS Operations**: ~50ms per encrypt/decrypt (network latency)

### Optimization

```python
# Batch encrypt multiple fields
fields = {"name": "John", "ssn": "123-45-6789", "diagnosis": "Condition XYZ"}
encrypted_fields = {
    name: manager.encrypt_field(value, name)
    for name, value in fields.items()
}

# Use connection pooling for KMS
import boto3
session = boto3.Session()
kms = session.client('kms', config=Config(max_pool_connections=50))
```

## Testing

Run encryption tests:

```bash
pytest backend/tests/test_encryption.py -v
```

**Coverage**:
- ✅ AES-256-GCM encryption/decryption
- ✅ Envelope encryption pattern
- ✅ Key rotation logic
- ✅ Database field encryption
- ✅ File encryption
- ✅ TLS 1.3 configuration
- ✅ Self-signed certificate generation
- ✅ Tampering detection

## References

- [NIST SP 800-57: Key Management](https://csrc.nist.gov/publications/detail/sp/800-57-part-1/rev-5/final)
- [AWS KMS Best Practices](https://docs.aws.amazon.com/kms/latest/developerguide/best-practices.html)
- [Azure Key Vault Security](https://docs.microsoft.com/en-us/azure/key-vault/general/security-features)
- [GCP Cloud KMS Overview](https://cloud.google.com/kms/docs/key-management-service)
- [OWASP Cryptographic Storage Cheat Sheet](https://cheatsheetseries.owasp.org/cheatsheets/Cryptographic_Storage_Cheat_Sheet.html)
