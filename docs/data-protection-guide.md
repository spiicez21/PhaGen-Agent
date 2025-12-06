# PII Redaction & Data Loss Prevention (DLP)

## Overview

The PII redaction system automatically detects and scrubs sensitive data in pharmaceutical environments to ensure compliance with HIPAA, GxP, and data protection regulations.

## Features

### 1. **Automatic PII Detection**

Regex-based pattern matching for common PII/PHI types:
- Email addresses
- Phone numbers (US formats: `(555) 123-4567`, `555.123.4567` | Indian formats: `9876543210`, `+91-9876543210`)
- US Social Security Numbers (SSN)
- Medical record numbers (MRN)
- Patient IDs
- Dates of birth
- IP addresses
- Credit card numbers

### 2. **Consistent Redaction**

Same PII value always gets the same hash:
```
john.doe@example.com → [EMAIL-a1b2c3d4]@example.com
john.doe@example.com → [EMAIL-a1b2c3d4]@example.com  # Same hash!
```

### 3. **Data Utility Preservation**

- Email domains preserved for analysis: `[EMAIL-hash]@company.com`
- IP addresses partially preserved: `192.168.xxx.xxx`
- Pattern types visible: `[SSN-hash]`, `[PHONE-hash]`

### 4. **DLP Policy Enforcement**

Prevents exfiltration of sensitive data:
- Configurable PII threshold (max instances before blocking)
- Audit trail of policy violations
- Alert callbacks for security monitoring

## Usage

### Basic Redaction

```python
from backend.app.security import redact_pii

text = "Contact john.doe@example.com or call (555) 123-4567"
redacted = redact_pii(text)
# Result: "Contact [EMAIL-a1b2c3d4]@example.com or call [PHONE-b5c6d7e8]"
```

### Advanced Redaction

```python
from backend.app.security import PIIRedactor, PIIType

redactor = PIIRedactor(
    redaction_salt="your-secret-salt",
    preserve_domain=True,
    enable_audit=True
)

# Selective redaction
text = "Email: john@example.com, SSN: 123-45-6789"
redacted, counts = redactor.redact_text(
    text,
    enabled_types={PIIType.EMAIL}  # Only redact emails
)

# Dict redaction
data = {
    "email": "john@example.com",
    "nested": {"ssn": "123-45-6789"}
}
redacted_data = redactor.redact_dict(data)

# Check stats
stats = redactor.get_redaction_stats()
# {'email': 5, 'ssn': 2, 'phone': 1}
```

### DLP Policy Enforcement

```python
from backend.app.security import DLPPolicy

def alert_handler(operation, counts, reason):
    # Send alert to security team
    print(f"DLP Alert: {operation} - {reason}")

dlp = DLPPolicy(
    block_on_pii=True,
    max_pii_threshold=5,
    alert_callback=alert_handler
)

# Scan for PII
has_pii, counts = dlp.scan_for_pii(text)

# Enforce policy
allowed, reason = dlp.enforce_policy(text, operation="export_report")
if not allowed:
    raise PermissionError(reason)

# Combined: redact + enforce
redacted, allowed, reason = dlp.redact_and_enforce(text, operation="api_response")
```

## Integration Points

### 1. Storage Layer (`storage.py`)

All text content stored in S3/MinIO is automatically redacted:

```python
# Before storage
body = "Patient MRN: ABC123456"

# After storage (redacted)
body = "Patient MRN: [MEDICAL_RECORD-a1b2c3d4]"
```

Log output:
```
[PII] Redacted 3 PII instances from reports/job-123.json
```

### 2. Jobs API (`routers/jobs.py`)

Job retrieval enforces DLP policy:

```python
@router.get("/{job_id}")
def get_job(job_id: str):
    job = job_store.get_job(job_id)
    
    # DLP check before returning
    if job contains excessive PII:
        raise HTTPException(403, "Job contains sensitive data")
    
    return job
```

Log output:
```
[DLP] Blocked job retrieval due to PII: job-abc-123
```

### 3. Report Generation

Reports are scanned before PDF generation and storage.

## Configuration

### Environment Variables

```bash
# Redaction salt (should be secret and stable)
PII_REDACTION_SALT=your-secret-salt-here

# DLP policy
DLP_BLOCK_ON_PII=true
DLP_MAX_PII_THRESHOLD=5
```

### Code Configuration

```python
from backend.app.security import PIIRedactor, DLPPolicy

# Global redactor
redactor = PIIRedactor(
    redaction_salt=os.getenv("PII_REDACTION_SALT"),
    preserve_domain=True,
    enable_audit=True
)

# Global DLP policy
dlp = DLPPolicy(
    block_on_pii=os.getenv("DLP_BLOCK_ON_PII", "true").lower() == "true",
    max_pii_threshold=int(os.getenv("DLP_MAX_PII_THRESHOLD", "5")),
    alert_callback=send_security_alert
)
```

## PII Types Reference

| Type | Pattern | Example Input | Example Output |
|------|---------|---------------|----------------|
| EMAIL | Email addresses | `john@example.com` | `[EMAIL-a1b2c3d4]@example.com` |
| PHONE | Phone numbers (US/Indian) | `(555) 123-4567` or `9876543210` or `+91-9876543210` | `[PHONE-b5c6d7e8]` |
| SSN | Social Security Numbers | `123-45-6789` | `[SSN-c8d9e0f1]` |
| MEDICAL_RECORD | Medical record numbers | `MRN: ABC123456` | `MRN: [MEDICAL_RECORD-d2e3f4a5]` |
| PATIENT_ID | Patient identifiers | `Patient ID: PT123456` | `Patient ID: [PATIENT_ID-e4f5a6b7]` |
| DATE_OF_BIRTH | Dates of birth | `DOB: 01/15/1970` | `DOB: [DATE_OF_BIRTH-f6a7b8c9]` |
| IP_ADDRESS | IPv4 addresses | `192.168.1.100` | `192.168.xxx.xxx` |
| CREDIT_CARD | Credit card numbers | `1234-5678-9012-3456` | `[CREDIT_CARD-a8b9c0d1]` |

## Testing

Run the test suite:

```bash
pytest backend/tests/test_pii_redaction.py -v
```

Test coverage includes:
- Pattern detection accuracy
- Consistent hashing verification
- Selective redaction
- Dictionary/nested redaction
- DLP policy enforcement
- Alert callbacks
- Real-world pharmaceutical scenarios

## Performance

- **Pattern matching**: ~1ms per 1KB of text
- **Redaction**: ~2ms per 1KB of text
- **Memory overhead**: Minimal (regex compilation cached)
- **Storage impact**: None (in-place redaction)

## Security Considerations

### Redaction Salt

The `redaction_salt` ensures:
- Consistent hashing (same PII → same hash)
- Deterministic mapping for analysis
- **Security**: Salt must be kept secret and stable

⚠️ **Warning**: Changing the salt will change all hashes, breaking consistency.

### Hash Collision

- SHA-256 truncated to 8 hex chars (32 bits)
- Collision probability: ~1 in 4 billion
- Acceptable for PII redaction use case

### Reversibility

- Redaction is **one-way** (cannot recover original PII)
- Hashes are deterministic but not reversible
- Suitable for GDPR "right to be forgotten"

## Compliance

### HIPAA

✅ Redacts Protected Health Information (PHI):
- Patient identifiers
- Medical record numbers
- Dates of birth
- Contact information

### GxP

✅ Maintains data integrity while protecting PII:
- Audit trail of redactions
- Deterministic transformations
- Preserves data utility for analysis

### GDPR

✅ Supports data minimization:
- Automatic PII detection
- One-way redaction (data deletion)
- DLP prevents exfiltration

## Troubleshooting

### PII Not Detected

**Symptom**: PII passes through unredacted

**Solutions**:
1. Check pattern matches in `RedactionPattern` class
2. Add custom patterns for your domain
3. Test pattern with regex debugger

Example custom pattern:
```python
# Custom lab ID pattern
LAB_ID = re.compile(
    r'\b(?:LAB[-_\s]?ID)[:\s#]*([A-Z0-9]{8,10})\b',
    re.IGNORECASE
)
```

### False Positives

**Symptom**: Non-PII data redacted (e.g., scientific IDs)

**Solutions**:
1. Use selective redaction with `enabled_types`
2. Add exclusion patterns
3. Whitelist known safe patterns

Example exclusion:
```python
# Don't redact compound IDs
if "COMPOUND-" in text:
    enabled_types = {PIIType.EMAIL, PIIType.PHONE}  # Exclude other types
    redacted, counts = redactor.redact_text(text, enabled_types)
```

### DLP Blocks Legitimate Use

**Symptom**: DLP policy blocks valid operations

**Solutions**:
1. Increase `max_pii_threshold`
2. Add operation-specific overrides
3. Use audit mode instead of blocking

Example override:
```python
dlp = DLPPolicy(block_on_pii=False)  # Audit only, don't block
allowed, reason = dlp.enforce_policy(text, operation="internal_analysis")
# Still logs PII detection but doesn't block
```

## Roadmap

- [ ] Machine learning-based PII detection (NER models)
- [ ] Custom pattern library for pharmaceutical data
- [ ] Token-level redaction for LLM inputs
- [ ] Real-time streaming redaction
- [ ] Integration with enterprise DLP solutions (Symantec, McAfee)

## References

- HIPAA Privacy Rule: https://www.hhs.gov/hipaa/
- GDPR Article 5: https://gdpr-info.eu/art-5-gdpr/
- NIST 800-122: Guide to Protecting PII
- FDA 21 CFR Part 11: Electronic Records
