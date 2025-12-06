# Container Image Signing and Supply-Chain Verification

This directory contains tools for signing container images and verifying their integrity in the supply chain.

## Overview

PhaGen uses **Cosign** (from Sigstore) to cryptographically sign container images and attach Software Bill of Materials (SBOM). This ensures:

1. **Image Integrity**: Images haven't been tampered with
2. **Provenance**: Images come from trusted build pipelines
3. **Supply Chain Security**: Dependencies are tracked via SBOM
4. **Vulnerability Management**: Automated scanning with Trivy

## Tools

### 1. `sign-containers.sh`

Bash script for local image signing and SBOM generation.

**Prerequisites:**
- [Cosign](https://github.com/sigstore/cosign) - Image signing
- [Syft](https://github.com/anchore/syft) - SBOM generation (optional)
- Docker - Container images must be built locally

**Installation:**
```bash
# Install Cosign
curl -O -L "https://github.com/sigstore/cosign/releases/latest/download/cosign-linux-amd64"
sudo mv cosign-linux-amd64 /usr/local/bin/cosign
sudo chmod +x /usr/local/bin/cosign

# Install Syft (optional)
curl -sSfL https://raw.githubusercontent.com/anchore/syft/main/install.sh | sh -s -- -b /usr/local/bin
```

**Usage:**
```bash
# Sign default PhaGen images
./infra/scripts/sign-containers.sh

# Sign specific images
./infra/scripts/sign-containers.sh phagen-backend:v1.0.0

# Sign multiple images
./infra/scripts/sign-containers.sh \
  phagen-backend:v1.0.0 \
  phagen-frontend:v1.0.0 \
  phagen-rdkit-service:v1.0.0
```

**What it does:**
1. Generates keypair (`cosign.key` + `cosign.pub`) if not present
2. Builds SBOM for each image using Syft
3. Signs image with Cosign
4. Attaches SBOM to image metadata
5. Generates provenance attestation (git commit, build date, etc.)
6. Verifies signature

**Key Management:**
- Private key: `cosign.key` (store securely, use password protection)
- Public key: `cosign.pub` (distribute to consumers)
- Use environment variable: `export COSIGN_PASSWORD=<password>`

### 2. `verify_images.py`

Python script for deployment-time image verification.

**Installation:**
```bash
# Install dependencies
pip install -r requirements.txt  # If needed

# Install Cosign (see above)
# Install Trivy (optional, for vulnerability scanning)
curl -sfL https://raw.githubusercontent.com/aquasecurity/trivy/main/contrib/install.sh | sh -s -- -b /usr/local/bin
```

**Usage:**
```bash
# Verify with public key
python infra/scripts/verify_images.py \
  --public-key cosign.pub \
  phagen-backend:latest

# Verify multiple images
python infra/scripts/verify_images.py \
  --public-key cosign.pub \
  phagen-backend:latest \
  phagen-frontend:latest \
  phagen-rdkit-service:latest

# Verify without vulnerability scanning
python infra/scripts/verify_images.py \
  --public-key cosign.pub \
  --no-vuln-check \
  phagen-backend:latest

# Set vulnerability severity threshold
python infra/scripts/verify_images.py \
  --public-key cosign.pub \
  --severity CRITICAL \
  phagen-backend:latest

# Verify keyless signatures (GitHub Actions)
python infra/scripts/verify_images.py \
  --organization your-org \
  ghcr.io/your-org/phagen-backend:latest
```

**What it does:**
1. Verifies image signature with Cosign
2. Checks SBOM attachment and parses it
3. Scans for vulnerabilities with Trivy (optional)
4. Generates verification report

**Exit codes:**
- `0`: All verifications passed
- `1`: At least one verification failed

### 3. `.github/workflows/sign-and-push.yml`

GitHub Actions workflow for automated signing in CI/CD.

**Features:**
- **Keyless signing**: Uses GitHub's OIDC tokens (no key management)
- **SBOM generation**: Automatic with Syft
- **Provenance attestation**: Build metadata attached to image
- **Multi-architecture**: Supports ARM64 + AMD64
- **Verification job**: Validates signatures after push
- **Security scanning**: Trivy integration with SARIF upload

**Workflow triggers:**
- Push to `main` or `release/*` branches
- Git tags matching `v*` (e.g., `v1.0.0`)
- Manual dispatch

**Configuration:**
```yaml
env:
  REGISTRY: ghcr.io  # GitHub Container Registry
  IMAGE_PREFIX: your-org/phagen  # Customize
```

**Matrix strategy:**
```yaml
matrix:
  service:
    - name: backend
      context: ./backend
      dockerfile: ./backend/Dockerfile
    - name: frontend
      context: ./frontend
      dockerfile: ./frontend/Dockerfile
    - name: rdkit-service
      context: ./backend
      dockerfile: ./backend/Dockerfile.rdkit
```

**Keyless verification:**
```bash
# Verify image signed by GitHub Actions
cosign verify \
  --certificate-identity-regexp="https://github.com/your-org/phagen" \
  --certificate-oidc-issuer="https://token.actions.githubusercontent.com" \
  ghcr.io/your-org/phagen-backend:latest
```

## Production Deployment Workflow

### Step 1: Build and Sign Images

**Option A: GitHub Actions (Recommended)**
```bash
# Push to main branch
git push origin main

# Or create release tag
git tag v1.0.0
git push origin v1.0.0
```

**Option B: Local Signing**
```bash
# Build images
docker build -t phagen-backend:v1.0.0 ./backend
docker build -t phagen-frontend:v1.0.0 ./frontend

# Sign images
./infra/scripts/sign-containers.sh \
  phagen-backend:v1.0.0 \
  phagen-frontend:v1.0.0

# Push to registry
docker push phagen-backend:v1.0.0
docker push phagen-frontend:v1.0.0
```

### Step 2: Distribute Public Key

**For key-based signing:**
```bash
# Copy public key to deployment environments
scp cosign.pub deployment-server:/etc/phagen/cosign.pub
```

**For keyless signing:**
No key distribution needed - verification uses GitHub's OIDC issuer.

### Step 3: Verify Before Deployment

**Add to deployment pipeline:**
```bash
# In Kubernetes admission controller
python infra/scripts/verify_images.py \
  --public-key /etc/phagen/cosign.pub \
  --fail-fast \
  phagen-backend:v1.0.0 \
  phagen-frontend:v1.0.0

# Exit code 0 = all passed, 1 = failed
if [ $? -eq 0 ]; then
  kubectl apply -f deployment.yaml
else
  echo "Image verification failed - aborting deployment"
  exit 1
fi
```

**In Docker Compose:**
```bash
# Pre-flight verification
python infra/scripts/verify_images.py \
  --public-key cosign.pub \
  phagen-backend:v1.0.0 \
  phagen-frontend:v1.0.0 \
  phagen-rdkit-service:v1.0.0

# Deploy only if verification passed
docker-compose up -d
```

### Step 4: Policy Enforcement

**Kubernetes Admission Controller** (recommended):
```yaml
# Use Sigstore Policy Controller
# https://docs.sigstore.dev/policy-controller/overview

apiVersion: policy.sigstore.dev/v1beta1
kind: ClusterImagePolicy
metadata:
  name: phagen-policy
spec:
  images:
  - glob: "ghcr.io/your-org/phagen-*"
  authorities:
  - keyless:
      url: https://fulcio.sigstore.dev
      identities:
      - issuer: https://token.actions.githubusercontent.com
        subject: https://github.com/your-org/phagen/.github/workflows/sign-and-push.yml@refs/heads/main
```

**Docker Content Trust** (alternative):
```bash
# Enable content trust
export DOCKER_CONTENT_TRUST=1

# Pull will fail if signature invalid
docker pull phagen-backend:v1.0.0
```

## SBOM Usage

### View SBOM

```bash
# Download SBOM from image
cosign download sbom phagen-backend:latest > sbom.json

# View with jq
cat sbom.json | jq '.packages[] | select(.name == "django")'

# Count packages
cat sbom.json | jq '.packages | length'
```

### Vulnerability Scanning

```bash
# Scan with Trivy
trivy image phagen-backend:latest

# Scan SBOM directly
trivy sbom sbom.json

# Generate report
trivy image --format json --output trivy-report.json phagen-backend:latest
```

### Dependency Tracking

```bash
# Extract all Python packages
cat sbom.json | jq -r '.packages[] | select(.type == "python") | "\(.name) \(.version)"'

# Check for specific CVE
trivy image --vuln-type os,library --severity HIGH,CRITICAL phagen-backend:latest | grep CVE-2024-12345
```

## Key Rotation

### Rotate Signing Keys

```bash
# Generate new keypair
cosign generate-key-pair

# This creates:
# - cosign.key (new private key)
# - cosign.pub (new public key)

# Backup old keys
mv cosign.key cosign.key.old
mv cosign.pub cosign.pub.old

# Sign new images with new key
./infra/scripts/sign-containers.sh phagen-backend:latest

# Update deployment environments with new public key
scp cosign.pub deployment-server:/etc/phagen/cosign.pub
```

### Key Expiration

Cosign keys don't expire by default. For production:

```bash
# Use short-lived keys (rotate monthly)
# Store old public keys for historical verification

# Example directory structure:
# /etc/phagen/keys/
#   2024-01/cosign.pub
#   2024-02/cosign.pub
#   2024-03/cosign.pub  # current
```

## Troubleshooting

### Signature Verification Fails

```bash
# Check if image has signature
cosign verify --key cosign.pub phagen-backend:latest

# Common issues:
# 1. Wrong public key
# 2. Image unsigned (no signature attached)
# 3. Image modified after signing
# 4. Registry doesn't support OCI artifacts
```

### SBOM Not Found

```bash
# Check if SBOM attached
cosign download sbom phagen-backend:latest

# If missing, regenerate:
syft phagen-backend:latest -o spdx-json > sbom.json
cosign attach sbom --sbom sbom.json phagen-backend:latest
```

### Trivy Scan Fails

```bash
# Update Trivy database
trivy image --download-db-only

# Clear cache
trivy image --clear-cache

# Retry scan
trivy image phagen-backend:latest
```

## Security Best Practices

1. **Private Key Protection**
   - Store `cosign.key` in secrets manager (AWS Secrets Manager, HashiCorp Vault)
   - Use password protection: `export COSIGN_PASSWORD=<strong-password>`
   - Rotate keys quarterly

2. **Keyless Signing (Production)**
   - Use GitHub Actions OIDC for automated builds
   - No key management required
   - Transparency log for audit trail

3. **Image Registry**
   - Use private registry for proprietary images
   - Enable vulnerability scanning at registry level
   - Configure immutable tags

4. **CI/CD Pipeline**
   - Sign images in CI, verify in CD
   - Block deployment of unsigned images
   - Archive SBOMs for compliance

5. **Monitoring**
   - Track signature verification failures
   - Alert on high/critical vulnerabilities
   - Audit SBOM for license compliance

## Compliance & Auditing

### Generate Compliance Report

```bash
# For all deployed images
./infra/scripts/verify_images.py \
  --public-key cosign.pub \
  phagen-backend:v1.0.0 \
  phagen-frontend:v1.0.0 \
  phagen-rdkit-service:v1.0.0 \
  > compliance-report.txt

# View SBOM for license audit
cosign download sbom phagen-backend:v1.0.0 | \
  jq '.packages[] | {name: .name, version: .version, license: .licenseConcluded}'
```

### Transparency Log

Keyless signatures are logged in Rekor (Sigstore's transparency log):

```bash
# View signature log entry
cosign verify phagen-backend:latest \
  --certificate-identity-regexp="https://github.com/your-org/phagen" \
  --certificate-oidc-issuer="https://token.actions.githubusercontent.com" \
  | jq '.[] | .optional.Bundle.Payload.body'
```

## References

- [Cosign Documentation](https://docs.sigstore.dev/cosign/overview/)
- [Syft SBOM Generator](https://github.com/anchore/syft)
- [Trivy Vulnerability Scanner](https://github.com/aquasecurity/trivy)
- [Sigstore Policy Controller](https://docs.sigstore.dev/policy-controller/overview/)
- [SLSA Framework](https://slsa.dev/)
