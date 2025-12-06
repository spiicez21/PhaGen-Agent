# Air-Gapped Deployment Guide

This guide covers deploying PhaGen in secure, offline environments (air-gapped) commonly found in pharmaceutical companies and regulated industries.

## Overview

Air-gapped deployment enables PhaGen to run completely offline with:
- No external API calls (PubMed, FDA, ClinicalTrials.gov disabled)
- Local ML models and vector indexes
- All dependencies bundled
- Container images pre-loaded
- No telemetry or automatic updates

## Use Cases

- **Secure pharma R&D labs**: No internet access due to IP protection
- **Clinical sites**: HIPAA/GDPR compliance requirements
- **Government facilities**: Air-gapped for security clearance
- **Remote locations**: Limited or unreliable network connectivity

## Architecture

```
Air-Gapped Environment
┌─────────────────────────────────────────┐
│                                         │
│  ┌──────────────┐   ┌───────────────┐  │
│  │   Frontend   │──▶│   Backend     │  │
│  │  Container   │   │   Container   │  │
│  └──────────────┘   └───────────────┘  │
│                            │            │
│                            ▼            │
│                     ┌──────────────┐   │
│                     │  PostgreSQL  │   │
│                     └──────────────┘   │
│                            │            │
│         Local Resources    ▼            │
│        ┌────────────────────────────┐  │
│        │ /opt/phagen/               │  │
│        │  ├── models/               │  │
│        │  │   ├── all-MiniLM-L6-v2 │  │
│        │  │   └── reranker/        │  │
│        │  ├── indexes/              │  │
│        │  │   └── chroma/          │  │
│        │  └── config/               │  │
│        └────────────────────────────┘  │
│                                         │
│  No External Network Access 🔒          │
└─────────────────────────────────────────┘
```

## Bundle Creation (Connected System)

### Step 1: Install Bundle Creation Tool

```bash
# Clone repository
git clone https://github.com/your-org/phagen.git
cd phagen

# Install dependencies
pip install -r requirements.txt
pip install sentence-transformers
```

### Step 2: Build Container Images

```bash
# Build all images
docker build -t phagen-backend:v1.0.0 ./backend
docker build -t phagen-frontend:v1.0.0 ./frontend
docker build -t phagen-rdkit-service:v1.0.0 -f ./backend/Dockerfile.rdkit ./backend
```

### Step 3: Create Bundle

```bash
# Full bundle with all components
python infra/airgap/create_bundle.py \
  --output airgap-bundle \
  --models sentence-transformers/all-MiniLM-L6-v2 \
  --indexes indexes/chroma \
  --requirements backend/requirements.txt requirements.txt \
  --images phagen-backend:v1.0.0 phagen-frontend:v1.0.0 phagen-rdkit-service:v1.0.0 \
  --config docker-compose.yml .env.example

# This will create:
# airgap-bundle/
#   ├── models/           # ML models
#   ├── indexes/          # Vector indexes
#   ├── dependencies/     # Python packages
#   ├── images/           # Docker images (tar files)
#   ├── config/           # Configuration files
#   ├── install.sh        # Installation script
#   ├── manifest.json     # Bundle manifest
#   └── README.md         # Installation instructions
```

### Step 4: Package Bundle

```bash
# Create compressed archive
tar -czf phagen-airgap-v1.0.0.tar.gz airgap-bundle/

# Calculate checksum
sha256sum phagen-airgap-v1.0.0.tar.gz > phagen-airgap-v1.0.0.sha256

# Bundle is now ready for transfer
ls -lh phagen-airgap-v1.0.0.*
```

## Bundle Transfer

Transfer the bundle to the air-gapped environment using approved secure channels:

### Option A: USB Drive
```bash
# Copy to USB drive
cp phagen-airgap-v1.0.0.tar.gz /media/usb-drive/
cp phagen-airgap-v1.0.0.sha256 /media/usb-drive/
```

### Option B: Secure File Transfer (if available)
```bash
# Via secure one-way transfer system
scp phagen-airgap-v1.0.0.tar.gz secure-gateway:/incoming/
```

### Option C: Physical Media (CD/DVD)
```bash
# Burn to DVD
brasero phagen-airgap-v1.0.0.tar.gz
```

## Installation (Air-Gapped System)

### Step 1: Verify Checksum

```bash
# Verify bundle integrity
sha256sum -c phagen-airgap-v1.0.0.sha256
# Should output: phagen-airgap-v1.0.0.tar.gz: OK
```

### Step 2: Extract Bundle

```bash
# Extract bundle
tar -xzf phagen-airgap-v1.0.0.tar.gz
cd airgap-bundle

# Verify manifest
cat manifest.json | jq '.bundle_version, .created_at'
```

### Step 3: Run Installation

```bash
# Run installation script (requires root)
sudo ./install.sh

# Script will:
# 1. Install Python dependencies from local wheels
# 2. Load Docker images from tar files
# 3. Copy models and indexes to /opt/phagen
# 4. Set offline mode environment variables
```

### Step 4: Configure Environment

```bash
# Navigate to installation directory
cd /opt/phagen/config

# Copy and edit configuration
cp docker-compose.yml.example docker-compose.yml
cp .env.example .env

# Edit .env for air-gapped mode
nano .env
```

Edit `.env`:
```bash
# Offline mode (REQUIRED)
PHAGEN_OFFLINE_MODE=true

# Local resource paths
PHAGEN_MODELS_DIR=/opt/phagen/models
PHAGEN_INDEXES_DIR=/opt/phagen/indexes

# Disable external APIs
PUBMED_API_ENABLED=false
CLINICALTRIALS_API_ENABLED=false
FDA_API_ENABLED=false

# Database configuration
POSTGRES_USER=phagen
POSTGRES_PASSWORD=<generate-secure-password>
POSTGRES_DB=phagen

# Security settings
SECRET_KEY=<generate-secure-key>
```

### Step 5: Start PhaGen

```bash
# Start all services
docker-compose up -d

# Check status
docker-compose ps

# View logs
docker-compose logs -f

# Wait for services to be healthy
docker-compose ps | grep "healthy"
```

### Step 6: Verify Installation

```bash
# Check backend health
curl http://localhost:8000/health
# Expected: {"status": "healthy", "offline_mode": true}

# Check frontend
curl http://localhost:3000
# Expected: HTML response

# Verify offline mode
curl http://localhost:8000/status | jq '.offline_mode'
# Expected: true
```

## Configuration

### Offline Mode Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PHAGEN_OFFLINE_MODE` | `false` | Enable offline mode |
| `PHAGEN_MODELS_DIR` | `/opt/phagen/models` | Local models directory |
| `PHAGEN_INDEXES_DIR` | `/opt/phagen/indexes` | Local indexes directory |
| `PUBMED_API_ENABLED` | `true` | Disable in offline mode |
| `CLINICALTRIALS_API_ENABLED` | `true` | Disable in offline mode |
| `FDA_API_ENABLED` | `true` | Disable in offline mode |

### Offline Feature Availability

| Feature | Online | Offline |
|---------|--------|---------|
| Patent search | ✅ | ✅ (local data only) |
| Clinical trials | ✅ | ✅ (bundled data) |
| Literature search | ✅ | ✅ (bundled data) |
| Market analysis | ✅ | ✅ (synthetic estimates) |
| Regulatory analysis | ✅ | ✅ (rule-based) |
| Vector similarity search | ✅ | ✅ (local indexes) |
| ML reranking | ✅ | ✅ (local models) |
| External API data | ✅ | ❌ (disabled) |
| Automatic updates | ✅ | ❌ (disabled) |
| Telemetry | ✅ | ❌ (disabled) |

## Updating Bundle

### Create Update Bundle

```bash
# On connected system, create incremental update
python infra/airgap/create_bundle.py \
  --output airgap-bundle-v1.1.0 \
  --models sentence-transformers/new-model \
  --indexes indexes/chroma \
  --images phagen-backend:v1.1.0

# Package update
tar -czf phagen-airgap-update-v1.1.0.tar.gz airgap-bundle-v1.1.0/
```

### Apply Update

```bash
# On air-gapped system
cd /opt/phagen

# Stop services
docker-compose down

# Extract update bundle
tar -xzf phagen-airgap-update-v1.1.0.tar.gz
cd airgap-bundle-v1.1.0

# Run installation script
sudo ./install.sh

# Restart services
cd /opt/phagen/config
docker-compose up -d
```

## Troubleshooting

### Issue: Services not starting

```bash
# Check logs
docker-compose logs

# Common causes:
# 1. Offline mode not enabled
grep PHAGEN_OFFLINE_MODE .env

# 2. Models directory not found
ls -la /opt/phagen/models

# 3. Database connection issue
docker-compose exec backend python -c "from app.main import test_db; test_db()"
```

### Issue: Models not loading

```bash
# Verify models exist
ls -la /opt/phagen/models

# Check permissions
sudo chown -R 1000:1000 /opt/phagen/models

# Test model loading
docker-compose exec backend python -c "
from sentence_transformers import SentenceTransformer
model = SentenceTransformer('/opt/phagen/models/sentence-transformers_all-MiniLM-L6-v2')
print('Model loaded successfully')
"
```

### Issue: Indexes not found

```bash
# Verify indexes exist
ls -la /opt/phagen/indexes/chroma

# Check ChromaDB can access
docker-compose exec backend python -c "
import chromadb
client = chromadb.PersistentClient(path='/opt/phagen/indexes/chroma')
print(client.list_collections())
"
```

### Issue: Feature requires online access

```bash
# Check if offline mode is enabled
curl http://localhost:8000/status | jq '.offline_mode'

# Enable offline mode if not set
echo "PHAGEN_OFFLINE_MODE=true" >> .env
docker-compose restart
```

## Data Updates

### Refresh Local Data (Periodic)

To update local data (clinical trials, patents, etc.) in offline mode:

1. **On connected system**: Export fresh data
   ```bash
   # Export latest data
   python scripts/export_data.py \
     --output data-export-$(date +%Y%m%d).tar.gz \
     --clinical-trials \
     --patents \
     --literature
   ```

2. **Transfer to air-gapped system** via secure channel

3. **Import data**:
   ```bash
   # On air-gapped system
   python scripts/import_data.py \
     --input data-export-20240315.tar.gz \
     --replace
   ```

### Rebuild Indexes

```bash
# On air-gapped system
docker-compose exec backend python -c "
from indexes.build_index import rebuild_all_indexes
rebuild_all_indexes()
"
```

## Security Considerations

### 1. Container Image Verification

```bash
# Verify image signatures before loading
python infra/scripts/verify_images.py \
  --public-key /secure/cosign.pub \
  phagen-backend:v1.0.0

# Load only if verification passes
docker load -i images/phagen-backend_v1.0.0.tar
```

### 2. Bundle Integrity

```bash
# Always verify checksums
sha256sum -c phagen-airgap-v1.0.0.sha256

# Verify manifest checksum
python -c "
import json
with open('manifest.json') as f:
    manifest = json.load(f)
print(f\"Bundle checksum: {manifest['checksum']}\")
"
```

### 3. Secure Credentials

```bash
# Generate secure passwords
openssl rand -base64 32  # For SECRET_KEY
openssl rand -base64 24  # For POSTGRES_PASSWORD

# Store in .env with restricted permissions
chmod 600 .env
```

### 4. Access Control

```bash
# Restrict access to installation directory
sudo chown -R root:phagen /opt/phagen
sudo chmod -R 750 /opt/phagen

# Create service user
sudo useradd -r -s /bin/false phagen
sudo usermod -aG docker phagen
```

## Compliance & Auditing

### Logging

```bash
# All operations are logged
docker-compose logs > /var/log/phagen/app.log

# Audit log for offline mode operations
tail -f /var/log/phagen/audit.log
```

### Manifest Tracking

```bash
# Record deployed bundle version
cat /opt/phagen/manifest.json | jq '{version, created_at, checksum}' \
  > /var/log/phagen/deployment-$(date +%Y%m%d).json
```

### Compliance Reports

```bash
# Generate SBOM for compliance
cosign download sbom phagen-backend:v1.0.0 > sbom-backend.json

# Check for vulnerabilities (use before deployment)
trivy sbom sbom-backend.json
```

## Performance Optimization

### Resource Allocation

```yaml
# docker-compose.yml
services:
  backend:
    deploy:
      resources:
        limits:
          cpus: '4'
          memory: 8G
        reservations:
          cpus: '2'
          memory: 4G
```

### Database Tuning

```bash
# Increase PostgreSQL shared buffers for offline mode
# In docker-compose.yml:
POSTGRES_SHARED_BUFFERS=2GB
POSTGRES_WORK_MEM=64MB
POSTGRES_MAINTENANCE_WORK_MEM=512MB
```

### Model Caching

```python
# Preload models on startup for faster queries
# In backend/app/main.py:
from app.offline_config import offline_config

@app.on_event("startup")
async def load_models():
    if offline_config.enabled:
        # Preload embedding model
        model_path = offline_config.get_model_path("sentence-transformers/all-MiniLM-L6-v2")
        global embedding_model
        embedding_model = SentenceTransformer(str(model_path))
```

## Support

For issues specific to air-gapped deployment:

1. **Collect logs**:
   ```bash
   docker-compose logs > logs.txt
   tar -czf phagen-logs.tar.gz logs.txt /var/log/phagen/
   ```

2. **Generate diagnostics**:
   ```bash
   python infra/airgap/diagnostics.py > diagnostics.txt
   ```

3. **Transfer to connected system** for analysis via secure channel

## Appendix

### Bundle Contents Example

```
airgap-bundle/
├── models/
│   ├── sentence-transformers_all-MiniLM-L6-v2/
│   │   ├── config.json
│   │   ├── pytorch_model.bin
│   │   └── tokenizer.json
│   └── reranker/
│       └── model.safetensors
├── indexes/
│   └── chroma/
│       ├── chroma.sqlite3
│       └── [collection data]
├── dependencies/
│   ├── requirements.txt
│   ├── packages/  # Source distributions
│   └── wheels/    # Binary wheels
├── images/
│   ├── phagen-backend_v1.0.0.tar (2.1 GB)
│   ├── phagen-frontend_v1.0.0.tar (1.5 GB)
│   └── phagen-rdkit-service_v1.0.0.tar (3.2 GB)
├── config/
│   ├── docker-compose.yml.example
│   └── .env.example
├── install.sh
├── manifest.json
└── README.md

Total size: ~8 GB compressed, ~15 GB uncompressed
```

### Manifest Schema

```json
{
  "bundle_version": "1.0.0",
  "created_at": "2024-03-15T10:30:00Z",
  "python_version": "3.10.12",
  "models": [
    {
      "name": "sentence-transformers/all-MiniLM-L6-v2",
      "type": "embedding",
      "path": "models/sentence-transformers_all-MiniLM-L6-v2",
      "size_mb": 90.5
    }
  ],
  "indexes": [
    {
      "name": "chroma",
      "path": "indexes/chroma",
      "size_mb": 1250.3
    }
  ],
  "dependencies": {
    "fastapi": "fastapi==0.109.0",
    "sentence-transformers": "sentence-transformers>=2.3.0"
  },
  "container_images": [
    {
      "name": "phagen-backend:v1.0.0",
      "filename": "phagen-backend_v1.0.0.tar",
      "size_mb": 2100.0
    }
  ],
  "checksum": "sha256:abcdef123456..."
}
```
