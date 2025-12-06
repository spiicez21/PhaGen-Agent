# PhaGen Agentic Architecture

## System Overview

```
┌────────────────────────────────────────────────────────────────────────┐
│                          PhaGen System Architecture                     │
├────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌──────────────┐         ┌──────────────┐         ┌───────────────┐  │
│  │   Frontend   │ ──────> │   Backend    │ ──────> │ Master Agent  │  │
│  │  (Next.js)   │         │  (FastAPI)   │         │               │  │
│  └──────────────┘         └──────────────┘         └───────────────┘  │
│         │                         │                    /  |  |  \      │
│         │                         │                   /   |  |   \     │
│         ▼                         ▼                  ▼    ▼  ▼    ▼    │
│  ┌──────────────┐         ┌──────────────┐    Clinical Patent Lit Mkt │
│  │ Evidence UI  │         │  Redis Cache │         Workers            │
│  │  Dashboard   │         │  (2hr/24hr)  │            │               │
│  └──────────────┘         └──────────────┘            ▼               │
│                                   │              ┌──────────────┐      │
│                                   │              │   Chroma DB  │      │
│                                   │              │   (Vectors)  │      │
│                                   ▼              └──────────────┘      │
│                           ┌──────────────┐                             │
│                           │  PostgreSQL  │                             │
│                           │  + pgBouncer │                             │
│                           └──────────────┘                             │
│                                                                         │
│  Supporting Services:                                                  │
│  ┌───────────┐  ┌────────────┐  ┌───────────┐  ┌──────────────┐     │
│  │  Celery   │  │   MinIO    │  │  Ollama   │  │ RDKit Service│     │
│  │  Workers  │  │  (S3 API)  │  │ (LLM API) │  │ (Rendering)  │     │
│  └───────────┘  └────────────┘  └───────────┘  └──────────────┘     │
└────────────────────────────────────────────────────────────────────────┘
```

## Core Components

### 1. **Frontend (Next.js)**
   - **Molecule Intake**: Form for SMILES/name input with validation
   - **Job Timeline**: Real-time progress visualization with worker status
   - **Evidence Dashboard**: Tabbed views for clinical/patent/literature/market evidence
   - **Report Viewer**: PDF export with structure rendering and citations
   - **Admin Console**: Crawler status, dataset management, index rebuilding
   - **Location**: `frontend/` directory
   - **Port**: 3000 (development), configurable for production

### 2. **Backend (FastAPI)**
   - **Job API**: 
     - `POST /api/jobs` - Create new analysis job (supports `?use_celery=true`)
     - `GET /api/jobs/{id}` - Retrieve job status and results
     - `GET /api/jobs` - List all jobs with filtering
   - **Feedback API**: Evidence rating for active learning (`/api/feedback`)
   - **Health Checks**: `/health`, `/ready`, `/metrics` endpoints
   - **Database**: PostgreSQL with SQLAlchemy ORM + pgBouncer pooling
   - **Cache**: Redis for retrieval (2hr TTL) and LLM responses (24hr TTL)
   - **Location**: `backend/` directory
   - **Port**: 8000 (development), configurable for production

### 3. **Agents Package**
   - **MasterAgent**: Orchestrates 4 worker agents, synthesizes innovation story
   - **Worker Agents**:
     - `ClinicalWorker`: Trials, phases, endpoints, populations from registries
     - `PatentWorker`: Assignees, claims, priority dates, blocking analysis
     - `LiteratureWorker`: Mechanism-of-action, DOI citations, evidence quality
     - `MarketWorker`: TAM, incidence, competition with market API integration
   - **Shared Components**:
     - Retriever with Chroma/FAISS backend
     - LLM client (Ollama Gemma2:2B default)
     - Synonym expansion for query enrichment
     - Source ranking scorecard (clinical > regulatory > literature > patent)
   - **Location**: `agents/` directory

### 4. **Crawler (Crawlee + TypeScript)**
   - API-first sources with robots-aware fallback
   - Normalization pipeline: strips boilerplate, redacts PII, chunks with metadata
   - Domain-level crawl budgets (default 2 pages/host)
   - Evidence deduplication (clinical wins over literature)
   - Outputs to `crawler/storage/datasets/` for indexing
   - **Continuous Crawler**: Change detection with SHA-256 hashing (`crawler/continuous_crawler.py`)
   - **Location**: `crawler/` directory
   - **Usage**: `npm run crawl`

### 5. **Vector Database (ChromaDB)**
   - Persistent storage: `indexes/chroma/`
   - Embeddings: sentence-transformers/all-MiniLM-L6-v2 (default)
   - Custom embedding training pipeline available (`backend/app/ml/embedding_trainer.py`)
   - Embedding cache: `.embedding_cache.json` for incremental updates
   - Index snapshots: Daily with manifests under `indexes/chroma_snapshots/`
   - **Location**: `indexes/` directory
   - **Build**: `python indexes/build_index.py`

### 6. **Infrastructure Services**

#### Docker Compose Stack (`infra/docker-compose.yml`)
   - **PostgreSQL**: Primary database (Supabase or self-hosted)
   - **pgBouncer**: Connection pooler (port 6432, 25 default/100 max connections)
   - **Redis**: Cache + Celery broker (port 6379)
   - **MinIO**: S3-compatible object storage (port 9000/9001)
   - **Ollama**: Local LLM hosting (port 11434)
   - **RDKit Service**: SMILES rendering (port 8001)

#### Kubernetes Deployment (`k8s/`)
   - **API Deployment**: 3 replicas with HPA (2-10 range)
   - **Celery Workers**: 2 replicas with HPA (1-8 range)
   - **Redis**: Persistent volume for cache/queue
   - **Ingress**: NGINX with TLS support
   - **ConfigMaps/Secrets**: Environment configuration

### 7. **Advanced Features**

#### Patent Claims Analysis (`agents/workers/patent_claims.py`)
   - Claim extraction (independent/dependent parsing)
   - Semantic blocking risk scoring
   - Freedom-to-operate (FTO) reports with ClaimMap visualization

#### Market APIs (`agents/workers/market_apis.py`)
   - Multi-provider aggregator (IQVIA, Evaluate Pharma)
   - Automatic fallback to synthetic data
   - Confidence scoring for data quality

#### Regulatory Engine (`agents/workers/regulatory_engine.py`)
   - Multi-region rules (US FDA, EU EMA, India CDSCO)
   - Accelerated pathway detection
   - Timeline/probability modeling

### 8. **Enterprise Security**

#### HSM Integration (`backend/app/security/hsm_manager.py`)
   - Hardware-backed key storage (SoftHSM/AWS CloudHSM/Azure Key Vault)
   - PKCS#11 interface with AES-GCM encryption
   - Database encryption key management

#### Container Signing (`infra/scripts/`)
   - Cosign signatures with keyless GitHub OIDC
   - SBOM generation via Syft
   - Trivy vulnerability scanning
   - Verification scripts for deployment

#### Air-Gapped Mode (`infra/airgap/`)
   - Bundle creator for offline deployment
   - Packages models, indexes, dependencies, images
   - Zero external network access
   - Installation automation

## Data Flow

### Standard Job Execution

```
1. User Input
   └─> Frontend form (molecule name + SMILES)
       └─> POST /api/jobs
           └─> Backend creates job record in PostgreSQL

2. Job Processing (Two Modes)
   
   A) Synchronous (default):
      └─> Background task spawns MasterAgent
          └─> Synonym expansion enriches query
          └─> 4 workers execute in parallel:
              ├─> ClinicalWorker → Chroma retrieval → LLM synthesis
              ├─> PatentWorker → Chroma retrieval → LLM synthesis  
              ├─> LiteratureWorker → Chroma retrieval → LLM synthesis
              └─> MarketWorker → Chroma/API → LLM synthesis
          └─> MasterAgent synthesizes innovation story
          └─> Results stored in PostgreSQL + Redis cache

   B) Distributed (Celery):
      └─> Job enqueued to Redis
          └─> Celery worker picks up task
              └─> (Same processing as synchronous)
              └─> Results stored via Celery backend

3. Result Retrieval
   └─> Frontend polls GET /api/jobs/{id}
       └─> Backend checks Redis cache (hit: 2hr TTL)
           └─> Cache miss: Query PostgreSQL
               └─> Return job status + evidence payload

4. User Interaction
   └─> Evidence dashboard renders tabbed views
       └─> User upvotes/downvotes evidence → POST /api/feedback
           └─> Feedback stored for reranker training
       └─> Export to PDF → Backend generates report
           └─> MinIO stores PDF artifact
           └─> Frontend downloads signed URL
```

### Data Processing Pipeline

```
Data Ingestion:
1. Crawler (Crawlee) fetches documents
   └─> Normalization: Strip boilerplate, redact PII, chunk
       └─> Output: crawler/storage/datasets/default/*.json

2. Index Builder (indexes/build_index.py)
   └─> Load documents from datasets
       └─> Embedding generation (cached in .embedding_cache.json)
           └─> SMILES normalization (RDKit canonicalization)
               └─> Deduplication (clinical > literature priority)
                   └─> ChromaDB persistence (indexes/chroma/)
                       └─> Snapshot creation (daily manifests)

Retrieval:
1. Worker query → Synonym expansion
   └─> ChromaDB semantic search (top-k passages)
       └─> Source ranking (clinical > regulatory > literature > patent)
           └─> Context budgeting per worker
               └─> Return ranked evidence chunks

LLM Synthesis:
1. Worker receives ranked passages
   └─> Prompt template with evidence context
       └─> Ollama LLM call (Gemma2:2B default)
           └─> Redis cache (24hr TTL) for LLM responses
               └─> Structured JSON output with citations
```

### Continuous Updates

```
Continuous Crawler:
1. Periodic scan of source documents
   └─> SHA-256 content hashing
       └─> Change detection (new/modified/deleted)
           └─> Freshness scoring (30-day decay)
               └─> Trigger incremental reindex for changed docs

Active Learning:
1. User feedback collection
   └─> Evidence upvotes/downvotes stored in evidence_feedback table
       └─> Aggregate statistics via GET /api/feedback/stats
           └─> Reranker training (ml_cli.py reranker)
               └─> Updated weights improve future retrieval
```

## Technology Stack

### Frontend
- **Framework**: Next.js 14 (React 18)
- **Styling**: Tailwind CSS
- **State Management**: React hooks + SWR for data fetching
- **Visualization**: Recharts for timeline/charts
- **PDF Generation**: React-PDF or server-side generation

### Backend
- **Framework**: FastAPI 0.109+
- **Database**: PostgreSQL 15+ with SQLAlchemy 2.0
- **Cache**: Redis 7.0+ with redis-py
- **Task Queue**: Celery 5.3+ with Redis broker
- **Object Storage**: MinIO (S3-compatible)
- **Connection Pool**: pgBouncer (transaction mode)

### AI/ML
- **LLM**: Ollama (Gemma2:2B, smollm:360m, or custom)
- **Embeddings**: sentence-transformers (all-MiniLM-L6-v2)
- **Vector DB**: ChromaDB 0.4+
- **Chemistry**: RDKit 2023.9+ for SMILES rendering
- **Custom Training**: PyTorch 2.0+ for fine-tuning

### Infrastructure
- **Containerization**: Docker + Docker Compose
- **Orchestration**: Kubernetes 1.24+ (AWS EKS, Azure AKS, GCP GKE)
- **CI/CD**: GitHub Actions
- **Monitoring**: Prometheus metrics endpoint
- **Security**: Cosign, Trivy, HSM providers

### Crawler
- **Framework**: Crawlee (TypeScript/Node.js)
- **Storage**: JSONL datasets with metadata

## Deployment Modes

### 1. Local Development
```bash
docker-compose up
# Services on localhost:
# - Frontend: 3000
# - Backend: 8000
# - Redis: 6379
# - MinIO: 9000/9001
# - Ollama: 11434
# - RDKit: 8001
```

### 2. Production (Kubernetes)
```bash
kubectl apply -f k8s/
# Features:
# - Horizontal autoscaling (HPA)
# - Load balancing (Ingress)
# - Health checks (liveness/readiness)
# - Secrets management
# - Persistent volumes
```

### 3. Air-Gapped (Offline)
```bash
# Create bundle on connected system:
python infra/airgap/create_bundle.py --output bundle

# Transfer to air-gapped system
# Extract and install:
sudo ./install.sh

# Set offline mode:
export PHAGEN_OFFLINE_MODE=true
docker-compose up
```

## Security Architecture

### Defense in Depth

1. **Network Layer**
   - Default-deny egress controls
   - VPC/subnet isolation in cloud deployments
   - TLS 1.3 for all external communication

2. **Application Layer**
   - Input validation with Pydantic models
   - SQL injection prevention via SQLAlchemy ORM
   - CORS configuration for API access
   - Rate limiting on sensitive endpoints

3. **Data Layer**
   - Database encryption at rest (AES-256)
   - HSM for key management (FIPS 140-2 ready)
   - PII redaction in crawled content
   - Audit logging for all data access

4. **Container Layer**
   - Signed images with Cosign
   - SBOM for dependency tracking
   - Vulnerability scanning with Trivy
   - Minimal base images (distroless where possible)

5. **Compliance**
   - Zero Data Retention (ZDR) mode available
   - Air-gapped deployment for regulated environments
   - Audit trail for all operations
   - RBAC ready (future enhancement)

## Performance Characteristics

### Typical Job Execution
- **Cold start**: 15-30 seconds (first job, model loading)
- **Warm execution**: 5-10 seconds (cached models)
- **Retrieval**: 1-2 seconds per worker (Chroma query)
- **LLM synthesis**: 2-5 seconds per worker (Ollama)
- **Total**: 8-15 seconds for complete analysis

### Scalability
- **Horizontal**: 2-10 API replicas, 1-8 Celery workers (HPA)
- **Throughput**: 10-50 jobs/minute (depending on cluster size)
- **Storage**: 10GB baseline + 1GB per 10,000 documents indexed
- **Memory**: 4-8GB per backend replica (LLM models in RAM)

### Optimization
- **Redis caching**: 80-90% cache hit rate for repeat queries
- **Embedding cache**: Speeds up index rebuilds by 5-10x
- **pgBouncer pooling**: Reduces DB connection overhead by 60-70%
- **Celery distribution**: Linear scaling up to 8 workers

## Monitoring & Observability

### Health Checks
- `/health` - Basic liveness (returns 200 if service running)
- `/ready` - Readiness check (validates DB/Redis connectivity)
- `/metrics` - Prometheus-style metrics (request counts, latencies)

### Key Metrics
- Job success/failure rate
- Worker execution time per agent
- Retrieval precision and coverage
- Cache hit rates (Redis)
- Database connection pool utilization
- LLM token usage

### Logging
- Structured JSON logs for all components
- SLA metrics per worker (latency, tokens, failures)
- Citation gap detection for hallucination monitoring
- API call budget tracking (NCBI, OpenFDA rate limits)

## Development Workflow

1. **Setup**: Clone repo, run `docker-compose up`
2. **Frontend**: `cd frontend && npm run dev`
3. **Backend**: `cd backend && uvicorn app.main:app --reload`
4. **Crawler**: `cd crawler && npm run crawl`
5. **Indexing**: `cd indexes && python build_index.py`
6. **Testing**: `pytest backend/tests/`
7. **CI**: GitHub Actions validates on every PR

See `02-GETTING-STARTED.md` for detailed setup instructions.
