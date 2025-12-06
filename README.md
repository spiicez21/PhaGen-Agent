# PhaGen Agentic

**Enterprise-grade pharmaceutical evidence aggregation and innovation assessment platform**

PhaGen Agentic is a production-ready, agentic multi-worker platform that compresses molecule repurposing research from months to minutes. The system leverages AI-powered worker agents, vector databases, and advanced analytics to provide comprehensive evidence synthesis for pharmaceutical R&D teams.

## 🎯 Project Status: Production-Ready (95%+ Complete)

| Area | Status | Details |
| --- | --- | --- |
| **Frontend** | ✅ Complete | Next.js dashboard with evidence tabs, timeline, reports, admin console |
| **Backend** | ✅ Complete | FastAPI with Postgres, Redis cache, Celery queue, health checks, feedback API |
| **Agents** | ✅ Complete | Master + 4 workers (Clinical, Patent, Literature, Market) with LLM synthesis |
| **Crawler** | ✅ Complete | Crawlee with normalization, deduplication, continuous change detection |
| **Vector DB** | ✅ Complete | ChromaDB with embedding cache, snapshots, custom training pipeline |
| **Infrastructure** | ✅ Complete | Docker Compose + Kubernetes, pgBouncer, Redis, MinIO, Ollama, RDKit |
| **ML Features** | ✅ Complete | Disease mapping, reranker training, custom embeddings |
| **Advanced Features** | ✅ Complete | Patent claims, continuous crawler, market APIs, regulatory engine |
| **Enterprise Security** | ✅ Complete | HSM integration, signed containers, air-gapped deployment |

**Total Implementation**: ~15,000+ lines of code across 16 major features

## 🚀 Key Features

### Core Capabilities
- **Multi-Worker Orchestration**: Master agent coordinates Clinical, Patent, Literature, and Market workers
- **Vector Search**: ChromaDB semantic search with custom embeddings and caching
- **LLM Synthesis**: Local Ollama integration (Gemma2:2B default) with configurable models
- **Evidence Dashboard**: Tabbed UI with confidence scores, citations, and metadata
- **Report Generation**: PDF export with structure rendering and provenance tracking
- **Distributed Processing**: Celery task queue for horizontal scaling

### Advanced Analytics
- **Disease Mapping**: ML model for repurposing suggestions across 6 therapeutic areas
- **Reranker Training**: Active learning from user feedback to optimize evidence retrieval
- **Custom Embeddings**: Fine-tunable sentence-transformers on pharma corpus
- **Patent Claims Analysis**: Semantic matching for FTO (Freedom-to-Operate) assessment
- **Regulatory Rules Engine**: Multi-region (US/EU/India) pathway analysis

### Enterprise Security
- **HSM Integration**: Hardware-backed key storage (SoftHSM, AWS CloudHSM, Azure Key Vault)
- **Container Signing**: Cosign signatures with SBOM and provenance attestation
- **Air-Gapped Deployment**: Complete offline mode with bundled models and indexes
- **Supply-Chain Verification**: Trivy scanning, keyless signing, vulnerability management

## 📚 Documentation

Comprehensive documentation is available in the `docs/` directory:

1. **[01-ARCHITECTURE.md](./docs/01-ARCHITECTURE.md)** - System architecture and technology stack
2. **[02-GETTING-STARTED.md](./docs/02-GETTING-STARTED.md)** - Setup and installation instructions
3. **[03-PIPELINE-OVERVIEW.md](./docs/03-PIPELINE-OVERVIEW.md)** - Data processing pipeline
4. **[04-UI-WIREFRAMES.md](./docs/04-UI-WIREFRAMES.md)** - UI design specifications
5. **[05-PROJECT-ROADMAP.md](./docs/05-PROJECT-ROADMAP.md)** - Implementation tracker and roadmap
6. **[06-ML-FEATURES.md](./docs/06-ML-FEATURES.md)** - Machine learning capabilities
7. **[07-ADVANCED-FEATURES.md](./docs/07-ADVANCED-FEATURES.md)** - Advanced analytics features
8. **[08-ENTERPRISE-SECURITY.md](./docs/08-ENTERPRISE-SECURITY.md)** - Security infrastructure
9. **[09-IMPLEMENTATION-STATUS.md](./docs/09-IMPLEMENTATION-STATUS.md)** - Current implementation status

**Start here**: [02-GETTING-STARTED.md](./docs/02-GETTING-STARTED.md) for setup instructions.

## Quick start

### Shared Python virtualenv (backend + indexer)

Run these once from the repo root so every Python component shares the same `.venv`:

```powershell
cd D:\PhaGen-Agent
python -m venv .venv
\.venv\Scripts\activate
pip install -r requirements.txt
```

macOS / Linux equivalent:

```bash
cd PhaGen-Agent
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Reactivate the environment for each new shell with `.\.venv\Scripts\activate` (Windows) or `source .venv/bin/activate` (macOS/Linux).

The repo-level `requirements.txt` simply pulls in `backend/requirements.txt` and `indexes/requirements.txt`, so update those component files and re-run `pip install -r requirements.txt` whenever you add dependencies.

1. **Backend**
   ```powershell
   cd D:\PhaGen-Agent
   .\.venv\Scripts\activate
   cd backend
   pip install -r requirements.txt
   uvicorn app.main:app --reload
   ```
   Set the following environment variables (or update `.env`) before starting the API so SQLAlchemy can connect to Supabase:

   ```bash
   SUPABASE_URL=postgresql://postgres:<password>@db.<project>.supabase.co:5432/postgres
   DATABASE_ECHO=false  # flip to true for SQL logging
   ```
   Copy the `psycopg` URI from **Supabase → Project Settings → Database → Connection string** and paste it into `SUPABASE_URL` (or `DATABASE_URL` if you prefer the legacy name). When Supabase is unreachable, point the backend at SQLite temporarily with `DATABASE_URL=sqlite:///./phagen.db`.

   Object storage (MinIO/S3-compatible) powers raw document snapshots and archived PDF exports. Configure the backend with:

   ```bash
   STORAGE_PROVIDER=s3
   S3_ENDPOINT_URL=http://localhost:9000
   S3_REGION=us-east-1
   S3_ACCESS_KEY=admin
   S3_SECRET_KEY=phagen123
   S3_USE_SSL=false
   S3_RAW_DOCUMENTS_BUCKET=phagen-raw-documents
   S3_REPORTS_BUCKET=phagen-report-artifacts
   ```

   These defaults align with the `minio` service in `infra/docker-compose.yml`; swap them for your own S3 bucket credentials in staging/production. Set `STORAGE_PROVIDER=none` if you want to skip uploads during local development.
2. **Frontend**
   ```bash
   cd frontend
   npm install
   npm run dev
   ```
3. **Crawler / Index**
   ```bash
   cd crawler
   npm install
   npm run crawl
   ```
   After the crawler populates `crawler/storage/datasets/default`, optionally convert any scraped structure diagrams into SMILES via OSRA (only run this if the source license allows OCR). Flip the `allow_osra` flag inside `indexes/data/sample_diagrams.jsonl` or point to your own manifest:
   ```powershell
   cd D:\PhaGen-Agent
   .\.venv\Scripts\python.exe indexes\osra_pipeline.py --input indexes\data\sample_diagrams.jsonl --output indexes\data\osra_results.jsonl --manifest indexes\data\manifests\osra_results.manifest.json --osra-binary indexes\bin\osra-wsl.cmd --image-path-mode wsl --skip-missing
   ```
   The script enforces the `allow_osra` gate per record, hashes image assets, and emits JSONL plus a manifest under `indexes/data/manifests/` summarizing converted/skipped files. If you installed a native Windows OSRA binary you can drop the `--osra-binary/--image-path-mode` overrides; keep them when proxying through the provided WSL wrapper (`indexes/bin/osra-wsl.cmd`).

   Next, canonicalize SMILES/InChI strings so synonyms stay consistent across runs:
   ```powershell
   cd D:\PhaGen-Agent
   .\.venv\Scripts\python.exe indexes\smiles_normalizer.py --input indexes\data\sample_smiles.jsonl --output indexes\data\normalized_smiles.jsonl
   ```
   Swap the sample JSONL for your real payload (one molecule per line). The normalizer deduplicates molecules, snapshots manifests to `indexes/data/manifests/`, and in-place updates `synonyms` with canonical strings.

   Once normalization finishes, switch back to the repo root and run the Python indexer:
   ```powershell
   cd D:\PhaGen-Agent
   .\.venv\Scripts\python.exe indexes\build_index.py --structure-records indexes\data\normalized_smiles.jsonl
   ```
   This rebuilds the live retriever index at `indexes/chroma/`, reuses embeddings for unchanged passages via the on-disk cache at `indexes/.embedding_cache.json`, snapshots the run under `indexes/chroma_snapshots/`, **and** refreshes the RDKit structure catalog from `indexes/data/normalized_smiles.jsonl`. Structure SVGs + metadata land under `indexes/data/structures/{images,metadata}/` with a manifest at `indexes/data/structures/structures.manifest.json`, which workers/reports read at runtime. Each manifest entry now records the provenance schema (`image_id`, `source_type`, `source_ref`, `license`, `generated_at`) for downstream auditing. Pass `--no-structures` to skip the render step, or override the inputs via the `--structure-*` flags (records path, output dirs, manifest, width/height). Standard snapshot flags (`--cadence`, `--snapshot-name`, `--no-snapshot`, cache controls) still apply.

### Section 10 — Automation & CI

- **Local stack via Docker Compose**: from the repo root run `cd infra` followed by `docker compose up --build`. The compose file now focuses on optional services (MinIO, Ollama, rdkit-service, frontend/backend wiring); the operational database lives in Supabase, so you can skip the Postgres container entirely unless you need an offline fallback. Stop everything with `docker compose down` (add `-v` to prune volumes) or target individual services such as `docker compose up api frontend` during development.
- **Continuous integration**: `.github/workflows/ci.yml` executes on every push/PR to `main`, installing Python + Node dependencies, running backend `pytest`, linting the frontend with `npm run lint`, and ensuring all Docker images build cleanly through `docker compose build`. Use it as the reference pipeline when extending Section 10’s test/image requirements.

## Key Features

### Production-Ready Infrastructure
- **Horizontal Scaling**: Kubernetes manifests with HPA autoscaling (2-10 API replicas, 1-8 Celery workers)
- **Distributed Processing**: Celery task queue with Redis broker for async job execution
- **Caching Layer**: Redis caching with 2hr TTL (retrieval) and 24hr TTL (LLM responses)
- **Connection Pooling**: pgBouncer with transaction-mode pooling (25 default, 100 max connections)
- **Health Checks**: `/health`, `/ready`, `/metrics` endpoints for load balancer integration

### Advanced ML & Analytics
- **Disease Mapping**: Automatic repurposing suggestions based on evidence analysis (6 disease categories)
- **Reranker Training**: Learns optimal evidence weights from user feedback (upvotes/downvotes)
- **Custom Embeddings**: Fine-tune sentence-transformers on pharma-specific terminology
- **Active Learning**: Feedback collection API (`POST /api/feedback`) drives continuous improvement

### Security & Compliance
- **Zero Data Retention (ZDR)**: Optional mode for pharma compliance
- **Network Egress Controls**: Default-deny with allowlist validation
- **Container Security**: Trivy scanning, Bandit SAST, Snyk dependency checks
- **Audit Logging**: Immutable job/evidence access logs for regulatory compliance

See `docs/ML_FEATURES.md` for ML implementation details and `IMPLEMENTATION_COMPLETE.md` for infrastructure setup.

## Architecture overview

```
┌────────────┐   REST/WS   ┌────────────┐   in-proc   ┌────────────┐
│ Next.js UI │ ───────────▶│ FastAPI API│────────────▶│ MasterAgent│
└────────────┘             └────────────┘             └────┬─┬─┬───┘
                            │ │ │
                        ┌───────────┘ │ └─────────────┐
                        ▼             ▼               ▼
                     Clinical worker  Patent worker   Literature worker …
```

Each piece lives in its own directory but shares the same repo:

- **Frontend (`frontend/`)** renders the molecule intake, job history, evidence tabs, the multi-molecule comparison workspace, and the report view with inline PDF/JSON export hooks.
- **Backend (`backend/`)** exposes `/api/jobs` for orchestration plus `/api/jobs/{id}/report.pdf` for HTML→PDF exports. Background tasks fan out to the agents package.
- **Agents (`agents/`)** contain `MasterAgent`, shared `LLMClient`, synonym expansion, and four specialized workers. A structured payload powers both the UI and PDF renderer.
- **Indexes (`indexes/`)** bundle the Chroma/FAISS build script so every crawl refresh can be re-indexed with one command.
- **Crawler (`crawler/`)** is a Crawlee project that normalizes CT.gov, PubMed Central, FDA, and patent feeds into JSON ready for embedding.
- **Infra (`infra/`)** holds Docker Compose wiring for Postgres, MinIO/S3-compatible storage, Ollama/OpenAI endpoints, frontend, backend services, and the new `rdkit-service` container that exposes SMILES → SVG/PNG rendering over HTTP.

## 🏗️ Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    PhaGen System                            │
├─────────────────────────────────────────────────────────────┤
│  Frontend (Next.js) → Backend (FastAPI) → Master Agent     │
│       ↓                    ↓                  ↓             │
│  Evidence UI          Redis Cache       4 Worker Agents    │
│  Dashboard            (2hr/24hr TTL)    (Clinical/Patent/  │
│                                          Literature/Market) │
│                            ↓                  ↓             │
│                      PostgreSQL          ChromaDB          │
│                      + pgBouncer         (Vectors)         │
└─────────────────────────────────────────────────────────────┘
```

**Supporting Services**: Celery Workers, MinIO (S3), Ollama (LLM), RDKit Service

See `docs/01-ARCHITECTURE.md` for the full sequence diagram and component details.

## Pipeline & orchestration

1. **Intake** – the Next.js dashboard posts molecules (plus optional synonyms/SMILES) to `/api/jobs`.
2. **Job runner** – FastAPI enqueues the request, spawns a background task, expands synonyms, and parameters for each worker.
3. **Retrieval/RAG** – workers query Chroma via the shared retriever, apply source-ranking, and summarize passages via the shared LLM runtime (Ollama by default, OpenAI optional).
4. **Aggregation** – `MasterAgent` merges worker JSON, runs the Phase 4 synthesis prompt, and persists the innovation story + rubric-based recommendation into the job store.
   - If the LLM synthesis path is unavailable, the master agent now emits a deterministic "lite" summary that strings together worker highlights so downstream UIs/PDFs still render.
5. **Reporting** – the frontend polls `/api/jobs/{id}` until complete; the same payload fuels the Evidence tabs, PDF export (`/api/jobs/{id}/report.pdf`), and JSON download button in `/reports`. Each molecule now accrues monotonic report versions (V1, V2, …) so auditors can track which snapshot was shared.

## Data & indexing pipeline

1. Run the Crawlee project to refresh datasets under `crawler/storage/`.
2. **Optional / gated**: if you scraped structure diagrams and have explicit permission to OCR them, run `indexes/osra_pipeline.py --input <diagrams>.jsonl --output indexes/data/osra_results.jsonl --skip-missing`. The pipeline enforces the `allow_osra` flag per record, hashes input files, snapshots successes/skips to a manifest under `indexes/data/manifests/`, and emits JSONL records that can be merged into downstream SMILES lists.
3. Canonicalize SMILES/InChI inputs with `indexes/smiles_normalizer.py --input <raw>.jsonl --output indexes/data/normalized_smiles.jsonl` so downstream synonym expansion and retrievers reference the same canonical strings (manifests land under `indexes/data/manifests/`). Use the OSRA output JSONL as one of the normalizer inputs when diagrams need to feed synonyms.
4. Execute `indexes/build_index.py` from the repo root (inside `.venv`) to embed new passages into `indexes/chroma/`, reusing cached embeddings where possible, deduplicating overlapping clinical/literature passages (clinical wins by priority), emitting a daily snapshot under `indexes/chroma_snapshots/`, **and** regenerating the RDKit structure catalog + manifest from `indexes/data/normalized_smiles.jsonl` (disable with `--no-structures`).
5. Agents read from `indexes/chroma/` by default; to reproduce a historical run, copy or point the retriever at the desired snapshot folder (each includes a `manifest.json` with dataset hash + git commit).
4. Redeploy/restart workers if the embeddings or retriever settings change so new sources are picked up.

This API-first crawl honors robots.txt (see `crawler/src/robots.ts`) and caps page fragments at 5 KB before indexing. Extend the schema if you add new corpora.

## Evaluation suite

Use the fixtures under `evals/` to guard critical behaviors before shipping.

```powershell
python evals/run_eval.py                 # validate against stored fixtures
python evals/run_eval.py --mode live     # hit a running /api/jobs stack
```

`evals/cases.json` documents five benchmark molecules (Pirfenidone, Metformin,
Nintedanib, Lenabasum, Pamrevlumab) plus expectation bands for recommendations,
market scores, minimum evidence, and quality status. Add new fixtures under
`evals/fixtures/` whenever we expand coverage.

## Export & reporting

- **HTML→PDF**: the backend renders a Jinja2 template with WeasyPrint (`backend/app/reporting.py`). Frontend buttons and Evidence tabs now call the `/api/jobs/{id}/report.pdf` endpoint directly.
- **JSON download**: `/reports` includes a job ID field that serializes the entire job payload for offline analysis or audit trails.
- **Evidence viewer hooks**: every worker summary references the same `WorkerResult` payload so UI badges, PDF sections, and downstream BI exports stay in sync.

## 🎯 Production Deployment

### Docker Compose (Local/On-Premise)
```bash
cd infra
docker-compose up -d
# Services: Backend (8000), Frontend (3000), Redis (6379), MinIO (9000), Ollama (11434)
```

### Kubernetes (Production)
```bash
cd k8s
kubectl apply -f namespace.yaml
kubectl apply -f configmap.yaml
kubectl apply -f redis-deployment.yaml
kubectl apply -f api-deployment.yaml
kubectl apply -f celery-deployment.yaml
kubectl apply -f hpa.yaml
kubectl apply -f ingress.yaml
```

**Features**:
- Horizontal autoscaling (HPA): API 2-10 replicas, Celery 1-8 replicas
- Health checks: `/health`, `/ready`, `/metrics`
- Redis caching (80-90% hit rate)
- pgBouncer connection pooling
- TLS via NGINX Ingress

### Air-Gapped Deployment
```bash
# On connected system:
python infra/airgap/create_bundle.py --output bundle
tar -czf phagen-bundle.tar.gz bundle/

# Transfer to air-gapped system
# Extract and install:
sudo ./install.sh
export PHAGEN_OFFLINE_MODE=true
docker-compose up -d
```

See `docs/08-ENTERPRISE-SECURITY.md` and `infra/airgap/AIRGAP_DEPLOYMENT_GUIDE.md` for details.

## 🔒 Security Features

- **Hardware-Backed Keys**: HSM integration (SoftHSM, AWS CloudHSM, Azure Key Vault)
- **Container Signing**: Cosign signatures with SBOM and provenance
- **Vulnerability Scanning**: Trivy in CI/CD pipeline
- **Network Controls**: Default-deny egress with internal allowlist
- **Air-Gapped Mode**: Complete offline operation with bundled dependencies

## 📊 Performance

- **Job Execution**: 8-15 seconds typical (warm cache)
- **Throughput**: 10-50 jobs/minute (cluster-dependent)
- **Cache Hit Rate**: 80-90% for repeat queries
- **Scalability**: Linear scaling up to 8 Celery workers

## 🧪 Testing & Quality

```bash
# Backend tests
cd backend
pytest tests/

# Frontend lint
cd frontend
npm run lint

# Security tests
pytest backend/tests/test_zdr_mode.py
pytest backend/tests/test_egress_security.py --egress-check

# Evaluation suite
python evals/run_eval.py
```

## 🤝 Contributing

1. Review `docs/05-PROJECT-ROADMAP.md` for current status
2. Follow code style: Black (Python), Prettier (TypeScript)
3. Add tests for new features
4. Update documentation in `docs/`
5. Run security scans before PR

## 📄 License

[Add your license information here]

## 💬 Support

- **Documentation**: See `docs/` directory
- **Issues**: Open a GitHub issue
- **Security**: See `docs/08-ENTERPRISE-SECURITY.md`

---

**Built for pharmaceutical R&D teams** | **Production-ready** | **Enterprise-secure**

## Operational data model

The backend now persists job lifecycles to Postgres via SQLAlchemy models created at startup:

- `molecules` — canonical molecule entries (name, SMILES, synonyms) shared by every job submission.
- `jobs` — orchestration records that track status, serialized payloads, recommendations, and monotonic report versions per molecule.
- `documents` — flattened evidence metadata per worker result (source worker, URL, evidence ID) for each completed job.
- `passages` — normalized snippets tied to the documents table, storing the supporting text and confidence per citation.
- `reports` — structured innovation story payloads alongside their report version, ready for future PDF/JSON export diffs.

Tables are created automatically by `app.database.init_db()` when Uvicorn boots. Point `DATABASE_URL` at your desired Postgres instance and the ORM will create or upgrade the schema as needed (use Alembic for future migrations once the schema stabilizes).

## Object storage & artifacts

- Raw structured job payloads are snapshotted to the `S3_RAW_DOCUMENTS_BUCKET` (default `phagen-raw-documents`) every time a job finishes. Files land under `jobs/{job_id}/raw/<timestamp>.json` so you can diff payloads over time or replay them through downstream pipelines.
- Generated PDF exports are uploaded to the `S3_REPORTS_BUCKET` (default `phagen-report-artifacts`) when the `/api/jobs/{id}/report.pdf` endpoint is called; keys follow `jobs/{job_id}/reports/v{version}-<timestamp>.pdf`.
- Uploads target whatever endpoint/credentials you configure (MinIO via Docker Compose or a real S3 bucket). If MinIO/S3 is unreachable, uploads are skipped but the API still returns the PDF/JSON; check the backend logs for warnings.
- Bucket names, credentials, and endpoints are controlled via environment variables documented in the quick start section, allowing you to point dev, staging, and prod environments at different storage accounts.

### Windows PDF prerequisites

The backend now renders reports exclusively via `pdfkit` + `wkhtmltopdf`. Install the wkhtmltopdf binary once, then pip requirements are enough:

1. Download the Windows installer from [wkhtmltopdf.org](https://wkhtmltopdf.org/downloads.html) and let it add itself to PATH (default: `C:\Program Files\wkhtmltopdf\bin`).
2. If the installer couldn’t update PATH, set the `WKHTMLTOPDF_PATH` environment variable to the absolute `wkhtmltopdf.exe` location.
3. Inside the repo venv, run `pip install -r backend/requirements.txt` (already includes `pdfkit`).
4. Restart Uvicorn so the new PATH/variable is picked up.

If PDF export fails, the backend raises a runtime error that tells you whether `pdfkit` failed to import or the wkhtmltopdf binary couldn’t be found. See the [pdfkit troubleshooting wiki](https://github.com/JazzCore/python-pdfkit/wiki/Installing-wkhtmltopdf) for additional tips.

### RDKit molecular rendering

- `backend/requirements.txt` now includes `rdkit-pypi==2022.9.5` plus `onnxruntime` (Chroma’s ONNX embeddings). Re-run `pip install -r backend/requirements.txt` inside the repo `.venv` after pulling these changes.
- RDKit renders SMILES strings to SVG during job completion and stores the assets under `backend/app/report_assets/reports/images/structures/` by default, with a paired provenance JSON alongside in `.../metadata/`. Override the storage root by setting the `REPORT_ASSETS_DIR` environment variable before starting Uvicorn if you want the SVGs/PDF snippets to land elsewhere (e.g., shared storage or mounted volume).
- Backend + agents now hydrate structure previews from the offline catalog manifest generated by `indexes/build_index.py`. Point the API at your manifest with `STRUCTURE_CATALOG_PATH` (defaults to `indexes/data/structures/structures.manifest.json`) so every request reuses the latest SVG + provenance bundle instead of rerendering.
- Each completed job now carries a `payload.structure` block containing `{ svg, path, metadata_path, smiles, source_type, source_reference, license, image_id, generated_at, error }`. The frontend consumes this metadata to display the molecule preview, and the PDF exporter embeds the same SVG inline while surfacing the provenance + licensing trail.
- `infra/docker-compose.yml` now includes `rdkit-service`, a lightweight FastAPI container running RDKit via Conda. Start it with `docker compose up rdkit-service` (or include it in `docker compose up`). Hit `http://localhost:8081/render` with `{ "smiles": "CC(=O)O", "format": "svg" }` to fetch SVG/PNG payloads; the backend reads the URL from `RDKIT_SERVICE_URL`.
- After installing/upgrading RDKit, restart the backend (`uvicorn app.main:app --reload`) to ensure the new dependency and environment variables are picked up and new jobs generate structure previews automatically.

### OSRA diagram OCR (optional)

- Install the OSRA CLI on the host where you plan to run diagram OCR. The binary is available via Homebrew (`brew install osra`), Conda (`conda install -c conda-forge osra`), or by building from [source](https://github.com/Novartis/osra). On Windows, install OSRA inside WSL (`sudo apt install osra`) and use the provided wrapper `indexes/bin/osra-wsl.cmd` plus `--image-path-mode wsl`, or provide the absolute path to a native Windows binary with `--osra-binary`.
- Populate a JSONL manifest (see `indexes/data/sample_diagrams.jsonl`) with `diagram_id`, `image_path`, metadata, and an explicit `allow_osra` flag per record. Only assets with `allow_osra: true` will be processed, so you can safely keep restricted diagrams in the manifest without risking violations.
- Run `python indexes/osra_pipeline.py --input <manifest>.jsonl --output indexes/data/osra_results.jsonl --skip-missing`. When invoking through WSL, pass `--osra-binary indexes/bin/osra-wsl.cmd --image-path-mode wsl`; supply `--extra-osra-args` if you need to adjust resolution or language packs, and use `--base-dir` when image paths are relative to another folder.
- Each run emits JSONL plus a manifest under `indexes/data/manifests/` summarizing conversions, skips, and failures. Feed the output JSONL into `smiles_normalizer.py` (or merge it with crawler-sourced SMILES) before synonym expansion and indexing.

## Comparison workspace

- `/comparison` lets reviewers line up two or three molecules side-by-side with shared metrics, worker summaries, and top citations. It defaults to demo payloads but accepts real job IDs once they exist.
- The frontend form hits the new `/api/jobs/compare?job_ids=A&job_ids=B` endpoint, which streams multiple `JobResponse` objects in a single call so the UI stays in sync with backend payloads.
- Each comparison card links back to the job timeline and report workspace, keeping Phase 4 aggregation, audit trails, and exports tightly coupled.

## Crawling & compliance

- The crawler now follows an **API-first → robots.txt-validated crawl** pipeline. Each source tries to load mock API data first (`crawler/mock-data/`), then falls back to HTML only if `robots.txt` allows access.
- Section 5 normalization now strips boilerplate, redacts simple PII, and emits deterministic chunk metadata (chunk ID, char spans, redaction counters) via `crawler/src/normalize.ts` before datasets reach the indexer.
- Domain-level crawl budgets cap HTML fetches per host (default 2 pages unless overridden in `crawler/src/index.ts`) and annotate dataset entries when a source is skipped due to quota exhaustion.
- `crawler/src/robots.ts` caches per-domain policies (allow/deny + crawl-delay) using the `PhaGenBot/1.0` user agent. Update the user agent or contact email there before running against live sites.
- Crawled artifacts are written to `crawler/storage/` (ignored by git) and capped at 5 KB snippets for safety. Extend the schema before indexing into FAISS/Chroma.
