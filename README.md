# PhaGen Agentic

PhaGen Agentic is a hackathon-ready, agentic multi-worker platform that compresses molecule repurposing research from months to minutes. The repo is organized as a monorepo so each domain (frontend, backend, agents, crawler, infrastructure, docs) shares one source of truth.

## Current status

| Area | State |
| --- | --- |
| Frontend | Basic Next.js scaffold with molecule input and job timeline placeholder |
| Backend | FastAPI service with job orchestration endpoints backed by Postgres persistence |
| Agents | Master agent plus four worker stubs wired to deterministic mock data and RAG hooks |
| Crawler | Crawlee TypeScript project with seed list plus Section 5 normalization (boilerplate stripping, PII redaction, chunk metadata) |
| Infra | Docker Compose file for local stack (Postgres, MinIO, Ollama, rdkit-service, services) |

The MVP mirrors the implementation plan in `ignore.md`. Build phases are mapped to repo folders so each team can work in parallel.

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
   Set the following environment variables (or update `.env`) before starting the API so SQLAlchemy can connect to Postgres:

   ```bash
   DATABASE_URL=postgresql+psycopg://phagen:phagen@localhost:5432/phagen
   DATABASE_ECHO=false  # flip to true for SQL logging
   ```
   The default URL matches the `postgres` service in `infra/docker-compose.yml`. Update the value if you run Postgres elsewhere (e.g., managed cloud instance) or temporarily point it at `sqlite:///./phagen.db` when developing without a database container.

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
   After the crawler populates `crawler/storage/datasets/default`, switch back to the repo root and run the Python indexer:
   ```powershell
   cd D:\PhaGen-Agent
   .\.venv\Scripts\python.exe indexes\build_index.py
   ```
   This rebuilds the live retriever index at `indexes/chroma/`, reuses embeddings for unchanged passages via the on-disk cache at `indexes/.embedding_cache.json`, and snapshots the run under `indexes/chroma_snapshots/` using a daily timestamp. Use `--cadence monthly` or `--snapshot-name my-run` for custom folders, `--no-snapshot` to skip copies, and `--no-cache`/`--cache-path` if you need to bypass or relocate the embedding cache.

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

See `docs/architecture.md` for the full sequence diagram and responsibilities per component.

## Pipeline & orchestration

1. **Intake** – the Next.js dashboard posts molecules (plus optional synonyms/SMILES) to `/api/jobs`.
2. **Job runner** – FastAPI enqueues the request, spawns a background task, expands synonyms, and parameters for each worker.
3. **Retrieval/RAG** – workers query Chroma via the shared retriever, apply source-ranking, and summarize passages via the shared LLM runtime (Ollama by default, OpenAI optional).
4. **Aggregation** – `MasterAgent` merges worker JSON, runs the Phase 4 synthesis prompt, and persists the innovation story + rubric-based recommendation into the job store.
   - If the LLM synthesis path is unavailable, the master agent now emits a deterministic "lite" summary that strings together worker highlights so downstream UIs/PDFs still render.
5. **Reporting** – the frontend polls `/api/jobs/{id}` until complete; the same payload fuels the Evidence tabs, PDF export (`/api/jobs/{id}/report.pdf`), and JSON download button in `/reports`. Each molecule now accrues monotonic report versions (V1, V2, …) so auditors can track which snapshot was shared.

## Data & indexing pipeline

1. Run the Crawlee project to refresh datasets under `crawler/storage/`.
2. Execute `indexes/build_index.py` from the repo root (inside `.venv`) to embed new passages into `indexes/chroma/`, reusing cached embeddings where possible, deduplicating overlapping clinical/literature passages (clinical wins by priority), and emitting a daily snapshot under `indexes/chroma_snapshots/`.
3. Agents read from `indexes/chroma/` by default; to reproduce a historical run, copy or point the retriever at the desired snapshot folder (each includes a `manifest.json` with dataset hash + git commit).
4. Redeploy/restart workers if the embeddings or retriever settings change so new sources are picked up.

This API-first crawl honors robots.txt (see `crawler/src/robots.ts`) and caps page fragments at 5 KB before indexing. Extend the schema if you add new corpora.

## Export & reporting

- **HTML→PDF**: the backend renders a Jinja2 template with WeasyPrint (`backend/app/reporting.py`). Frontend buttons and Evidence tabs now call the `/api/jobs/{id}/report.pdf` endpoint directly.
- **JSON download**: `/reports` includes a job ID field that serializes the entire job payload for offline analysis or audit trails.
- **Evidence viewer hooks**: every worker summary references the same `WorkerResult` payload so UI badges, PDF sections, and downstream BI exports stay in sync.

## Validation & traceability

- Every evidence snippet now receives a deterministic ID (e.g., `clinical-1`) when the master agent serializes results.
- The innovation story is split into sentence-level claims, each linked to one or more evidence IDs; the payload exposes this under `validation` with pass/fail status plus linked counts.
- `/comparison`, `/results`, and the PDF report highlight the validation summary so reviewers can confirm every claim is grounded before sharing deliverables.

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
- Each completed job now carries a `payload.structure` block containing `{ svg, path, metadata_path, smiles, source_type, source_reference, image_id, generated_at, error }`. The frontend consumes this metadata to display the molecule preview, and the PDF exporter embeds the same SVG inline while surfacing the provenance trail.
- `infra/docker-compose.yml` now includes `rdkit-service`, a lightweight FastAPI container running RDKit via Conda. Start it with `docker compose up rdkit-service` (or include it in `docker compose up`). Hit `http://localhost:8081/render` with `{ "smiles": "CC(=O)O", "format": "svg" }` to fetch SVG/PNG payloads; the backend reads the URL from `RDKIT_SERVICE_URL`.
- After installing/upgrading RDKit, restart the backend (`uvicorn app.main:app --reload`) to ensure the new dependency and environment variables are picked up and new jobs generate structure previews automatically.

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
