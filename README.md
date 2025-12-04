# PhaGen Agentic

PhaGen Agentic is a hackathon-ready, agentic multi-worker platform that compresses molecule repurposing research from months to minutes. The repo is organized as a monorepo so each domain (frontend, backend, agents, crawler, infrastructure, docs) shares one source of truth.

## Current status

| Area | State |
| --- | --- |
| Frontend | Basic Next.js scaffold with molecule input and job timeline placeholder |
| Backend | FastAPI service with job orchestration endpoints and mock job store |
| Agents | Master agent plus four worker stubs wired to deterministic mock data and RAG hooks |
| Crawler | Crawlee TypeScript project with seed list + document normalization utilities |
| Infra | Docker Compose file for local stack (Postgres, MinIO, Ollama, services) |

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
   This writes embeddings to `indexes/chroma/` (git-ignored). Re-run the script any time new crawl output lands so FAISS/Chroma stays current.

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
- **Infra (`infra/`)** holds Docker Compose wiring for Postgres, MinIO/S3-compatible storage, Ollama/OpenAI endpoints, frontend, and backend services.

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
2. Execute `indexes/build_index.py` from the repo root (inside `.venv`) to embed new passages into `indexes/chroma/`.
3. Point the agents retriever at the fresh Chroma snapshot (default path already aligns with `indexes/chroma/`).
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

### Windows PDF prerequisites

The backend now renders reports exclusively via `pdfkit` + `wkhtmltopdf`. Install the wkhtmltopdf binary once, then pip requirements are enough:

1. Download the Windows installer from [wkhtmltopdf.org](https://wkhtmltopdf.org/downloads.html) and let it add itself to PATH (default: `C:\Program Files\wkhtmltopdf\bin`).
2. If the installer couldn’t update PATH, set the `WKHTMLTOPDF_PATH` environment variable to the absolute `wkhtmltopdf.exe` location.
3. Inside the repo venv, run `pip install -r backend/requirements.txt` (already includes `pdfkit`).
4. Restart Uvicorn so the new PATH/variable is picked up.

If PDF export fails, the backend raises a runtime error that tells you whether `pdfkit` failed to import or the wkhtmltopdf binary couldn’t be found. See the [pdfkit troubleshooting wiki](https://github.com/JazzCore/python-pdfkit/wiki/Installing-wkhtmltopdf) for additional tips.

### RDKit molecular rendering

- `backend/requirements.txt` now includes `rdkit-pypi==2022.9.5` plus `onnxruntime` (Chroma’s ONNX embeddings). Re-run `pip install -r backend/requirements.txt` inside the repo `.venv` after pulling these changes.
- RDKit renders SMILES strings to SVG during job completion and stores the assets under `backend/app/report_assets/structures/` by default. Override the storage root by setting the `REPORT_ASSETS_DIR` environment variable before starting Uvicorn if you want the SVGs/PDF snippets to land elsewhere (e.g., shared storage).
- Each completed job now carries a `payload.structure` block containing `{ svg, path, smiles, error }`. The frontend consumes this metadata to display the molecule preview, and the PDF exporter embeds the same SVG inline.
- After installing/upgrading RDKit, restart the backend (`uvicorn app.main:app --reload`) to ensure the new dependency and environment variables are picked up and new jobs generate structure previews automatically.

## Comparison workspace

- `/comparison` lets reviewers line up two or three molecules side-by-side with shared metrics, worker summaries, and top citations. It defaults to demo payloads but accepts real job IDs once they exist.
- The frontend form hits the new `/api/jobs/compare?job_ids=A&job_ids=B` endpoint, which streams multiple `JobResponse` objects in a single call so the UI stays in sync with backend payloads.
- Each comparison card links back to the job timeline and report workspace, keeping Phase 4 aggregation, audit trails, and exports tightly coupled.

## Crawling & compliance

- The crawler now follows an **API-first → robots.txt-validated crawl** pipeline. Each source tries to load mock API data first (`crawler/mock-data/`), then falls back to HTML only if `robots.txt` allows access.
- `crawler/src/robots.ts` caches per-domain policies (allow/deny + crawl-delay) using the `PhaGenBot/1.0` user agent. Update the user agent or contact email there before running against live sites.
- Crawled artifacts are written to `crawler/storage/` (ignored by git) and capped at 5 KB snippets for safety. Extend the schema before indexing into FAISS/Chroma.
