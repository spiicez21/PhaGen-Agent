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

See `docs/architecture.md` for a diagram plus detailed flow.

## Crawling & compliance

- The crawler now follows an **API-first → robots.txt-validated crawl** pipeline. Each source tries to load mock API data first (`crawler/mock-data/`), then falls back to HTML only if `robots.txt` allows access.
- `crawler/src/robots.ts` caches per-domain policies (allow/deny + crawl-delay) using the `PhaGenBot/1.0` user agent. Update the user agent or contact email there before running against live sites.
- Crawled artifacts are written to `crawler/storage/` (ignored by git) and capped at 5 KB snippets for safety. Extend the schema before indexing into FAISS/Chroma.
