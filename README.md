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

1. **Backend**
   ```bash
   cd backend
   uv venv .venv ; .venv\Scripts\activate
   uv pip install -r requirements.txt
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

See `docs/architecture.md` for a diagram plus detailed flow.
