# PhaGen Agentic Architecture

```
┌───────────┐     ┌───────────┐     ┌─────────────┐
│ Frontend  │ --> │  Backend  │ --> │ Master Agent│
└───────────┘     └───────────┘     └─────────────┘
       \                               / |  |  \
        \                             /  |  |   \
         v                           v   v  v    v
      RAG Viewer            Clinical  Patent Literature Market Workers
```

1. **Frontend (Next.js)**
   - Molecule input form.
   - Job poller visualizes orchestration timeline.
   - Evidence payload viewer for hackathon demos.

2. **Backend (FastAPI)**
   - `/api/jobs` enqueue endpoint.
   - Background task triggers the master agent and persists results in the in-memory store.
   - `/api/jobs/{id}` returns job status and payload.

3. **Agents package**
   - `MasterAgent` orchestrates four workers.
   - Each worker currently uses a mock `Retriever`; swap with FAISS/Chroma once indexes exist.
   - Outputs structured summaries + metadata + evidence for PDF/report use.

4. **Crawler**
   - Crawlee TypeScript job for harvesting seed datasets.
   - Pushes normalized docs to dataset or disk for indexing.

5. **Infra**
   - Docker Compose (coming soon) wires Postgres, MinIO, Ollama, backend, and frontend.

## Data flow checklist

1. User submits molecule via UI.
2. Backend enqueues job and responds with `job_id`.
3. Background worker runs master agent → hits retriever → uses worker prompts.
4. Aggregated result stored in job store → frontend polls until completion.
5. UI renders innovation story and evidence JSON. PDF worker (Phase 4) will consume same payload.
