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
   - Section 5 normalization path strips boilerplate, redacts light PII, and chunks sources into metadata-rich passages before indexing.
   - Domain budgets enforce per-host page caps (with logs + dataset annotations) so we stay polite beyond robots.txt.
   - Index builder deduplicates overlapping clinical/literature passages so the retriever surfaces a single canonical chunk per evidence item.

5. **Infra**
   - `infra/docker-compose.yml` wires Postgres, MinIO, Ollama, backend, frontend, and the RDKit renderer onto a single bridge network so the full stack spins up with one `docker compose up`.
   - `.github/workflows/ci.yml` (Section 10) runs backend pytest, frontend lint, and `docker compose build` to guarantee the same images that ship locally also build in CI.
   - Includes `rdkit-service` (FastAPI + RDKit) for on-demand SMILES rendering consumed by backend workers.

## Data flow checklist

1. User submits molecule via UI.
2. Backend enqueues job and responds with `job_id`.
3. Background worker runs master agent → hits retriever → uses worker prompts.
4. Aggregated result stored in job store → frontend polls until completion.
5. UI renders innovation story and evidence JSON. PDF worker (Phase 4) will consume same payload.
