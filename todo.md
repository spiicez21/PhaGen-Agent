# PhaGen Agentic — Execution Tracker

This to-do list is distilled from the implementation roadmap in `ignore.md`. Tasks we have already finished are captured under **Completed**, while remaining items stay under **Pending** so the team can pick them up quickly.

## ✅ Completed
- **Phase 0 (Prep):** Monorepo scaffold (`frontend/`, `backend/`, `agents/`, `crawler/`, `indexes/`, `infra/`, `docs/`) plus executive-summary README are in place.
- **Phase 1 (Core orchestration & mocks):** Master agent skeleton with four worker stubs, deterministic mock data, and the basic Next.js dashboard + FastAPI backend flow are running end-to-end.
- **Phase 2 (Crawling foundation):** Crawlee project configured with API-first sources, robots-aware fallback, and datasets refreshed via `npm run crawl`.
- **Phase 2 (Indexer):** Python Chroma/FAISS pipeline (`indexes/build_index.py`) now ingests crawler output inside the shared `.venv` and persists embeddings to `indexes/chroma/`.
- **Phase 2 (Retriever):** Agents call a real Chroma-backed retriever with per-worker context budgeting instead of the previous mock snippets.
- **Phase UI (Core screens):** Landing hero, molecule intake, job status timeline, innovation story summary, and multi-panel evidence dashboard now reflect ~70% of the UI blueprint.

## 🔜 Pending
### Phase 3 — Worker intelligence & grounding
- [ ] ClinicalTrials worker: parse trials (phase, status, endpoints, populations) from indexed passages.
- [ ] Patent & Regulatory worker: surface assignees, priority dates, blocking claims, and contraindication notes.
- [ ] Literature worker: summarize mechanism-of-action evidence and cite passages.
- [ ] Market worker: produce market score based on incidence/prevalence proxies plus competitor context.
- [ ] Integrate chosen LLM runtime (Ollama Gemma2:2B or remote) with the worker prompts and retriever output.

### Phase 4 — Aggregation & reporting
- [ ] Master agent final synthesis prompt that merges worker JSON into the "Innovation Story" + recommendation rubric.
- [ ] PDF generation flow (HTML → PDF) with evidence viewer hooks and download button in the frontend.
- [ ] Evidence dashboard polish: tabs per worker, confidence badges, citation links back to sources.

### Phase UI — Experience & admin
- [ ] Add Full Report viewer (PDF/HTML) with download/share/export controls.
- [ ] Build history/saved runs workspace with filters + download actions.
- [ ] Create crawler status dashboard and dataset/index manager for admins.
- [ ] Implement stretch visualizations (knowledge graph, citation trace, molecule comparison) after MVP.

### Data & infra follow-ups
- [ ] Flesh out Postgres schema (`molecules`, `jobs`, `documents`, `passages`, `reports`) and wire persistence into backend.
- [ ] Add MinIO/S3 storage wiring for raw documents and PDF artifacts.
- [ ] Extend crawler normalization (chunking, boilerplate stripping, PII redaction) per Section 5 guidance.

### Quality, guardrails, and ops
- [ ] Implement retrieval precision checks / coverage metrics plus guardrails (evidence thresholds, anomaly detection).
- [ ] Set up Docker Compose & GitHub Actions workflow described in Section 10 (tests + image builds).
- [ ] Prepare the hackathon demo script, slides, and fallback assets outlined in Section 11.

### Stretch / future
- [ ] Track stretch ideas from Section 12 (active learning loop, patent semantic matching, continuous crawler, etc.) once the MVP is stable.
