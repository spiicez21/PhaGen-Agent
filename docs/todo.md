# PhaGen Agentic — Execution Tracker

This to-do list is distilled from the implementation roadmap in `ignore.md`. Tasks we have already finished are captured under **Completed**, while remaining items stay under **Pending** so the team can pick them up quickly.

## ✅ Completed
### Phase 0 — Prep
- Monorepo scaffold (`frontend/`, `backend/`, `agents/`, `crawler/`, `indexes/`, `infra/`, `docs/`) plus executive-summary README are in place.

### Phase 1 — Core orchestration & mocks
- Master agent skeleton with four worker stubs, deterministic mock data, and the basic Next.js dashboard + FastAPI backend flow are running end-to-end.

### Phase 2 — Data plumbing
- Crawlee project configured with API-first sources, robots-aware fallback, and datasets refreshed via `npm run crawl`.
- Python Chroma/FAISS pipeline (`indexes/build_index.py`) now ingests crawler output inside the shared `.venv` and persists embeddings to `indexes/chroma/`.
- Agents call a real Chroma-backed retriever with per-worker context budgeting instead of the previous mock snippets.

### Phase 3 — Worker intelligence
- ClinicalTrials agent now parses registry feeds into structured trials with phase, status, endpoints, and populations surfaced in metadata.
- Patent agent now surfaces assignees, priority dates, blocking claim summaries, and FDA contraindication notes directly from retrieved passages.
- Literature worker now summarizes mechanism-of-action evidence, cites DOI-backed passages, and exposes structured metadata for the dashboard.

### Phase 4 — Aggregation & reporting
- Evidence tabs ship with confidence badges, metadata grids, and citation links for each worker panel so reviewers can drill into sources.

### Phase UI — Experience & admin
- Landing hero, molecule intake, job status timeline, innovation story summary, and multi-panel evidence dashboard now reflect ~70% of the UI blueprint.
- Full report workspace renders PDF/JSON export controls plus sectionized content, matching the deliverable spec.
- Saved-run workspace lists molecules with status/recommendation, quick open, and PDF download actions.
- Crawler status console and dataset/index manager expose queue metrics, robots summaries, rebuild, and purge controls.

## 🔜 Pending
### Phase 3 — Worker intelligence & grounding
- [ ] Market worker: produce market score based on incidence/prevalence proxies plus competitor context.
- [ ] Integrate chosen LLM runtime (Ollama Gemma2:2B or remote) with the worker prompts and retriever output.
- [ ] Synonym expansion worker: convert SMILES → canonical name → synonyms/aliases for better retrieval.
- [ ] Source-ranking scorecard: prioritize clinical > regulatory > literature > patent passages when grounding workers.
- [ ] Worker timeouts & retries: add per-worker timeout window plus retry budget for flaky sources.
- [ ] Confidence calibration: normalize each worker's 0–1 confidence into Low/Med/High bands before surfacing.

### Phase 4 — Aggregation & reporting
- [ ] Master agent final synthesis prompt that merges worker JSON into the "Innovation Story" + recommendation rubric.
- [ ] PDF generation flow (HTML → PDF) with evidence viewer hooks and download button in the frontend.
- [ ] Multi-molecule comparison mode for side-by-side evidence (stretch demo goal).
- [ ] Validation pass linking every claim in the story to evidence IDs.
- [ ] Fallback summarizer that emits a lite summary if the LLM path fails.
- [ ] Report versioning so each molecule keeps V1/V2 snapshots for auditability.

### Phase UI — Experience & admin
- [ ] Implement stretch visualizations (knowledge graph, citation trace, molecule comparison) after MVP.
- [ ] Pagination + filtering for saved runs (status, date range, recommendation).
- [ ] Standardize confidence color palette across all evidence components.
- [ ] Worker-specific error states (timeout, missing data, robots block) in the UI.
- [ ] Visual diff viewer to compare report versions (V1 vs V2 changes).

### Data & infra follow-ups
- [ ] Flesh out Postgres schema (`molecules`, `jobs`, `documents`, `passages`, `reports`) and wire persistence into backend.
- [ ] Add MinIO/S3 storage wiring for raw documents and PDF artifacts.
- [ ] Extend crawler normalization (chunking, boilerplate stripping, PII redaction) per Section 5 guidance.
- [ ] Index snapshotting cadence (daily/monthly builds) for reproducibility.
- [ ] Embedding cache to avoid re-embedding unchanged passages.
- [ ] Domain-level crawl budgets to limit pages per source.
- [ ] Evidence deduplication across clinical/literature corpora.

### Quality, guardrails, and ops
- [ ] Implement retrieval precision checks / coverage metrics plus guardrails (evidence thresholds, anomaly detection).
- [ ] Set up Docker Compose & GitHub Actions workflow described in Section 10 (tests + image builds).
- [ ] Prepare the hackathon demo script, slides, and fallback assets outlined in Section 11.
- [ ] Model hallucination detector or citation-gap checker for final stories.
- [ ] API call budget monitor for NCBI/OpenFDA rate limits.
- [ ] Standardize LLM temperature per worker for consistent tone.
- [ ] Build an evals suite (5–10 sample molecules with expected outputs).
- [ ] SLA-style logging: latency, token usage, and failure counts per worker.

### Stretch / future
- [ ] Track stretch ideas from Section 12 (active learning loop, patent semantic matching, continuous crawler, etc.) once the MVP is stable.
- [ ] Horizontal autoscaling for agent workers.
- [ ] External `/api/analyze-molecule` endpoint for partners.
- [ ] Collaboration surface (annotations, shared notes, mentions).
- [ ] Molecule–disease mapping model for proactive repurposing suggestions.
