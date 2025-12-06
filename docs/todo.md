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
- Market worker now synthesizes TAM/incidence/competition context into a scored summary with structured metadata for the dashboard.
- Workers now call the shared LLM runtime (Ollama Gemma2:2B by default) to synthesize summaries directly from retriever output with graceful fallbacks.
- Synonym expansion worker now normalizes SMILES/canonical names into alias lists that feed every retriever query.
- Source-ranking scorecard now prioritizes clinical > regulatory > literature > patent passages before grounding summaries.
- Worker execution now enforces per-worker timeouts plus a small retry budget for flaky sources.
- Confidence calibration now maps worker confidence into Low/Med/High bands before surfacing results.

### Phase 4 — Aggregation & reporting
- Evidence tabs ship with confidence badges, metadata grids, and citation links for each worker panel so reviewers can drill into sources.
- Master agent now uses an LLM prompt that merges worker JSON into the Innovation Story and rubric-based recommendation.
- HTML-to-PDF export path now streams backend-generated reports with frontend download hooks.
- Multi-molecule comparison workspace now renders side-by-side evidence plus a `/api/jobs/compare` endpoint for shared payloads.
- Validation pass now assigns deterministic evidence IDs and links every innovation story claim to its supporting citations, surfacing the status in UI + PDF exports.
- Report versioning now tracks V1/V2 snapshots per molecule for downstream audits.
- Lite fallback summarizer now emits a multi-line innovation story when the LLM synthesis path is unavailable.
- RDKit rendering now materializes SMILES into SVG assets saved under `backend/app/report_assets/reports/images/structures/` with paired provenance JSON in `.../metadata/`, exposes `payload.structure` metadata (svg/path/metadata_path/smiles/source fields/image_id), and embeds the preview in both the frontend dashboard and PDF export.
- RDKit canonical rendering now has pytest-based visual regression tests with baseline SVGs to guard against future RDKit changes.
- PubChem REST fallback now backstops RDKit failures by pulling SVGs via PUG APIs, persisting them under the same reports/images schema with provenance metadata.

### Phase UI — Experience & admin
- Landing hero, molecule intake, job status timeline, innovation story summary, and multi-panel evidence dashboard now reflect ~70% of the UI blueprint.
- Full report workspace renders PDF/JSON export controls plus sectionized content, matching the deliverable spec.
- Saved-run workspace lists molecules with status/recommendation, quick open, and PDF download actions.
- Crawler status console and dataset/index manager expose queue metrics, robots summaries, rebuild, and purge controls.

### Data & infra follow-ups
- Postgres schema now covers `molecules`, `jobs`, `documents`, `passages`, and `reports`, and the FastAPI backend persists jobs plus evidence into the database via SQLAlchemy models.
- MinIO/S3 wiring now stores raw job payloads and generated PDF artifacts via the new object storage client, with configurable buckets and docker-compose defaults.
- Crawlee normalization pipeline now strips boilerplate, redacts PII, and chunks sources into metadata-rich passages per Section 5 guidance, feeding indexes with deterministic chunk IDs + redaction stats.
- Index snapshotter now runs via `indexes/build_index.py` (daily by default, monthly optional) and writes timestamped copies with manifests (dataset hash, git commit, record counts) under `indexes/chroma_snapshots/` while keeping the live retriever under `indexes/chroma/`.
- Embedding cache now stores MiniLM vectors keyed by doc hash in `indexes/.embedding_cache.json`, letting the index builder reuse embeddings for unchanged passages and only encode the delta.
- Domain-level crawl budgets now cap HTML fetches per host (default 2 unless overridden), log usage, and annotate dataset entries when a source is skipped for exceeding its quota.
- Evidence deduplication now drops lower-priority duplicates between clinical and literature corpora (clinical wins), annotating kept/dropped entries via `metadata.dedup` and surfacing stats inside snapshot manifests.
- RDKit rendering now has a dedicated `rdkit-service` container (Conda + RDKit + Pillow + FastAPI) wired into `docker-compose.yml` for on-demand SMILES → SVG/PNG rendering via `POST /render`.
- SMILES normalizer job (`python indexes/smiles_normalizer.py`) canonicalizes SMILES/InChI inputs via RDKit, dedupes aliases, and emits JSONL + manifest files before indexing and synonym expansion.
- Optional OSRA pipeline (`python indexes/osra_pipeline.py`) walks approved structure diagrams, enforces `allow_osra` gates, hashes image artifacts, converts diagrams to SMILES via the OSRA CLI, and emits JSONL + manifest outputs for downstream normalization.
- RDKit structure rendering now runs automatically inside `indexes/build_index.py`, consuming `indexes/data/normalized_smiles.jsonl`, emitting SVG/metadata assets + manifests under `indexes/data/structures/`, and wiring summary stats into snapshot manifests so workers/reports can reuse the catalog without re-rendering.
- Structure manifest entries now carry full provenance metadata (`image_id`, `source_type`, `source_ref`, `license`, `generated_at`) so downstream workers/reports can trace each asset.
- Docker Compose deployment templates now support local/on-prem hosting with `infra/docker-compose.yml` orchestrating MinIO, Ollama, RDKit service, Redis cache, pgBouncer connection pooler, and application services without external cloud dependencies.
- Private model hosting with Ollama now runs locally with configurable LLM models (smollm:360m, gemma2:2b, etc.) via `LLM_MODEL` and `LLM_BASE_URL` environment variables, with no external API calls required.
- Kubernetes deployment manifests now provide production-ready orchestration with `k8s/` directory containing Deployments, Services, ConfigMaps, Secrets, HPA autoscaling (API 2-10 replicas, Celery 1-8 replicas), NGINX Ingress with TLS, and health check endpoints (/health, /ready, /metrics).
- Redis caching layer now accelerates retrieval and LLM calls with `backend/app/cache.py` providing decorators for automatic caching (retrieval: 2hr TTL, LLM: 24hr TTL), graceful degradation when Redis unavailable, and pattern-based cache invalidation.
- Celery distributed task queue now enables horizontal scaling via `backend/app/celery_tasks.py` with Redis broker/backend, async job execution across multiple workers, periodic cleanup tasks via Celery Beat, and job API supports `?use_celery=true` flag for distributed execution.
- pgBouncer connection pooling now optimizes database connections via `infra/pgbouncer/` config with transaction-mode pooling (25 default pool, 100 max connections), reducing database load for high-concurrency scenarios.
- Evidence feedback API now collects user ratings via `POST /api/feedback` for active learning, storing upvote/downvote/flag feedback in `evidence_feedback` table with aggregate statistics endpoint (`GET /api/feedback/stats`) for reranker training.

### Quality, guardrails, and ops
- Retrieval precision coverage is now measured per worker (queries attempted, passages gathered, evidence counts, precision proxy), and guardrails raise alerts for low evidence, weak coverage, or anomalous retrieval gaps. Alerts and metrics ship in each job payload under the `quality` block so the UI/API can flag risky runs.
- Docker Compose stack + GitHub Actions workflow (Section 10) now mirror the ops plan: `infra/docker-compose.yml` boots the optional services (MinIO/Ollama/RDKit/backend/frontend) while Supabase hosts Postgres, and `.github/workflows/ci.yml` runs backend pytest, frontend lint, and `docker compose build` so image builds and tests gate every PR.
- Citation-gap detector now scans every innovation-story claim against its cited evidence, tracking lexical overlap and numeric support so hallucination-prone sentences trigger explicit alerts in the quality payload + UI.
- API call budget monitor now tracks NCBI/OpenFDA usage per minute/day, escalates status into the job quality payload, and exposes a UI panel warning reviewers before rate-limit breaches.
- LLM temperature map is now standardized per worker + master synthesis step, with env overrides (`LLM_TEMP_DEFAULT`, `LLM_TEMP_<WORKER>`, `LLM_TEMP_MASTER`) to keep tone consistent across runs.
- Eval suite now lives under `evals/` with five benchmark molecules, JSON fixtures, and `run_eval.py` so CI/local devs can verify recommendations, guardrails, and evidence coverage before merging.
- SLA-style logging now tracks latency (retrieval/LLM breakdown), token usage, and failure counts per worker via `WorkerSLAMetrics` class, with metrics exposed in job payloads under `result.metrics["sla"]` for observability dashboards.
- Security test suite now validates ZDR mode (`test_zdr_mode.py` - temp file purging, S3 write prevention), network egress controls (`test_egress_security.py` - external host blocking, internal allowlist), and container vulnerabilities (GitHub Actions workflow with Trivy scanning on rdkit-service and backend images, static analysis with Bandit, dependency scanning with Snyk/Safety).
- Network egress controls now enforce default-deny for external connections with allowlist for internal services (localhost, Docker network), validated via pytest with `--egress-check` flag that fails CI if external access succeeds.

## 🔜 Pending
### Phase UI — Experience & admin
- [ ] Implement stretch visualizations (knowledge graph, citation trace, molecule comparison) after MVP.
- [ ] Pagination + filtering for saved runs (status, date range, recommendation).
- [ ] Standardize confidence color palette across all evidence components.
- [ ] Worker-specific error states (timeout, missing data, robots block) in the UI.
- [ ] Visual diff viewer to compare report versions (V1 vs V2 changes).
- [ ] Molecule structure viewer component (SVG viewer + zoom + download) powered by RDKit-generated SVGs.
- [ ] Add "Request structure" action in molecule intake to create SMILES → RDKit render.

### Quality, guardrails, and ops
- [ ] Prepare the hackathon demo script, slides, and fallback assets outlined in Section 11.

#### ➕ Added — Security & Compliance (Core)
- [ ] Zero Data Retention (ZDR) mode flag that disables persistence of input documents/user uploads (clear memory + no S3 writes).
- [ ] VPC / On-Prem / Self-Host deployment option (Terraform + Docker Compose templates).
- [ ] Tenant isolation & RBAC with tenant-scoped DB schema/row-level permissions plus UI/API roles.
- [ ] End-to-end encryption: TLS 1.3 in transit, AES-256 at rest (DB/S3) with KMS-managed keys.
- [ ] Centralized secrets management (Vault or cloud KMS) with rotating keys.
- [ ] Immutable audit logs for worker runs, data access, and report generation (timestamps, user, job_id) with retention policy.
- [ ] PII redaction & DLP to scrub sensitive fields before storage and block exfiltration of PHI/email identifiers.
- [x] Network egress controls enforcing default-deny for customer-hosted agents (see Stretch section).
- [x] Integrated SAST/DAST & dependency scanning (Snyk/Dependabot) in CI (Trivy, Bandit, Safety in security-tests.yml).
- [ ] Pen testing & vulnerability management cadence (quarterly pen tests, critical CVE patch SLA).
- [ ] Incident response + breach communication plan with notify list and SLA.
- [ ] Compliance artifacts: SOC2 controls, ISO27001 checklist, HIPAA/GxP mapping for customer audits.
- [ ] Data residency options plus DPA/BAA templates for enterprise deals.
- [ ] Monitoring & SIEM integration forwarding audit logs to Splunk/Datadog/Elastic.
- [ ] Periodic access reviews & enforced MFA for admin roles.
- [ ] Privacy & legal review pack (privacy FAQ, data use statement, legal review dossier).

#### ➕ Added — Security Tests & Validation
- [x] Automated test verifying ZDR mode purges temp files and prevents S3 writes.
- [x] CI job running static analysis plus container vulnerability scan (Trivy) on `rdkit-service` and other images.
- [x] E2E security smoke test that attempts egress from agents; CI fails if egress succeeds in locked deployments.

### Stretch / future
#### Section 12 — Advanced Features
- [x] Active learning loop: users can upvote/downvote evidence to fine-tune reranker weights and improve retrieval relevance over time (feedback API implemented).
- [ ] Integrate paid market APIs (IQVIA, Evaluate Pharma) for real market sizing data instead of synthetic TAM estimates.
- [ ] Patent claim-level semantic matching & visual claim maps to identify specific blocking claims and freedom-to-operate risks.
- [ ] Continuous crawler with change detection & data freshness scoring to auto-refresh evidence when source documents update.
- [ ] Multi-region regulatory rule engine (US/EU/India) for localized compliance checks and approval pathway analysis.
- [ ] Team collaboration features (annotations, shared notes, task assignments, @mentions) for multi-user workflows.

#### Infrastructure & Scaling
- [x] Horizontal autoscaling for agent workers (Kubernetes HPA configured for API 2-10 replicas, Celery workers 1-8 replicas).
- [x] Load balancer with health checks for multi-instance backend deployments (/health and /ready endpoints, liveness/readiness probes).
- [x] Redis/Celery task queue for async worker distribution across multiple nodes (Celery tasks with Redis broker/backend).
- [x] Distributed caching layer (Redis) for retrieval results and LLM responses (cache.py with decorators, 2hr/24hr TTL).
- [x] Database connection pooling optimization (pgBouncer configured in docker-compose.yml with transaction pooling).
- [ ] RDKit GPU acceleration path for large batch image rendering plus sandbox for untrusted SMILES inputs.
- [x] Kubernetes manifests or Helm charts for cloud-native deployment (k8s/ directory with Deployment, Service, HPA, Ingress).

#### Enterprise & Security Infrastructure
- [x] External `/api/analyze-molecule` endpoint for partners (implemented as `POST /api/jobs` with public access).
- [x] Private model hosting with strict access controls for local LLMs (Ollama with configurable models via .env).
- [x] Network egress controls enforcing default-deny (tested via `test_egress_security.py` with CI enforcement).
- [ ] Hardware-backed key storage (HSM) for high-assurance customers.
- [ ] Signed container images and supply-chain verification (in-toto, Cosign, Notary v2).
- [ ] Air-gapped deployment mode with offline model bundles and pre-indexed datasets.

#### Advanced Analytics & ML
- [x] Molecule–disease mapping model for proactive repurposing suggestions (`backend/app/ml/disease_mapper.py` with MoleculeDiseaseMapper class, category-based prediction from clinical/literature evidence, integrated into job results via `repurposing_suggestions` field).
- [x] Reranker fine-tuning pipeline using historical job feedback for retrieval optimization (`backend/app/ml/reranker_trainer.py` with gradient descent training on feedback scores, learns evidence type weights from upvotes/downvotes, CLI entry point `train_reranker_from_feedback`).
- [x] Custom embedding model training on domain-specific pharma corpus (`backend/app/ml/embedding_trainer.py` with sentence-transformers fine-tuning, contrastive learning on pharma term pairs, export for Chroma integration).

## Notes & recommended priorities
- Security is non-negotiable for pharma customers — prioritize VPC/Self-host, ZDR mode, RBAC, and audit logs early.
- RDKit integration is low-risk/high-value: land `rdkit-service` soon so UI and reports embed canonical structures; add unit tests to avoid regressions.
- CI/CD must include container and dependency scanning from day one to cover RDKit/OSRA binaries.
- Implement egress controls and PII redaction before any customer demo that touches real data.
