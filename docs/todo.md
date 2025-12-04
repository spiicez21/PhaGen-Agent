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

## 🔜 Pending
### Phase UI — Experience & admin
- [ ] Implement stretch visualizations (knowledge graph, citation trace, molecule comparison) after MVP.
- [ ] Pagination + filtering for saved runs (status, date range, recommendation).
- [ ] Standardize confidence color palette across all evidence components.
- [ ] Worker-specific error states (timeout, missing data, robots block) in the UI.
- [ ] Visual diff viewer to compare report versions (V1 vs V2 changes).
- [ ] Molecule structure viewer component (SVG viewer + zoom + download) powered by RDKit-generated SVGs.
- [ ] Add "Request structure" action in molecule intake to create SMILES → RDKit render.

### Data & infra follow-ups
- [ ] Flesh out Postgres schema (`molecules`, `jobs`, `documents`, `passages`, `reports`) and wire persistence into backend.
- [ ] Add MinIO/S3 storage wiring for raw documents and PDF artifacts.
- [ ] Extend crawler normalization (chunking, boilerplate stripping, PII redaction) per Section 5 guidance.
- [ ] Index snapshotting cadence (daily/monthly builds) for reproducibility.
- [ ] Embedding cache to avoid re-embedding unchanged passages.
- [ ] Domain-level crawl budgets to limit pages per source.
- [ ] Evidence deduplication across clinical/literature corpora.
- [ ] Add `rdkit-service` Docker image (conda + rdkit + pillow) to `docker-compose.yml` for on-demand rendering.
- [ ] Add `smiles-normalizer` job that canonicalizes SMILES/InChI via RDKit before indexing/synonym expansion.
- [ ] Add optional OSRA image pipeline for OCR-to-SMILES conversion of scraped diagrams (only when allowed).
- [ ] Hook RDKit rendering into indexes and reports pipelines so structure images are available to workers.
- [ ] Image provenance schema storing `image_id`, `source_type` (rdkit/pubchem/scrape), `source_ref` (SMILES/URL/CID), license, and `generated_at`.

### Quality, guardrails, and ops
- [ ] Implement retrieval precision checks / coverage metrics plus guardrails (evidence thresholds, anomaly detection).
- [ ] Set up Docker Compose & GitHub Actions workflow described in Section 10 (tests + image builds).
- [ ] Prepare the hackathon demo script, slides, and fallback assets outlined in Section 11.
- [ ] Model hallucination detector or citation-gap checker for final stories.
- [ ] API call budget monitor for NCBI/OpenFDA rate limits.
- [ ] Standardize LLM temperature per worker for consistent tone.
- [ ] Build an evals suite (5–10 sample molecules with expected outputs).
- [ ] SLA-style logging: latency, token usage, and failure counts per worker.

#### ➕ Added — Security & Compliance (Core)
- [ ] Zero Data Retention (ZDR) mode flag that disables persistence of input documents/user uploads (clear memory + no S3 writes).
- [ ] VPC / On-Prem / Self-Host deployment option (Terraform + Docker Compose templates).
- [ ] Tenant isolation & RBAC with tenant-scoped DB schema/row-level permissions plus UI/API roles.
- [ ] End-to-end encryption: TLS 1.3 in transit, AES-256 at rest (DB/S3) with KMS-managed keys.
- [ ] Centralized secrets management (Vault or cloud KMS) with rotating keys.
- [ ] Immutable audit logs for worker runs, data access, and report generation (timestamps, user, job_id) with retention policy.
- [ ] PII redaction & DLP to scrub sensitive fields before storage and block exfiltration of PHI/email identifiers.
- [ ] Network egress controls enforcing default-deny for customer-hosted agents.
- [ ] Integrated SAST/DAST & dependency scanning (Snyk/Dependabot) in CI.
- [ ] Pen testing & vulnerability management cadence (quarterly pen tests, critical CVE patch SLA).
- [ ] Incident response + breach communication plan with notify list and SLA.
- [ ] Compliance artifacts: SOC2 controls, ISO27001 checklist, HIPAA/GxP mapping for customer audits.
- [ ] Data residency options plus DPA/BAA templates for enterprise deals.
- [ ] Monitoring & SIEM integration forwarding audit logs to Splunk/Datadog/Elastic.
- [ ] Periodic access reviews & enforced MFA for admin roles.
- [ ] Privacy & legal review pack (privacy FAQ, data use statement, legal review dossier).

#### ➕ Added — Security Tests & Validation
- [ ] Automated test verifying ZDR mode purges temp files and prevents S3 writes.
- [ ] CI job running static analysis plus container vulnerability scan (Trivy) on `rdkit-service` and other images.
- [ ] E2E security smoke test that attempts egress from agents; CI fails if egress succeeds in locked deployments.

### Stretch / future
- [ ] Track stretch ideas from Section 12 (active learning loop, patent semantic matching, continuous crawler, etc.) once the MVP is stable.
- [ ] Horizontal autoscaling for agent workers.
- [ ] External `/api/analyze-molecule` endpoint for partners.
- [ ] Collaboration surface (annotations, shared notes, mentions).
- [ ] Molecule–disease mapping model for proactive repurposing suggestions.
- [ ] Hardware-backed key storage (HSM) for high-assurance customers.
- [ ] Signed container images and supply-chain verification (in-toto).
- [ ] Private model hosting with strict access controls for local LLMs.
- [ ] RDKit GPU acceleration path for large batch image rendering plus sandbox for untrusted SMILES inputs.

## Notes & recommended priorities
- Security is non-negotiable for pharma customers — prioritize VPC/Self-host, ZDR mode, RBAC, and audit logs early.
- RDKit integration is low-risk/high-value: land `rdkit-service` soon so UI and reports embed canonical structures; add unit tests to avoid regressions.
- CI/CD must include container and dependency scanning from day one to cover RDKit/OSRA binaries.
- Implement egress controls and PII redaction before any customer demo that touches real data.
