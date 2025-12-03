# PhaGen Agentic — Professional AI Pipeline Architecture
**Enterprise-Grade, End-to-End Drug Repurposing Intelligence Pipeline**

PhaGen Agentic automates early-stage molecule scouting using agentic AI, evidence-grounded RAG, domain APIs, and compliant crawling. This document defines the full data, processing, reasoning, and reporting pipeline for enterprise pharma use.

---

## 🧩 Top-Level Pipeline Overview

User Input → Data Acquisition → Normalization → Indexing → Agentic Orchestration →  
LLM Grounding → Evidence Graph → Story Engine → UI & PDF → Monitoring & Governance

---

## 0. User Input Layer

Users provide:
- Molecule name
- Synonyms, SMILES, InChI Key (optional)
- Target indication (optional)
- Priority flags (high/medium/low)

System creates:
- Job ID  
- Input metadata record  
- Execution plan  

---

## 1. Data Acquisition Layer — **API-First (Legal & Stable)**

Sources:
- ClinicalTrials.gov API  
- NCBI E-utilities (PubMed, PMC OA subset)  
- OpenFDA (drug labels, adverse events)  
- EMA public guidelines  
- WHO publications  
- OpenAlex scientific metadata  
- USPTO patent bulk XML  
- Optional: GBD, OECD datasets  

Principles:
- Prefer structured APIs over HTML.
- Perform fallback crawling only if API missing data.

---

## 2. Crawlee Robot Layer (robots-aware)

Purpose:
- Fill evidence gaps not covered by APIs.

Capabilities:
- robots.txt compliance  
- Crawl-delay + rate limiting  
- Per-host budget management  
- Polite User-Agent (`PhaGenBot/1.0`)  
- HTML capture → cleaned text → document registry  
- Store raw HTML → MinIO/S3

Do NOT crawl:
- Elsevier, Springer, Wiley, NEJM, IEEE  
- Paid datasets (IQVIA, Clarivate)

---

## 3. Normalization Pipeline (ETL for AI)

Steps:
1. HTML boilerplate stripping  
2. Text extraction  
3. PII redaction  
4. Sentence splitting  
5. Chunking (300–600 token windows)  
6. Metadata injection:  
   - URL  
   - Source type  
   - Crawl date  
   - Trust score  
7. Output as `normalized_passages.jsonl`

---

## 4. Vector Store Indexing Layer

Vector DB:
- FAISS (performance) or Chroma (simple)

Stored per passage:
- Embedding  
- Metadata  
- Source URL  
- Chunk text  
- Recency score  
- Trust score  

Index Management:
- Snapshot versioning (`index_v1`, `index_v2`, `index_v3`)  
- Incremental rebuilds  
- Deduplication  
- Embedding cache  

---

## 5. Agentic Orchestration Layer (LangGraph/CrewAI)

### Agents
#### Master Agent
- Manages workflow  
- Delegates tasks  
- Enforces timeouts  
- Validates outputs  

#### Worker Agents
- ClinicalTrials Worker  
- Literature Worker  
- Patent/IP Worker  
- Regulatory Worker  
- Market Worker  
- Safety/Toxicology Worker  
- Synonym Expansion Worker  

Each worker:
- Queries retriever  
- Gets Top-k evidence  
- Runs LLM extraction prompts  
- Outputs structured JSON  

---

## 6. LLM Reasoning & Grounding Engine

Tasks:
- Summarize retrieved data  
- Extract numeric outcomes  
- Identify signals & blockers  
- Generate grounded JSON structures  
- Calibrate confidence scores  

Pipeline:
1. Query → embeddings  
2. Semantic search + reranker  
3. Shortlisting top-k passages  
4. Context shaping / compression  
5. LLM worker prompt  
6. Output validation  

Models:
- Local: Ollama Gemma2:2B  
- Remote: GPT-4.1 / Claude 3.5  

---

## 7. Evidence Graph & Innovation Story Engine

### Evidence Graph:
- Merge clinical, literature, regulatory, market nodes  
- Cross-link citations  
- Assign weights  
- Visualizable in UI  

### Story Engine:
- Converts worker JSON →  
  - Innovation Story  
  - Patent/Risk Map  
  - Market Viability Score  
  - Recommendation (Go / Investigate / No-Go)  

Validation:
- Every claim must map to ≥1 evidence item.  
- Claims without evidence → removed.  

---

## 8. Reporting & UI Layer

Dashboard Includes:
- Molecule intake  
- Job status timeline  
- Worker execution logs  
- Evidence panels  
- Confidence badges  
- Citation viewer  

Reports:
- HTML full report  
- PDF export (Puppeteer)  
- JSON export  
- Version control (v1, v2, v3…)  

Advanced UI:
- Knowledge graph  
- Molecule comparison view  
- Citation trace explorer  

---

## 9. Monitoring, Audit & Governance

Monitoring:
- RAG precision metrics  
- Retrieval coverage  
- Worker latency  
- LLM token usage  
- Crawl status  

Audit Logs:
- Raw retrieved passages  
- Worker JSON output  
- Final story generations  
- Timestamps + models used  

Compliance:
- PII stripping  
- Evidence attribution  
- Data license enforcement  
- robots.txt logs  

Drift Detection:
- Data freshness checks  
- Weekly evaluation suite  
- Alerting on unusual agent behavior  

---

## 🎯 Summary

This pipeline provides:
- Legal & reproducible data ingestion  
- Structured AI-ready normalization  
- Strong retrieval grounding  
- Multi-agent reasoning  
- Enterprise-grade observability  
- Export-ready PDF/JSON reports  
- Compliance, auditability & governance  

It is suitable for:
- Pharma R&D teams  
- Innovation labs  
- Regulatory-compliant research environments  
- Hackathons & VC technical due diligence  

---
