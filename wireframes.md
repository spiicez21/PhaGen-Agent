# PhaGen Agentic — Complete UI Blueprint & Full Wireframes

> Purpose: This document defines the complete UI blueprint and ASCII wireframes for the **PhaGen Agentic** platform. It includes phased UI plans, feature breakdowns for every page, and detailed ASCII wireframes ready for developers, designers, or an agentic copilot.

---

# Table of Contents
1. UI Development Phases  
2. Full Page List & Feature Breakdown  
3. UI Technical Guidelines  
4. Navigation Flow  
5. UI QA Checklist  
6. Full ASCII Wireframes (All pages + Bonus visualizations)

---

# 1. UI Development Phases

## Phase 1 — Core Functional UI
Pages:
- Molecule Search  
- Job Status / Agent Pipeline  
- Results Overview (Innovation Story)

Focus: Deliver the primary molecule → analysis → result workflow.

## Phase 2 — Evidence Dashboards
Pages:
- Clinical Evidence Panel  
- Literature Evidence Panel  
- Patent & Regulatory Panel  
- Market Viability Panel

Focus: Transparency, credibility, and evidence-backed AI reasoning.

## Phase 3 — Reporting & User Workspace
Pages:
- Full Report Viewer (PDF/HTML)  
- History / Saved Reports

Focus: Professional deliverables, saved runs, repeatability.

## Phase 4 — System Intelligence & Admin
Pages:
- Crawler Status  
- Dataset / Index Manager  
- (Bonus) Molecule Knowledge Graph  
- (Bonus) Citation Trace Map  
- (Bonus) Molecule Comparison Page

---

# 2. Full Page List & Feature Breakdown

## 1. Landing / Login Page (Optional)
**Purpose:** Branding + entry point  
**Features**
- Title + tagline: *“PhaGen Agentic — Molecule Scouting in Minutes”*  
- Sign In (optional)  
- “Try Demo Molecule” button  
- Recent runs preview (optional)

## 2. Molecule Search Page (Main Input)
**Purpose:** Start new analysis  
**Core Components**
- Molecule search bar (autocomplete + synonyms)  
- Advanced inputs:  
  - SMILES  
  - InChI Key  
  - Upload (.json/.txt)  
- Analysis depth selector: Quick Scan / Full Repurposing Evaluation  
- Run Analysis button  
- Error handling (invalid molecule, missing fields)

**Outputs**
- Creates Job ID → navigates to Job Status page.

## 3. Job Status / Agent Pipeline Page
**Purpose:** Visualize multi-agent workflow  
**Core Features**
- Graph view or step timeline: Clinical, Literature, Patent, Market, Master Agent  
- Status indicators (Pending / Running / Completed / Failed)  
- Real-time logs (LLM, crawler, RAG outputs)  
- Cancel run button  
- Estimated completion time

## 4. Results Overview / Innovation Story
**Purpose:** Provide high-level insight summary  
**Components**
- Summary headline  
- 5 core evidence bullets  
- Recommendation tag: GO / INVESTIGATE / NO-GO  
- Market viability score  
- Confidence indicator  
- Buttons: View Evidence Dashboard / Download Report / View Citation Map

## 5. Clinical Evidence Panel
**Purpose:** Show trial-level signals from RAG + Crawlee  
**Components**
- Trial table: NCT ID, Phase, Status, Condition, Outcome summary  
- Expandable cards with raw text snippet, citation source, confidence  
- Filters: Phase, Status, Disease

## 6. Literature Evidence Panel
**Purpose:** Mechanistic & clinical literature intelligence  
**Components**
- Evidence list: Title, Mechanism-of-action bullets, Snippets, DOI / PMC links  
- Evidence strength meter  
- Filters: Preclinical, Mechanistic, Clinical, Meta-review

## 7. Patent & Regulatory Panel
**Purpose:** IP barriers + regulatory signals  
**Sections**
- Patent Summary: Assignee, Priority Date, Blocking Claims  
- Regulatory Notes: FDA guidance, Contraindications, Safety flags  
- Risk Indicator: Low / Medium / High  
- Expandable patent viewer with excerpt

## 8. Market Viability Panel
**Purpose:** Commercial opportunity assessment  
**Components**
- Market viability score (0–100)  
- Incidence / prevalence charts  
- Competitor map  
- Unmet need indicators  
- AI summary text

## 9. Full Report Viewer (PDF/HTML)
**Purpose:** Provide a professional deliverable  
**Components**
- Cover: Molecule name + logo  
- TOC  
- Detailed evidence sections  
- Citations list  
- Buttons: Download PDF / Export JSON / Share link (optional)

## 10. History / Saved Reports Page
**Purpose:** User workspace for past analyses  
**Components**
- Table: Molecule, Run date, Status, Recommendation  
- Search bar, Filters, Download links for each run

## 11. Crawler Status Page (Admin)
**Purpose:** Monitor ingestion pipeline  
**Components**
- Active jobs, Completed jobs, Failed jobs  
- Robots.txt status per domain  
- Crawl latency / delays  
- Retry button, Manual URL enqueue

## 12. Dataset / Index Manager (Admin)
**Purpose:** Manage FAISS/Chroma datasets  
**Components**
- List of index versions  
- Rebuild index button  
- Embedding statistics  
- Purge / clean datasets, Versioning

---

# 3. UI Technical Guidelines

**Recommended Tech Stack**
- Frontend: Next.js (React)  
- Styling: Tailwind CSS + shadcn/ui  
- Graphs: D3.js or Cytoscape.js  
- Charts: Chart.js or Recharts  
- PDF Rendering: Puppeteer / wkhtmltopdf  
- State: React Query or Zustand

**Design System**
- Use consistent spacing and typography  
- Evidence badges: colored by confidence (green/yellow/red)  
- Buttons: Primary / Secondary / Tertiary styles  
- Accessibility: keyboard navigation + screen reader labels

---

# 4. Navigation Flow

```
Landing → Molecule Search → Job Status → Results Overview
                              ↓
                     Evidence Dashboard
               (Clinical / Literature / Patent / Market)
                              ↓
                       Full Report Viewer
                              ↓
                        Saved Reports
```

---

# 5. UI QA Checklist

- [ ] All citations are clickable  
- [ ] Evidence items show source snippet  
- [ ] PDF export matches on-screen data  
- [ ] Agent pipeline updates in real-time  
- [ ] Confidence badges render correctly  
- [ ] Filters work across evidence panels

---

# 6. Full ASCII Wireframes (All pages + Bonus visualizations)

---

## 1️⃣ Landing / Login Page (Optional)

```
+-----------------------------------------------------------+
|                         PhaGen Agentic                    |
|              "Molecule Scouting in Minutes"               |
|-----------------------------------------------------------|
|                      [  Sign In Button  ]                |
|                                                           |
|                     - OR -                                |
|                                                           |
|                 [  Try Demo Molecule  ]                   |
|                                                           |
| Recent Analyses                                           |
| --------------------------------------------------------- |
| | Molecule | Date | Recommendation | [Open]             | |
| --------------------------------------------------------- |
+-----------------------------------------------------------+
```

---

## 2️⃣ Molecule Search Page

```
+-----------------------------------------------------------+
|  < Logo >                               [ Settings ]      |
|-----------------------------------------------------------|
|  Molecule Name / Synonyms                                 |
|  [______________________________________________]         |
|                                                           |
|  Advanced Inputs ▼                                         |
|  ---------------------------------------------            |
|  SMILES:  [_____________________________]                 |
|  InChI:   [_____________________________]                 |
|  Upload:  [Choose File]                                   |
|                                                           |
|  Analysis Depth:   ( ) Quick   (•) Full                   |
|                                                           |
|                  [ Run Analysis ]                         |
+-----------------------------------------------------------+
```

---

## 3️⃣ Job Status / Agent Pipeline Page

```
+------------------------------------------------------------+
| <Back>               Job #A1F234     Status: RUNNING        |
|------------------------------------------------------------|
|   Agent Pipeline                                           |
|   -------------------------------------------------------  |
|   [ Clinical Worker ] ---- Running (●●●○○)                  |
|   -------------------------------------------------------  |
|   [ Literature Worker ] ---- Pending                       |
|   -------------------------------------------------------  |
|   [ Patent Worker ] ---- Pending                           |
|   -------------------------------------------------------  |
|   [ Market Worker ] ---- Pending                           |
|   -------------------------------------------------------  |
|   [ Master Agent ] ---- Awaiting Workers                   |
|   -------------------------------------------------------  |
|                                                            |
|   Live Logs                                                |
|   -------------------------------------------------------  |
|   > Clinical Worker fetching trial data...                 |
|   > Evidence retrieved: 12 trials                          |
|   > ...                                                    |
|                                                            |
|                   [ Cancel Job ]                           |
+------------------------------------------------------------+
```

---

## 4️⃣ Results Overview / Innovation Story

```
+----------------------------------------------------------------+
| Molecule: ASPIRIN                [ Download Report ]           |
|----------------------------------------------------------------|
| Recommendation:   [ INVESTIGATE ]  | Score: 78/100             |
|----------------------------------------------------------------|
|  Innovation Story                                                 |
|----------------------------------------------------------------|
|  Aspirin shows strong mechanistic and preliminary clinical       |
|  signals for Disease X. Evidence from Phase 2 trials suggests    |
|  potential repurposing. Market viability is moderate with room   |
|  for unmet-need positioning.                                     |
|----------------------------------------------------------------|
|  Core Evidence Bullets                                           |
|  - Phase 2 trial shows 18% improvement...                        |
|  - Mechanism aligns with inflammatory pathway inhibition...      |
|  - Patent landscape low-risk post-2025...                        |
|  - Market shows moderate unmet need...                           |
|----------------------------------------------------------------|
| [ View Clinical ]   [ View Literature ]   [ View Patents ]      |
| [ View Market ]     [ Citation Map ]                            |
+----------------------------------------------------------------+
```

---

## 5️⃣ Clinical Evidence Panel

```
+---------------------------------------------------------------+
| <Back> Clinical Evidence                                       |
|---------------------------------------------------------------|
| Filters:  Phase ▼   Status ▼   Condition ▼                     |
|---------------------------------------------------------------|
| Trials                                                         |
| ------------------------------------------------------------- |
| | NCT ID     | Phase | Status   | Condition     | Outcomes  | |
| ------------------------------------------------------------- |
| | NCT01234   | 2     | Complete | Disease X     | +18%     | |
| | NCT05678   | 1     | Recruit  | Disease Y     | Preliminary | |
| ------------------------------------------------------------- |
|                                                             v |
| Expand Card:                                                  |
| ------------------------------------------------------------- |
| Title: Efficacy of Aspirin in Disease X                       |
| Sample Size: 120                                              |
| Outcome: 18% improvement                                      |
| Snippet: “Participants treated with Aspirin showed...”        |
| Source URL: https://clinicaltrials.gov/...                    |
+---------------------------------------------------------------+
```

---

## 6️⃣ Literature Evidence Panel

```
+---------------------------------------------------------------+
| <Back> Literature Evidence                                     |
|---------------------------------------------------------------|
| Filters: Type ▼ (Clinical / Preclinical / Mechanistic / Review)|
|---------------------------------------------------------------|
| Evidence List                                                  |
| ------------------------------------------------------------- |
| [Study Title]                                                 |
| Mechanism: COX inhibition reduces inflammatory cascade        |
| Snippet: “Aspirin exhibits suppression of…”                   |
| Strength: █████████░░  78%                                    |
| DOI: https://doi.org/...                                     |
| ------------------------------------------------------------- |
| [Another Study]                                               |
| Mechanism: Platelet modulation...                             |
| Snippet: “The study demonstrated…”                            |
| Strength: ██████░░░░  60%                                     |
| ------------------------------------------------------------- |
+---------------------------------------------------------------+
```

---

## 7️⃣ Patent & Regulatory Panel

```
+---------------------------------------------------------------+
| <Back> Patent & Regulatory                                     |
|---------------------------------------------------------------|
| Patent Summary                                                 |
| ------------------------------------------------------------- |
| Assignee: Bayer AG                                            |
| Priority Date: 2005                                            |
| Blocking Claims:                                               |
|  - Claim 1: Use in treating Disease X                         |
|  - Claim 5: Specific dosing regimen                           |
| Risk Level: [ LOW ]                                           |
| ------------------------------------------------------------- |
| Regulatory Notes                                               |
| ------------------------------------------------------------- |
| - FDA warning for GI bleeding risk                            |
| - Contraindications: XYZ                                      |
| - Label Source: link                                          |
+---------------------------------------------------------------+
```

---

## 8️⃣ Market Viability Panel

```
+---------------------------------------------------------------+
| <Back> Market Viability                                       |
|---------------------------------------------------------------|
| Market Score: 78/100                                          |
|---------------------------------------------------------------|
| Incidence Chart (ASCII approximation)                         |
| ████████████░░░░                                               |
| Condition Incidence: 14M/year                                 |
|---------------------------------------------------------------|
| Competitor Landscape                                          |
| ------------------------------------------------------------- |
| | Drug      | Company | Status | Market Share               | |
| ------------------------------------------------------------- |
| | Drug A    | ABC Co. | Approved | 35%                      | |
| | Drug B    | XYZ Co. | Phase 2  | -                        | |
| ------------------------------------------------------------- |
| Market Notes                                                  |
| - Moderate unmet need                                         |
| - Growing incidence trends                                    |
+---------------------------------------------------------------+
```

---

## 9️⃣ Full Report Viewer

```
+----------------------------------------------------------------+
| <Back> Full Report — ASPIRIN                                   |
|----------------------------------------------------------------|
| [ Download PDF ] [ Export JSON ]                               |
|----------------------------------------------------------------|
| TABLE OF CONTENTS                                              |
| - Summary                                                       |
| - Clinical Evidence                                             |
| - Literature Analysis                                           |
| - Patent & Regulatory                                           |
| - Market Viability                                              |
| - Citations                                                     |
|----------------------------------------------------------------|
| Report Body (Rendered HTML/PDF)                                |
|----------------------------------------------------------------|
| ... full narrative ...                                          |
+----------------------------------------------------------------+
```

---

## 🔟 Saved Reports / History Page

```
+--------------------------------------------------------------+
| Saved Reports                                                |
|--------------------------------------------------------------|
| Search: [________________________]                            |
|--------------------------------------------------------------|
| | Molecule | Date       | Recommendation | Actions         | |
|--------------------------------------------------------------|
| | Aspirin  | Jan 4, 25  | Investigate    | [Open] [PDF]    | |
| | Metformin| Jan 3, 25  | Go             | [Open] [PDF]    | |
+--------------------------------------------------------------+
```

---

## 1️⃣1️⃣ Crawler Status Page

```
+--------------------------------------------------------------+
| Crawler Status                                                |
|--------------------------------------------------------------|
| Active Jobs: 4     Completed: 320     Failed: 2              |
|--------------------------------------------------------------|
| Queue Status                                                 |
| ------------------------------------------------------------ |
| | URL                                  | Status   | Retries | |
| ------------------------------------------------------------ |
| | https://clinicaltrials.gov/...       | DONE     | 0      | |
| | https://pmc/articles/...             | RUNNING  | 1      | |
| ------------------------------------------------------------ |
| Robots.txt Summary                                             |
| - clinicaltrials.gov: Allowed (Delay: 1s)                    |
| - ncbi.nlm.nih.gov: Allowed                                  |
| - elsevier.com: BLOCKED                                      |
|--------------------------------------------------------------|
| [ Add URL ] [ Retry Failed ]                                 |
+--------------------------------------------------------------+
```

---

## 1️⃣2️⃣ Dataset / Index Manager

```
+--------------------------------------------------------------+
| Dataset / Index Manager                                       |
|--------------------------------------------------------------|
| Index Versions                                                |
| ------------------------------------------------------------ |
| | Version | Date       | Size | Status     | Actions       | |
| ------------------------------------------------------------ |
| | v1.0    | Jan 3,25   | 32MB | Active     | [Rebuild]     | |
| | v0.9    | Jan 2,25   | 29MB | Archived   | [Restore]     | |
| ------------------------------------------------------------ |
| Embedding Stats                                               |
| - Documents: 4,200                                            |
| - Passages: 15,300                                            |
| - Avg token per passage: 420                                  |
+--------------------------------------------------------------+
```

---

## ⭐ 13️⃣ Molecule Knowledge Graph

```
+--------------------------------------------------------------+
| Molecule Knowledge Graph                                      |
|--------------------------------------------------------------|
|        (Molecule)                                             |
|            ●                                                  |
|         /     \                                               |
|       ●         ●                                             |
|  (Target A)   (Disease X)                                     |
|      |             \                                          |
|      |              ● (Trial)                                 |
|      |                                                     ●  |
|      |-----------------------------------------------● (Patent)|
|                                                                |
| Click any node → Right panel shows evidence & citations        |
+--------------------------------------------------------------+
```

---

## ⭐ 14️⃣ Citation Trace Map

```
+--------------------------------------------------------------+
| Citation Trace Map                                            |
|--------------------------------------------------------------|
| Claim → Worker → Passage → Source URL                         |
|--------------------------------------------------------------|
| [Claim] ----> [Clinical Agent] ----> [Snippet #12] ----> URL  |
| [Claim] ----> [Literature Agent] --> [Snippet #07] ----> URL |
| [Claim] ----> [Patent Agent] ------> [Patent Claim 4] ---> URL|
|--------------------------------------------------------------|
| Click any block → expands to right panel                      |
+--------------------------------------------------------------+
```

---

## ⭐ 15️⃣ Molecule Comparison Page

```
+---------------------------------------------------------------+
| Compare Molecules                                              |
|---------------------------------------------------------------|
| Molecule A: [ Aspirin ▼ ]   vs   Molecule B: [ Ibuprofen ▼ ]  |
|---------------------------------------------------------------|
| Clinical Evidence Strength                                     |
| Aspirin   ██████████░░  80%                                    |
| Ibuprofen ███████░░░░  60%                                     |
|---------------------------------------------------------------|
| Mechanism Comparison                                           |
| COX1 inhibition vs COX1/COX2                                  |
|---------------------------------------------------------------|
| Patent Risk                                                    |
| Aspirin: Low   |   Ibuprofen: Medium                           |
|---------------------------------------------------------------|
| Market Viability                                               |
| Aspirin: 78/100 | Ibuprofen: 65/100                            |
+---------------------------------------------------------------+
```

---

# End of Document
