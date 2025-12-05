import type {
  ClinicalTrial,
  CompetitorRow,
  DemoJobSnapshot,
  HistoryRun,
  IndexVersion,
  LiteratureEntry,
  MarketMetric,
  MasterPayload,
  ComparisonSlot,
  PatentSummary,
  QueueItem,
  RegulatoryNote,
  ReportSection,
  RobotsStatus,
  TimelineStep
} from "./types";

const PLACEHOLDER_STRUCTURE_SVG = `
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 120">
  <rect width="200" height="120" fill="#f5f5f5" stroke="#d0d0d0" />
  <text x="100" y="64" text-anchor="middle" font-size="14" fill="#333">Structure preview</text>
</svg>`;

export const SAMPLE_PAYLOAD: MasterPayload = {
  smiles: "O=C1NC(=O)N(C)C(=O)N1C",
  innovation_story:
    "Pirfenidone maintains a consistent anti-fibrotic signal across PF-ILD cohorts with manageable liver monitoring, opening a viable path for label expansion once follow-up registrational data matures.",
  recommendation: "Investigate",
  market_score: 72,
  report_version: 1,
  structure: {
    svg: PLACEHOLDER_STRUCTURE_SVG,
    path: "backend/app/report_assets/reports/images/structures/pirfenidone-demo.svg",
    metadata_path: "backend/app/report_assets/reports/images/metadata/pirfenidone-demo.json",
    smiles: "O=C1NC(=O)N(C)C(=O)N1C",
    source_type: "smiles",
    source_reference: "O=C1NC(=O)N(C)C(=O)N1C",
    generated_at: "2025-12-04T10:00:00Z",
    image_id: "pirfenidone-demo",
    license: "internal-use-only"
  },
  workers: {
    clinical: {
      summary:
        "Phase 2b study NCT01907900 slowed FVC decline by ~45 mL vs placebo over 24 weeks; benefit also seen in registry patients failing nintedanib.",
      confidence: 0.78,
      confidence_band: "high",
      evidence: [
        {
          type: "Trial",
          text: "Pirfenidone slowed the decline in forced vital capacity versus placebo across progressive fibrosing interstitial lung disease patients.",
          url: "https://clinicaltrials.gov/study/NCT01907900",
          confidence: 0.84,
          evidence_id: "clinical-1"
        },
        {
          type: "Registry",
          text: "Real-world compassionate use program reported stabilised FVC over 6 months in refractory IPF cases.",
          url: "https://example.org/registry",
          confidence: 0.6,
          evidence_id: "clinical-2"
        }
      ],
      metadata: {
        trials: '[{"nct_id":"NCT01907900","phase":"Phase 2","status":"Completed"}]',
        population: "PF-ILD"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 12,
        final_passages: 5,
        evidence_count: 2,
        unique_sources: 2,
        coverage_ratio: 1,
        precision_proxy: 0.5,
        high_conf_evidence: 1,
        retriever_top_k: 5
      }
    },
    literature: {
      summary: "Multiple PMC articles describe TGF-beta modulation and anti-inflammatory properties supporting use beyond IPF.",
      confidence: 0.74,
      confidence_band: "medium",
      evidence: [
        {
          type: "Mechanism",
          text: "Pirfenidone downregulates profibrotic cytokines and reduces collagen deposition in murine bleomycin models.",
          url: "https://pubmed.ncbi.nlm.nih.gov/33411160/",
          confidence: 0.72,
          evidence_id: "literature-1"
        }
      ],
      metadata: {
        sources: "3 peer-reviewed articles"
      },
      metrics: {
        query_terms: 3,
        retrieved_passages: 9,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.6,
        high_conf_evidence: 1,
        retriever_top_k: 5
      },
      alerts: [
        "Low evidence coverage (1/2)",
        "Evidence monoculture (single unique source)"
      ]
    },
    patent: {
      summary: "Key Genentech filings on pirfenidone combinations expire 2028-2030; no blocking claims on PF-ILD indication found.",
      confidence: 0.63,
      confidence_band: "medium",
      evidence: [
        {
          type: "Patent",
          text: "US20190234567 lists pirfenidone + nintedanib combos but indication claims are limited to IPF.",
          url: "https://patents.google.com/patent/US20190234567",
          confidence: 0.58,
          evidence_id: "patent-1"
        }
      ],
      metadata: {
        risk: "Low",
        regulatory_notes: "Label monitoring: hepatotoxicity"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 6,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.58,
        high_conf_evidence: 0,
        retriever_top_k: 5
      }
    },
    market: {
      summary: "~120k PF-ILD patients across US/EU with limited options; opportunity >$1.4B with modest uptake from pulmonologists.",
      confidence: 0.69,
      confidence_band: "medium",
      evidence: [
        {
          type: "Market",
          text: "IQVIA analogue sizing suggests $1.4B PF-ILD TAM by 2030 with 6% CAGR.",
          url: "https://example.org/market",
          confidence: 0.65,
          evidence_id: "market-1"
        }
      ],
      metadata: {
        score_breakdown: "Demand 0.8 | Competition 0.5 | Access 0.7"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 7,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.65,
        high_conf_evidence: 1,
        retriever_top_k: 5
      }
    }
  },
  quality: {
    status: "needs_attention",
    metrics: {
      clinical: {
        query_terms: 2,
        retrieved_passages: 12,
        final_passages: 5,
        evidence_count: 2,
        unique_sources: 2,
        coverage_ratio: 1,
        precision_proxy: 0.5,
        high_conf_evidence: 1,
        retriever_top_k: 5
      },
      literature: {
        query_terms: 3,
        retrieved_passages: 9,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.6,
        high_conf_evidence: 1,
        retriever_top_k: 5
      },
      patent: {
        query_terms: 2,
        retrieved_passages: 6,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.58,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      market: {
        query_terms: 2,
        retrieved_passages: 7,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.65,
        high_conf_evidence: 1,
        retriever_top_k: 5
      }
    },
    alerts: {
      literature: [
        "Low evidence coverage (1/2)",
        "Evidence monoculture (single unique source)"
      ]
    },
    story_checks: {
      status: "needs_attention",
      claims_total: 3,
      claims_flagged: 1,
      citation_gap_ratio: 0.333,
      flagged_claims: [
        {
          claim_id: "claim-2",
          claim_text: "Mechanism-of-action literature reinforces TGF-beta modulation for broader label expansion.",
          reason: "Low lexical overlap (0.18)",
          support_score: 0.18,
          evidence_ids: ["literature-1"]
        }
      ]
    }
  },
  validation: {
    status: "pass",
    claims_total: 3,
    claims_linked: 3,
    claim_links: [
      {
        claim_id: "claim-1",
        claim_text: "Pirfenidone maintains a consistent anti-fibrotic signal across PF-ILD cohorts with manageable liver monitoring.",
        worker: "clinical",
        evidence_ids: ["clinical-1", "clinical-2"],
        status: "linked"
      },
      {
        claim_id: "claim-2",
        claim_text: "Mechanism-of-action literature reinforces TGF-beta modulation for broader label expansion.",
        worker: "literature",
        evidence_ids: ["literature-1"],
        status: "linked"
      },
      {
        claim_id: "claim-3",
        claim_text: "Market guardrails show a $1.4B TAM with manageable competition once registrational data matures.",
        worker: "market",
        evidence_ids: ["market-1"],
        status: "linked"
      }
    ]
  }
};

export const METFORMIN_PAYLOAD: MasterPayload = {
  smiles: "CN(C)C(=N)N=C(N)N",
  innovation_story:
    "Metformin shows emerging anti-fibrotic and anti-inflammatory activity in NASH and SSc cohorts, but controlled data remains limited outside metabolic endpoints.",
  recommendation: "Monitor",
  market_score: 61,
  report_version: 2,
  structure: {
    svg: PLACEHOLDER_STRUCTURE_SVG,
    path: "backend/app/report_assets/reports/images/structures/metformin-demo.svg",
    metadata_path: "backend/app/report_assets/reports/images/metadata/metformin-demo.json",
    smiles: "CN(C)C(=N)N=C(N)N",
    source_type: "smiles",
    source_reference: "CN(C)C(=N)N=C(N)N",
    generated_at: "2025-12-02T18:45:00Z",
    image_id: "metformin-demo",
    license: "internal-use-only"
  },
  workers: {
    clinical: {
      summary:
        "Small investigator-sponsored trials report mild liver fat reduction in NASH with metformin add-on; metabolic endpoints outperform fibrosis surrogates.",
      confidence: 0.55,
      confidence_band: "medium",
      evidence: [
        {
          type: "Trial",
          text: "Metformin + lifestyle yielded a 12% relative reduction in liver fat over 36 weeks (NCT04511234).",
          url: "https://clinicaltrials.gov/study/NCT04511234",
          confidence: 0.52,
          evidence_id: "clinical-1"
        }
      ],
      metadata: {
        trials: '[{"nct_id":"NCT04511234","phase":"Phase 2","status":"Active"}]',
        population: "NASH"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 8,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.52,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      alerts: [
        "Low evidence coverage (1/2)",
        "Evidence monoculture (single unique source)"
      ]
    },
    literature: {
      summary:
        "PMC studies highlight AMPK-driven inhibition of stellate cell activation plus glycemic control benefits relevant to fibrotic settings.",
      confidence: 0.58,
      confidence_band: "medium",
      evidence: [
        {
          type: "Mechanism",
          text: "Metformin activated AMPK and reduced α-SMA expression in hepatic stellate cells in vitro.",
          url: "https://pubmed.ncbi.nlm.nih.gov/33551234/",
          confidence: 0.57,
          evidence_id: "literature-1"
        }
      ],
      metadata: {
        sources: "2 peer-reviewed articles"
      },
      metrics: {
        query_terms: 3,
        retrieved_passages: 7,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.57,
        high_conf_evidence: 1,
        retriever_top_k: 5
      }
    },
    patent: {
      summary: "Crowded generics space; combo IP around metabolic-fibrotic overlap is narrow but enforceable.",
      confidence: 0.48,
      confidence_band: "low",
      evidence: [
        {
          type: "Patent",
          text: "WO2023099981 covers metformin + GLP-1 combos targeting metabolic fibrosis but lacks organ-specific claims.",
          url: "https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2023099981",
          confidence: 0.44,
          evidence_id: "patent-1"
        }
      ],
      metadata: {
        risk: "Medium",
        regulatory_notes: "Monitor lactic acidosis risk"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 5,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.44,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      alerts: ["Low precision proxy (0.44 < 0.50)"]
    },
    market: {
      summary:
        "$0.9B upside if metabolic fibrosis segmentation succeeds; payer pushback expected without biopsy-backed endpoints.",
      confidence: 0.5,
      confidence_band: "low",
      evidence: [
        {
          type: "Market",
          text: "Consultant panel estimated $0.9B TAM assuming 15% adoption in NASH clinics by 2029.",
          url: "https://example.org/market/metformin",
          confidence: 0.48,
          evidence_id: "market-1"
        }
      ],
      metadata: {
        score_breakdown: "Demand 0.6 | Competition 0.5 | Access 0.4"
      },
      metrics: {
        query_terms: 2,
        retrieved_passages: 6,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.48,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      alerts: ["Low precision proxy (0.48 < 0.50)"]
    }
  },
  quality: {
    status: "investigate",
    metrics: {
      clinical: {
        query_terms: 2,
        retrieved_passages: 8,
        final_passages: 3,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.6,
        precision_proxy: 0.52,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      literature: {
        query_terms: 3,
        retrieved_passages: 7,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.57,
        high_conf_evidence: 1,
        retriever_top_k: 5
      },
      patent: {
        query_terms: 2,
        retrieved_passages: 5,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.44,
        high_conf_evidence: 0,
        retriever_top_k: 5
      },
      market: {
        query_terms: 2,
        retrieved_passages: 6,
        final_passages: 2,
        evidence_count: 1,
        unique_sources: 1,
        coverage_ratio: 0.4,
        precision_proxy: 0.48,
        high_conf_evidence: 0,
        retriever_top_k: 5
      }
    },
    alerts: {
      clinical: ["Low evidence coverage (1/2)", "Evidence monoculture (single unique source)"],
      patent: ["Low precision proxy (0.44 < 0.50)"],
      market: ["Low precision proxy (0.48 < 0.50)"]
    },
    story_checks: {
      status: "investigate",
      claims_total: 2,
      claims_flagged: 1,
      citation_gap_ratio: 0.5,
      flagged_claims: [
        {
          claim_id: "claim-2",
          claim_text: "Commercial upside depends on payer acceptance of metabolic-fibrosis surrogates.",
          reason: "No citations linked to claim",
          support_score: 0,
          evidence_ids: []
        }
      ]
    }
  },
  validation: {
    status: "needs_review",
    claims_total: 2,
    claims_linked: 1,
    claim_links: [
      {
        claim_id: "claim-1",
        claim_text: "Metformin shows modest anti-fibrotic signals in small NASH cohorts but lacks robust endpoints.",
        worker: "clinical",
        evidence_ids: ["clinical-1"],
        status: "linked"
      },
      {
        claim_id: "claim-2",
        claim_text: "Commercial upside depends on payer acceptance of metabolic-fibrosis surrogates.",
        worker: "market",
        evidence_ids: [],
        status: "missing"
      }
    ]
  }
};

export const DEMO_JOB: DemoJobSnapshot = {
  id: "JOB-4821",
  status: "RUNNING",
  depth: "full",
  startedAt: "2025-12-03T09:24:00Z",
  etaMinutes: 6
};

export const JOB_TIMELINE: TimelineStep[] = [
  {
    label: "Queued",
    description: "Validation + synonym expansion",
    status: "completed"
  },
  {
    label: "Workers orchestrated",
    description: "Clinical, literature, patent, market agents gathering evidence",
    status: "running"
  },
  {
    label: "Innovation story",
    description: "Master agent synthesising recommendation",
    status: "pending"
  }
];

export const JOB_LOGS: string[] = [
  "[09:24:12] Molecule normalised (Pirfenidone)",
  "[09:24:13] Synonym expansion added 6 variants",
  "[09:24:20] Clinical worker fetching trial data (CT.gov, EUCTR)",
  "[09:24:56] Clinical worker emitted 12 trial records",
  "[09:25:02] Literature worker retrieved 24 abstracts",
  "[09:25:18] Patent worker queued 4 USPTO filings",
  "[09:25:45] Market worker completed IQVIA query"
];

export const CLINICAL_TRIALS: ClinicalTrial[] = [
  {
    nctId: "NCT01907900",
    phase: "2b",
    status: "Completed",
    condition: "Progressive Fibrosing ILD",
    outcome: "+45 mL FVC vs placebo",
    population: "PF-ILD",
    sampleSize: "231",
    source: "https://clinicaltrials.gov/study/NCT01907900"
  },
  {
    nctId: "NCT05734218",
    phase: "2",
    status: "Recruiting",
    condition: "Systemic sclerosis ILD",
    outcome: "Primary: FVC slope",
    population: "Diffuse SSc",
    sampleSize: "148",
    source: "https://clinicaltrials.gov/study/NCT05734218"
  },
  {
    nctId: "NCT01262001",
    phase: "3",
    status: "Completed",
    condition: "Idiopathic pulmonary fibrosis",
    outcome: "FVC preservation",
    population: "IPF",
    sampleSize: "555",
    source: "https://www.clinicaltrialsregister.eu/ctr-search/trial/2010-023871-11"
  }
];

export const LITERATURE_EVIDENCE: LiteratureEntry[] = [
  {
    title: "Pirfenidone attenuates pro-fibrotic signalling",
    mechanism: "TGF-beta modulation + collagen suppression",
    snippet: "In murine bleomycin models, pirfenidone reduced collagen deposition by 42% and downregulated CTGF pathways.",
    strength: 0.82,
    doi: "10.1001/example.2024.1182"
  },
  {
    title: "Real-world safety of pirfenidone",
    mechanism: "Inflammation dampening in PF-ILD",
    snippet: "Compassionate use programme showed stable liver enzymes with monthly monitoring across 86 patients.",
    strength: 0.71,
    doi: "10.1016/example.2025.044"
  },
  {
    title: "Combination antifibrotic strategy",
    mechanism: "Synergy with nintedanib",
    snippet: "Dual therapy resulted in additive inhibition of fibroblast proliferation in vitro.",
    strength: 0.63,
    doi: "10.1056/example.2025.007"
  }
];

export const PATENT_SUMMARY: PatentSummary = {
  assignee: "Genentech / Roche",
  priorityDate: "2009-08-14",
  blockingClaims: [
    "Claim 1: Use of pirfenidone combinations for fibrotic lung disease",
    "Claim 5: Dosing regimen tied to liver enzyme thresholds"
  ],
  risk: "Low"
};

export const REGULATORY_NOTES: RegulatoryNote[] = [
  { label: "FDA", detail: "Monitor liver enzymes monthly for first 3 months" },
  { label: "EMA", detail: "Contraindicated with severe hepatic impairment" },
  { label: "Label Watch", detail: "Hepatotoxicity and photosensitivity remain key guardrails" }
];

export const MARKET_METRICS: MarketMetric[] = [
  { label: "Market score", value: "72 / 100", description: "Composite of demand, competition, access" },
  { label: "TAM (2030)", value: "$1.4B", description: "PF-ILD incidence x weighted adoption" },
  { label: "Incidence", value: "120k pts", description: "US + EU5" },
  { label: "Unmet need", value: "High", description: "Limited SOC beyond nintedanib" }
];

export const MARKET_COMPETITORS: CompetitorRow[] = [
  { drug: "Nintedanib", company: "BI", status: "Approved", share: "45%" },
  { drug: "Ziritaxestat", company: "Galapagos", status: "Hold", share: "-" },
  { drug: "Pamrevlumab", company: "FibroGen", status: "Phase 3", share: "-" }
];

export const HISTORY_RUNS: HistoryRun[] = [
  { molecule: "Pirfenidone", date: "Dec 03, 2025", recommendation: "Investigate", status: "Completed" },
  { molecule: "Metformin", date: "Dec 01, 2025", recommendation: "Go", status: "Completed" },
  { molecule: "Nintedanib", date: "Nov 28, 2025", recommendation: "No-Go", status: "Completed" },
  { molecule: "Lenabasum", date: "Nov 24, 2025", recommendation: "Investigate", status: "Running" }
];

export const COMPARISON_SLOTS: ComparisonSlot[] = [
  {
    jobId: DEMO_JOB.id,
    molecule: "Pirfenidone",
    lastUpdated: "Dec 03, 2025",
    payload: SAMPLE_PAYLOAD,
    reportVersion: SAMPLE_PAYLOAD.report_version
  },
  {
    jobId: "JOB-4810",
    molecule: "Metformin",
    lastUpdated: "Nov 30, 2025",
    payload: METFORMIN_PAYLOAD,
    reportVersion: METFORMIN_PAYLOAD.report_version
  }
];

export const CRAWLER_QUEUE: QueueItem[] = [
  { url: "https://clinicaltrials.gov/study/NCT01907900", status: "Done", retries: 0 },
  { url: "https://pmc.ncbi.nlm.nih.gov/article/PMC1182", status: "Running", retries: 1 },
  { url: "https://patents.google.com/patent/US20190234567", status: "Queued", retries: 0 }
];

export const ROBOTS_STATUS: RobotsStatus[] = [
  { domain: "clinicaltrials.gov", access: "Allowed", note: "Delay 1s" },
  { domain: "ncbi.nlm.nih.gov", access: "Allowed", note: "Key required" },
  { domain: "elsevier.com", access: "Blocked", note: "Paywalled" }
];

export const INDEX_VERSIONS: IndexVersion[] = [
  { version: "v1.2", date: "Dec 02, 2025", size: "44 MB", status: "Active" },
  { version: "v1.1", date: "Nov 20, 2025", size: "41 MB", status: "Archived" }
];

export const REPORT_SECTIONS: ReportSection[] = [
  {
    title: "Executive summary",
    body:
      "Pirfenidone demonstrates consistent antifibrotic benefit across PF-ILD subpopulations with an acceptable monitoring burden. Recommendation: Investigate for label expansion." 
  },
  {
    title: "Clinical evidence",
    body:
      "12 trial and registry assets contribute to the signal. Phase 2b data shows 45 mL FVC preservation vs placebo." 
  },
  {
    title: "Market viability",
    body: "TAM estimated $1.4B with moderate competition and high unmet need." 
  }
];
