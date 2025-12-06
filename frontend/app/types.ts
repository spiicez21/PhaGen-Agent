export type WorkerKind = "clinical" | "literature" | "patent" | "market";

export interface WorkerEvidence {
  type: string;
  text: string;
  url: string;
  confidence: number;
  evidence_id?: string;
}

export interface WorkerResultPayload {
  summary: string;
  evidence: WorkerEvidence[];
  confidence: number;
  confidence_band: "low" | "medium" | "high";
  metadata: Record<string, string>;
  metrics?: WorkerQualityMetrics;
  alerts?: string[];
}

export interface StoryClaimLink {
  claim_id: string;
  claim_text: string;
  worker: string;
  evidence_ids: string[];
  status: "linked" | "missing";
}

export interface ValidationSummary {
  status: "pass" | "needs_review";
  claims_total: number;
  claims_linked: number;
  claim_links: StoryClaimLink[];
}

export type StructureSource = "smiles" | "inchi" | "pubchem";

export interface StructureAsset {
  svg: string;
  path: string;
  smiles: string;
  metadata_path?: string;
  source_type?: StructureSource;
  source_reference?: string;
  inchi?: string;
  generated_at?: string;
  image_id?: string;
  error?: string;
  license?: string;
}

export interface MasterPayload {
  innovation_story: string;
  recommendation: string;
  market_score: number;
  workers: Record<string, WorkerResultPayload>;
  validation?: ValidationSummary;
  quality?: QualitySummary;
  report_version?: number;
  structure?: StructureAsset;
  smiles?: string;
}

export type TimelineStatus = "pending" | "running" | "completed" | "failed";

export interface TimelineStep {
  label: string;
  status: TimelineStatus;
  description?: string;
}

export type JobStatus = "PENDING" | "RUNNING" | "COMPLETED" | "FAILED";

export interface DemoJobSnapshot {
  id: string;
  status: JobStatus;
  depth: "quick" | "full";
  startedAt: string;
  etaMinutes: number;
}

export interface ClinicalTrial {
  nctId: string;
  phase: string;
  status: string;
  condition: string;
  outcome: string;
  population: string;
  sampleSize: string;
  source: string;
}

export interface LiteratureEntry {
  title: string;
  mechanism: string;
  snippet: string;
  strength: number;
  doi: string;
}

export interface PatentSummary {
  assignee: string;
  priorityDate: string;
  blockingClaims: string[];
  risk: "Low" | "Medium" | "High";
}

export interface RegulatoryNote {
  label: string;
  detail: string;
}

export interface MarketMetric {
  label: string;
  value: string;
  description: string;
}

export interface CompetitorRow {
  drug: string;
  company: string;
  status: string;
  share: string;
}

export interface HistoryRun {
  id: string;
  molecule: string;
  date: string;
  recommendation: string;
  status: string;
}

export interface QueueItem {
  url: string;
  status: string;
  retries: number;
}

export interface RobotsStatus {
  domain: string;
  access: string;
  note: string;
}

export interface IndexVersion {
  version: string;
  date: string;
  size: string;
  status: string;
}

export interface ReportSection {
  title: string;
  body: string;
}

export type MasterPayloadWithMeta = MasterPayload & { molecule?: string };

export interface ComparisonSlot {
  jobId: string;
  molecule: string;
  lastUpdated: string;
  payload: MasterPayload;
  reportVersion?: number;
}

export interface JobApiResponse {
  job_id: string;
  status: JobStatus;
  created_at: string;
  updated_at: string;
  payload?: MasterPayloadWithMeta;
  recommendation?: string;
  report_version?: number;
}

export interface WorkerQualityMetrics {
  query_terms: number;
  retrieved_passages: number;
  final_passages: number;
  evidence_count: number;
  unique_sources: number;
  coverage_ratio: number;
  precision_proxy: number;
  high_conf_evidence: number;
  retriever_top_k: number;
}

export type QualityStatus = "pass" | "needs_attention" | "investigate";

export interface StoryQualityFlag {
  claim_id: string;
  claim_text: string;
  reason: string;
  support_score: number;
  evidence_ids: string[];
}

export interface StoryQualitySummary {
  status: QualityStatus;
  claims_total: number;
  claims_flagged: number;
  citation_gap_ratio: number;
  flagged_claims: StoryQualityFlag[];
}

export type BudgetHealth = "ok" | "warning" | "exceeded";

export interface ApiBudgetSnapshot {
  label: string;
  per_minute_limit: number;
  per_day_limit: number;
  minute_used: number;
  minute_remaining: number;
  day_used: number;
  day_remaining: number;
  status: BudgetHealth;
  last_call_at?: string | null;
  reset_in_seconds: {
    minute: number;
    day: number;
  };
  lifetime_calls: number;
}

export interface QualitySummary {
  status: QualityStatus;
  metrics: Record<string, WorkerQualityMetrics>;
  alerts: Record<string, string[]>;
  story_checks?: StoryQualitySummary;
  api_budgets?: Record<string, ApiBudgetSnapshot>;
}
