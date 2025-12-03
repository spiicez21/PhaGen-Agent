export type WorkerKind = "clinical" | "literature" | "patent" | "market";

export interface WorkerEvidence {
  type: string;
  text: string;
  url: string;
  confidence: number;
}

export interface WorkerResultPayload {
  summary: string;
  evidence: WorkerEvidence[];
  confidence: number;
  confidence_band: "low" | "medium" | "high";
  metadata: Record<string, string>;
}

export interface MasterPayload {
  innovation_story: string;
  recommendation: string;
  market_score: number;
  workers: Record<string, WorkerResultPayload>;
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
}

export interface JobApiResponse {
  job_id: string;
  status: JobStatus;
  created_at: string;
  updated_at: string;
  payload?: MasterPayloadWithMeta;
  recommendation?: string;
}
