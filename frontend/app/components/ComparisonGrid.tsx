"use client";

import Link from "next/link";
import type { ComparisonSlot, WorkerResultPayload } from "../types";

interface ComparisonGridProps {
  slots: ComparisonSlot[];
}

const CONFIDENCE_COPY: Record<WorkerResultPayload["confidence_band"], string> = {
  low: "Low confidence",
  medium: "Medium confidence",
  high: "High confidence"
};

const formatWorkerLabel = (key: string): string => key.charAt(0).toUpperCase() + key.slice(1);

export function ComparisonGrid({ slots }: ComparisonGridProps) {
  if (!slots.length) {
    return (
      <div className="empty-state">
        <p>No molecules selected yet.</p>
        <p className="empty-state__hint">Enter at least two job IDs to compare their evidence side by side.</p>
      </div>
    );
  }

  return (
    <div className="comparison-grid">
      {slots.map((slot) => (
        <article key={slot.jobId} className="section-card comparison-card space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <p className="eyebrow">{slot.lastUpdated}</p>
              <h2 className="text-xl font-semibold">{slot.molecule}</h2>
            </div>
            <div className="flex flex-col items-end gap-2 text-right">
              <span className="chip">{slot.payload.recommendation}</span>
              {(slot.reportVersion ?? slot.payload.report_version) && (
                <span className="badge">
                  Report V{slot.reportVersion ?? slot.payload.report_version}
                </span>
              )}
              <span className="subtle-text">Job ID · {slot.jobId}</span>
              {slot.payload.validation && (
                <span
                  className={`validation-chip validation-chip--${slot.payload.validation.status}`}
                >
                  {slot.payload.validation.claims_linked}/{slot.payload.validation.claims_total} claims linked
                </span>
              )}
            </div>
          </div>
          <p className="subtle-text">{slot.payload.innovation_story}</p>
          <div className="comparison-metrics">
            <div>
              <p className="eyebrow">Market score</p>
              <p className="text-2xl font-semibold">{slot.payload.market_score}</p>
            </div>
            <div>
              <p className="eyebrow">Recommendation</p>
              <p className="text-2xl font-semibold">{slot.payload.recommendation}</p>
            </div>
            <div>
              <p className="eyebrow">Evidence panels</p>
              <p className="text-2xl font-semibold">{Object.keys(slot.payload.workers).length}</p>
            </div>
          </div>
          <div className="comparison-workers">
            {Object.entries(slot.payload.workers).map(([workerKey, worker]) => (
              <div key={`${slot.jobId}-${workerKey}`} className="comparison-worker">
                <div className="comparison-worker__header">
                  <p className="eyebrow">{formatWorkerLabel(workerKey)}</p>
                  <span className="confidence-pill">{CONFIDENCE_COPY[worker.confidence_band]}</span>
                </div>
                <p className="comparison-worker__summary">{worker.summary}</p>
                {worker.evidence?.[0] && (
                  <p className="comparison-worker__evidence">
                    <span className="badge">Top cite</span> {worker.evidence[0].text}
                  </p>
                )}
              </div>
            ))}
          </div>
          <div className="flex flex-wrap gap-2">
            <Link className="btn-secondary" href="/reports">
              Open in report workspace
            </Link>
            <Link className="btn-secondary" href={`/job?jobId=${slot.jobId}`}>
              View job timeline
            </Link>
          </div>
        </article>
      ))}
    </div>
  );
}
