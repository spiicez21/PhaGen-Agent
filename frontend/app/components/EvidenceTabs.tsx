"use client";

import { useEffect, useMemo, useState } from "react";
import type { WorkerKind, WorkerResultPayload } from "../types";

const PANEL_META: Record<WorkerKind, { label: string; description: string; empty: string }> = {
  clinical: {
    label: "Clinical",
    description: "Signals from trials, registries, and observational programs.",
    empty: "No clinical snippets yet. Run the crawler + index builder, then re-run this molecule."
  },
  literature: {
    label: "Literature",
    description: "Mechanism-of-action insights across preclinical and clinical papers.",
    empty: "Literature worker has no evidence for this run."
  },
  patent: {
    label: "Patent & Regulatory",
    description: "Freedom-to-operate, priority dates, and regulatory guardrails.",
    empty: "Patent worker metadata is empty. Add USPTO/FDA corpora to populate it."
  },
  market: {
    label: "Market",
    description: "Commercial sizing, competitor snapshots, and unmet-need commentary.",
    empty: "Market worker requires incidence / prevalence data to render cards."
  }
};

interface EvidenceTabsProps {
  workers?: Record<string, WorkerResultPayload>;
}

const formatMetadataValue = (value: string): string => {
  const trimmed = value?.trim?.() ?? "";
  if (!trimmed) return "--";
  if ((trimmed.startsWith("[") && trimmed.endsWith("]")) || (trimmed.startsWith("{") && trimmed.endsWith("}"))) {
    try {
      const parsed = JSON.parse(trimmed);
      if (Array.isArray(parsed)) return `${parsed.length} entries`;
      if (typeof parsed === "object" && parsed !== null) return `${Object.keys(parsed).length} fields`;
    } catch (err) {
      return value;
    }
  }
  return value;
};

const confidenceLabel = (score?: number) => {
  if (score === undefined) return "";
  if (score >= 0.8) return "High confidence";
  if (score >= 0.6) return "Moderate confidence";
  return "Signal needs validation";
};

export const EvidenceTabs = ({ workers = {} }: EvidenceTabsProps) => {
  const workerKeys = useMemo(() => (Object.keys(PANEL_META) as WorkerKind[]), []);
  const available = useMemo(
    () => workerKeys.filter((key) => Boolean(workers[key])),
    [workerKeys, workers]
  );
  const [active, setActive] = useState<WorkerKind>(available[0] ?? "clinical");

  useEffect(() => {
    if (available.length && !available.includes(active)) {
      setActive(available[0]);
    }
  }, [active, available]);

  const meta = PANEL_META[active];
  const worker = workers[active];

  return (
    <section className="panel" aria-labelledby="evidence-heading">
      <header className="panel__header">
        <div>
          <p className="eyebrow">Evidence dashboard</p>
          <h2 id="evidence-heading" className="text-2xl font-semibold">{meta.label}</h2>
          <p className="panel__description">{meta.description}</p>
        </div>
        <div className="tablist" role="tablist">
          {workerKeys.map((key) => (
            <button
              key={key}
              type="button"
              className={`tab ${active === key ? "tab--active" : ""}`}
              aria-pressed={active === key}
              onClick={() => setActive(key)}
            >
              {PANEL_META[key].label}
            </button>
          ))}
        </div>
      </header>

      {!worker ? (
        <div className="empty-state">
          <p>{meta.empty}</p>
          <p className="empty-state__hint">Tip: ingest fresh data then rerun the agent orchestration.</p>
        </div>
      ) : (
        <div className="evidence-layout">
          <div className="summary-card">
            <p className="eyebrow">Worker summary</p>
            <p className="summary-card__text">{worker.summary}</p>
            <p className="summary-card__confidence">{confidenceLabel(worker.confidence)}</p>
            {Object.keys(worker.metadata ?? {}).length > 0 && (
              <dl className="metadata-grid">
                {Object.entries(worker.metadata).map(([key, value]) => (
                  <div key={key}>
                    <dt>{key.replaceAll("_", " ")}</dt>
                    <dd>{formatMetadataValue(value)}</dd>
                  </div>
                ))}
              </dl>
            )}
          </div>
          <div className="evidence-stack">
            {worker.evidence?.length ? (
              worker.evidence.map((item, idx) => (
                <article key={`${active}-${idx}`} className="evidence-card">
                  <header>
                    <span className="badge">{item.type}</span>
                    <span className="confidence-pill">{Math.round(item.confidence * 100)}% confidence</span>
                  </header>
                  <p>{item.text}</p>
                  {item.url && (
                    <a href={item.url} target="_blank" rel="noreferrer" className="evidence-link">
                      View source ↗
                    </a>
                  )}
                </article>
              ))
            ) : (
              <div className="empty-state">
                <p>No evidence bullets yet for this worker.</p>
              </div>
            )}
          </div>
        </div>
      )}
    </section>
  );
};
