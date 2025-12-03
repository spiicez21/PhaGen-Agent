"use client";

import { useState } from "react";
import { DEMO_JOB, REPORT_SECTIONS, SAMPLE_PAYLOAD } from "../sample-data";

const downloadBlob = (blob: Blob, filename: string) => {
  const url = window.URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  window.URL.revokeObjectURL(url);
};

export default function ReportsPage() {
  const [jobId, setJobId] = useState<string>(DEMO_JOB.id);
  const [isDownloading, setIsDownloading] = useState(false);
  const [isExportingJson, setIsExportingJson] = useState(false);
  const [status, setStatus] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const normalizedJobId = jobId.trim();

  const requireJobId = (): string | null => {
    if (!normalizedJobId) {
      setError("Enter a job ID before exporting.");
      return null;
    }
    setError(null);
    return normalizedJobId;
  };

  const handlePdfDownload = async () => {
    const id = requireJobId();
    if (!id) return;
    try {
      setIsDownloading(true);
      setStatus("Preparing PDF...");
      const response = await fetch(`/api/jobs/${id}/report.pdf`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const pdfBlob = await response.blob();
      downloadBlob(pdfBlob, `phagen-report-${id}.pdf`);
      setStatus("PDF downloaded successfully.");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to download PDF.");
      setStatus(null);
    } finally {
      setIsDownloading(false);
    }
  };

  const handleJsonExport = async () => {
    const id = requireJobId();
    if (!id) return;
    try {
      setIsExportingJson(true);
      setStatus("Fetching run payload...");
      const response = await fetch(`/api/jobs/${id}`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const payload = await response.json();
      const blob = new Blob([JSON.stringify(payload, null, 2)], {
        type: "application/json",
      });
      downloadBlob(blob, `phagen-run-${id}.json`);
      setStatus("JSON exported successfully.");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to export JSON.");
      setStatus(null);
    } finally {
      setIsExportingJson(false);
    }
  };

  const validation = SAMPLE_PAYLOAD.validation;

  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Report workspace</p>
        <h1 className="text-2xl font-semibold">Full deliverable layout</h1>
        <p className="subtle-text">
          Review, edit, and export the structured report. PDF + JSON export hooks connect directly to the backend jobs API.
        </p>
        <div className="grid gap-3 md:grid-cols-[2fr_auto_auto] md:items-end">
          <label className="space-y-2">
            <span className="eyebrow">Job ID</span>
            <input
              className="input"
              placeholder="JOB-XXXX"
              value={jobId}
              onChange={(event) => setJobId(event.target.value)}
            />
          </label>
          <button
            className="btn-primary"
            type="button"
            onClick={handlePdfDownload}
            disabled={isDownloading}
          >
            {isDownloading ? "Generating..." : "Download PDF"}
          </button>
          <button
            className="btn-secondary"
            type="button"
            onClick={handleJsonExport}
            disabled={isExportingJson}
          >
            {isExportingJson ? "Exporting..." : "Export JSON"}
          </button>
        </div>
        {status && <p className="text-sm text-emerald-400">{status}</p>}
        {error && <p className="text-sm text-rose-400">{error}</p>}
      </section>

      <section className="section-card space-y-4">
        <p className="eyebrow">Table of contents</p>
        <ol className="space-y-2 list-decimal pl-6">
          {REPORT_SECTIONS.map((section) => (
            <li key={section.title}>{section.title}</li>
          ))}
        </ol>
      </section>

      {validation && (
        <section className="section-card space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <p className="eyebrow">Validation snapshot</p>
              <h2 className="text-xl font-semibold">Claims ↔ Evidence</h2>
            </div>
            <span className={`validation-chip validation-chip--${validation.status}`}>
              {validation.claims_linked}/{validation.claims_total} linked
            </span>
          </div>
          <div className="validation-claims">
            {validation.claim_links.map((claim) => (
              <article key={claim.claim_id} className="validation-claim space-y-2">
                <p className="eyebrow">{claim.worker}</p>
                <p>{claim.claim_text}</p>
                <p className="subtle-text">
                  Evidence IDs · {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                </p>
              </article>
            ))}
          </div>
        </section>
      )}

      {REPORT_SECTIONS.map((section) => (
        <section key={section.title} className="section-card space-y-2">
          <h2 className="text-xl font-semibold">{section.title}</h2>
          <p>{section.body}</p>
        </section>
      ))}
    </div>
  );
}
