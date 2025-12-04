"use client";

import { FormEvent, useState } from "react";
import { ComparisonGrid } from "../components/ComparisonGrid";
import { COMPARISON_SLOTS } from "../sample-data";
import type { ComparisonSlot, JobApiResponse, MasterPayload, MasterPayloadWithMeta } from "../types";

const formatDate = (isoString: string): string =>
  new Date(isoString).toLocaleDateString("en-US", {
    month: "short",
    day: "2-digit",
    year: "numeric"
  });

export default function ComparisonPage() {
  const [jobIds, setJobIds] = useState<string[]>(COMPARISON_SLOTS.map((slot) => slot.jobId));
  const [slots, setSlots] = useState<ComparisonSlot[]>(COMPARISON_SLOTS);
  const [isLoading, setIsLoading] = useState(false);
  const [status, setStatus] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const handleJobIdChange = (index: number, value: string) => {
    setJobIds((prev) => prev.map((id, idx) => (idx === index ? value : id)));
  };

  const addJobField = () => {
    setJobIds((prev) => (prev.length >= 3 ? prev : [...prev, ""]));
  };

  const mapJobToSlot = (job: JobApiResponse): ComparisonSlot | null => {
    if (!job.payload) return null;
    const payloadWithMeta = job.payload as MasterPayloadWithMeta;
    const { molecule: moleculeLabel, ...rest } = payloadWithMeta;
    const normalizedPayload = rest as MasterPayload;

    return {
      jobId: job.job_id,
      molecule: moleculeLabel ?? `Job ${job.job_id}`,
      lastUpdated: formatDate(job.updated_at),
      payload: normalizedPayload
    };
  };

  const handleCompare = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    const trimmed = jobIds.map((id) => id.trim()).filter(Boolean);
    if (trimmed.length < 2) {
      setError("Enter at least two job IDs to compare.");
      setStatus(null);
      return;
    }

    try {
      setIsLoading(true);
      setError(null);
      setStatus("Fetching comparison data...");
      const params = new URLSearchParams();
      trimmed.forEach((id) => params.append("job_ids", id));
      const response = await fetch(`/api/jobs/compare?${params.toString()}`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const payload = (await response.json()) as JobApiResponse[];
      const nextSlots = payload
        .map(mapJobToSlot)
        .filter((slot): slot is ComparisonSlot => Boolean(slot));

      if (!nextSlots.length) {
        setError("Jobs have not produced comparison-ready payloads yet.");
        setStatus(null);
        return;
      }

      setSlots(nextSlots);
      setStatus("Comparison refreshed from backend jobs API.");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to load comparison data.");
      setStatus(null);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Phase 4 · Aggregation</p>
        <h1 className="text-2xl font-semibold">Multi-molecule comparison</h1>
        <p className="subtle-text">
          Line up two or three molecules to compare confidence bands, innovation stories, and per-worker evidence in a single view. Use live job IDs or stick with the demo payloads to explore the layout.
        </p>
        <form className="comparison-form" onSubmit={handleCompare}>
          <div className="comparison-form__fields">
            {jobIds.map((id, index) => (
              <label key={`job-field-${index}`} className="space-y-2">
                <span className="eyebrow">Molecule slot {index + 1}</span>
                <input
                  className="input"
                  placeholder="JOB-XXXX"
                  value={id}
                  onChange={(event) => handleJobIdChange(index, event.target.value)}
                />
              </label>
            ))}
          </div>
          <div className="comparison-form__actions">
            <button
              className="btn-secondary"
              type="button"
              onClick={addJobField}
              disabled={jobIds.length >= 3}
            >
              Add molecule slot
            </button>
            <button className="btn-primary" type="submit" disabled={isLoading}>
              {isLoading ? "Comparing..." : "Compare jobs"}
            </button>
          </div>
        </form>
        {status && <p className="text-sm text-emerald-400">{status}</p>}
        {error && <p className="text-sm text-rose-400">{error}</p>}
      </section>

      <ComparisonGrid slots={slots} />
    </div>
  );
}
