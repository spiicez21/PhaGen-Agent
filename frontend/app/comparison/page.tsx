"use client";

import { FormEvent, useState } from "react";
import { ComparisonGrid } from "../components/ComparisonGrid";
import { COMPARISON_SLOTS } from "../sample-data";
import type { ComparisonSlot, JobApiResponse, MasterPayload, MasterPayloadWithMeta } from "../types";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Plus, ArrowRightLeft, Loader2, AlertCircle, CheckCircle2 } from "lucide-react";

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
    const reportVersion = job.report_version ?? payloadWithMeta.report_version;

    return {
      jobId: job.job_id,
      molecule: moleculeLabel ?? `Job ${job.job_id}`,
      lastUpdated: formatDate(job.updated_at),
      payload: normalizedPayload,
      reportVersion
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
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2 text-sm text-muted-foreground mb-2">
            <span className="px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">Phase 4</span>
            <span>Aggregation</span>
          </div>
          <CardTitle className="text-2xl">Multi-molecule Comparison</CardTitle>
          <CardDescription>
            Line up two or three molecules to compare confidence bands, innovation stories, and per-worker evidence in a single view. 
            Use live job IDs or stick with the demo payloads to explore the layout.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <form onSubmit={handleCompare} className="space-y-6">
            <div className="grid gap-6 md:grid-cols-3">
              {jobIds.map((id, index) => (
                <div key={`job-field-${index}`} className="space-y-2">
                  <Label htmlFor={`job-${index}`}>Molecule Slot {index + 1}</Label>
                  <Input
                    id={`job-${index}`}
                    placeholder="JOB-XXXX"
                    value={id}
                    onChange={(event) => handleJobIdChange(index, event.target.value)}
                  />
                </div>
              ))}
              {jobIds.length < 3 && (
                <div className="flex items-end">
                  <Button
                    type="button"
                    variant="outline"
                    className="w-full border-dashed"
                    onClick={addJobField}
                  >
                    <Plus className="mr-2 h-4 w-4" /> Add Slot
                  </Button>
                </div>
              )}
            </div>

            <div className="flex items-center justify-between pt-4 border-t">
              <div className="flex items-center gap-2 text-sm">
                {status && (
                  <span className="flex items-center text-emerald-600">
                    <CheckCircle2 className="mr-2 h-4 w-4" /> {status}
                  </span>
                )}
                {error && (
                  <span className="flex items-center text-destructive">
                    <AlertCircle className="mr-2 h-4 w-4" /> {error}
                  </span>
                )}
              </div>
              <Button type="submit" disabled={isLoading}>
                {isLoading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" /> Comparing...
                  </>
                ) : (
                  <>
                    <ArrowRightLeft className="mr-2 h-4 w-4" /> Compare Jobs
                  </>
                )}
              </Button>
            </div>
          </form>
        </CardContent>
      </Card>

      <ComparisonGrid slots={slots} />
    </div>
  );
}
