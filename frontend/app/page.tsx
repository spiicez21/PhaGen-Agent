"use client";

import { useCallback, useMemo, useState } from "react";
import type { ChangeEvent, FormEvent } from "react";
import { JobTimeline, type TimelineStep } from "./components/JobTimeline";

const API_BASE = process.env.NEXT_PUBLIC_API_URL ?? "http://localhost:8000";

interface JobResponse {
  job_id: string;
  status: "PENDING" | "RUNNING" | "COMPLETED" | "FAILED";
  recommendation?: string;
  payload?: Record<string, unknown>;
}

export default function Home() {
  const [molecule, setMolecule] = useState("Pirfenidone");
  const [job, setJob] = useState<JobResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const steps: TimelineStep[] = useMemo(() => {
    const status = job?.status ?? "PENDING";
    return [
      { label: "Queued", status: status === "PENDING" ? "running" : "completed" },
      { label: "Workers running", status: status === "RUNNING" ? "running" : status === "COMPLETED" ? "completed" : status === "FAILED" ? "failed" : "pending" },
      { label: "Report synthesized", status: status === "COMPLETED" ? "completed" : status === "FAILED" ? "failed" : "pending" }
    ];
  }, [job?.status]);

  const pollJob = useCallback(async (jobId: string) => {
    const res = await fetch(`${API_BASE}/api/jobs/${jobId}`);
    if (!res.ok) throw new Error("Failed to fetch job status");
    const payload: JobResponse = await res.json();
    setJob(payload);
    if (payload.status === "RUNNING") {
      setTimeout(() => pollJob(jobId), 1500);
    }
  }, []);

  const submit = useCallback(
    async (event: FormEvent<HTMLFormElement>) => {
      event.preventDefault();
      setLoading(true);
      setError(null);
      try {
        const res = await fetch(`${API_BASE}/api/jobs`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ molecule })
        });
        if (!res.ok) throw new Error("Unable to enqueue job");
        const payload: JobResponse = await res.json();
        setJob(payload);
        pollJob(payload.job_id);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Unknown error");
      } finally {
        setLoading(false);
      }
    },
    [molecule, pollJob]
  );

  return (
    <section className="space-y-6">
      <header className="space-y-2">
        <p className="uppercase tracking-widest text-xs text-slate-400">Repurposing engine</p>
        <h1 className="text-3xl font-semibold">PhaGen Agentic</h1>
        <p className="text-slate-300">
          Enter a molecule to trigger the master agent orchestration. The current implementation returns mocked-yet-structured worker outputs.
        </p>
      </header>

      <form className="card space-y-4" onSubmit={submit}>
        <label className="text-sm text-slate-300" htmlFor="molecule">
          Molecule or InChI Key
        </label>
        <input
          id="molecule"
          value={molecule}
          onChange={(e: ChangeEvent<HTMLInputElement>) => setMolecule(e.target.value)}
          placeholder="Ex: Pirfenidone"
        />
        <button
          type="submit"
          className="bg-emerald-500 text-slate-950 font-semibold rounded-lg px-4 py-2 disabled:opacity-50"
          disabled={loading}
        >
          {loading ? "Enqueuing..." : "Run orchestration"}
        </button>
        {error && <p className="text-rose-400 text-sm">{error}</p>}
      </form>

      <JobTimeline steps={steps} />

      {job?.payload && (
        <div className="card space-y-4">
          <div className="flex items-center justify-between">
            <h3 className="text-lg font-semibold">Recommendation</h3>
            <span className="text-emerald-400 font-bold">{job.recommendation ?? "Pending"}</span>
          </div>
          <pre className="text-xs bg-slate-900 p-4 rounded-lg overflow-auto max-h-96">
            {JSON.stringify(job.payload, null, 2)}
          </pre>
        </div>
      )}
    </section>
  );
}
