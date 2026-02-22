"use client";

import { useEffect, useState } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import {
  ArrowRight, CheckCircle2, FileText, Loader2,
  AlertTriangle, Clock, Activity, Terminal
} from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

const AGENTS = [
  { key: "clinical",   label: "Clinical",   color: "bg-blue-500",   textColor: "text-blue-500" },
  { key: "literature", label: "Literature", color: "bg-purple-500", textColor: "text-purple-500" },
  { key: "market",     label: "Market",     color: "bg-teal-500",   textColor: "text-teal-500" },
  { key: "patent",     label: "Patent",     color: "bg-amber-500",  textColor: "text-amber-500" },
] as const;

function JobStatusPage() {
  const searchParams = useSearchParams();
  const router       = useRouter();
  const urlJobId     = searchParams.get("id");

  const [loading, setLoading] = useState(true);
  const [error,   setError]   = useState<string | null>(null);
  const [job,     setJob]     = useState<any>(null);
  const [jobId,   setJobId]   = useState<string | null>(urlJobId);

  useEffect(() => {
    if (!urlJobId) {
      const recentJobs = jobStore.getAll();
      let target =
        recentJobs.find(j => j.status.toLowerCase() === "running") ||
        recentJobs.find(j => j.status.toLowerCase() === "pending");

      if (!target && recentJobs.length > 0) {
        target = recentJobs[0];
        if (target.status.toLowerCase() === "completed") {
          router.push(`/results?id=${target.job_id}`);
          return;
        }
      }

      if (target) {
        setJobId(target.job_id);
        router.replace(`/job?id=${target.job_id}`);
        return;
      }

      setError("No jobs found. Start a new analysis.");
      setLoading(false);
      return;
    }
    setJobId(urlJobId);
  }, [urlJobId, router]);

  useEffect(() => {
    if (!jobId) return;
    let pollTimeout: NodeJS.Timeout;
    let redirectTimeout: NodeJS.Timeout;

    const pollJob = async () => {
      try {
        const jobData = await api.getJob(jobId);
        setJob(jobData);
        const status = jobData.status.toLowerCase();

        jobStore.update(jobId, {
          status: jobData.status,
          updated_at: jobData.updated_at,
          recommendation: jobData.payload?.recommendation,
        });

        if (status === "completed") {
          redirectTimeout = setTimeout(() => router.push(`/results?id=${jobId}`), 1500);
        } else if (status === "failed") {
          setError("Job failed: " + (jobData.payload?.error || "Unknown error"));
        } else if (status === "pending" || status === "running") {
          pollTimeout = setTimeout(pollJob, 3000);
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch job");
      } finally {
        setLoading(false);
      }
    };

    pollJob();
    return () => {
      clearTimeout(pollTimeout);
      clearTimeout(redirectTimeout);
    };
  }, [jobId, router]);

  /* ── Loading state ── */
  if (loading && !job) {
    return (
      <div className="flex items-center justify-center min-h-[50vh]">
        <div className="flex items-center gap-3 text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin text-primary" />
          <span className="text-xs">Connecting to pipeline...</span>
        </div>
      </div>
    );
  }

  /* ── Error state ── */
  if (error && !job) {
    return (
      <div className="max-w-sm mx-auto mt-20">
        <div className="card-flush p-6 text-center space-y-4">
          <AlertTriangle className="h-8 w-8 text-destructive mx-auto" />
          <div>
            <p className="text-sm font-semibold">Pipeline Error</p>
            <p className="text-xs text-muted-foreground mt-1">{error}</p>
          </div>
          <div className="flex gap-2 justify-center">
            <Link href="/molecule">
              <Button size="sm" className="text-xs h-7 bg-primary text-primary-foreground">New Analysis</Button>
            </Link>
            <Link href="/history">
              <Button variant="outline" size="sm" className="text-xs h-7">History</Button>
            </Link>
          </div>
        </div>
      </div>
    );
  }

  const status    = job?.status?.toLowerCase() || "unknown";
  const workers   = job?.payload?.workers || {};
  const doneCount = Object.keys(workers).length;
  const molName   = job?.payload?.molecule || job?.molecule || "Unknown";
  const progress  =
    status === "completed" ? 100 :
    status === "running"   ? Math.min(Math.round((doneCount / 4) * 100) || 20, 90) :
    status === "pending"   ? 5 : 0;

  const statusConfig: Record<string, { label: string; color: string; spin?: boolean }> = {
    pending:   { label: "Queued",     color: "text-amber-500" },
    running:   { label: "Processing", color: "text-blue-500",   spin: true },
    completed: { label: "Complete",   color: "text-green-500" },
    failed:    { label: "Failed",     color: "text-destructive" },
  };
  const sc = statusConfig[status] || { label: status, color: "text-muted-foreground" };

  const duration = job?.created_at && job?.updated_at
    ? `${Math.round((new Date(job.updated_at).getTime() - new Date(job.created_at).getTime()) / 1000)}s`
    : "—";

  return (
    <div className="max-w-3xl mx-auto space-y-4">

      {/* Page header */}
      <div className="flex items-start justify-between gap-4">
        <div>
          <h1 className="text-base font-semibold text-foreground">{molName}</h1>
          <p className="font-mono text-[10px] text-muted-foreground mt-0.5">
            JOB/{jobId?.slice(0, 8).toUpperCase()}
          </p>
        </div>
        <div className={`flex items-center gap-1.5 text-xs font-semibold ${sc.color}`}>
          {sc.spin ? (
            <Loader2 className="h-3.5 w-3.5 animate-spin" />
          ) : status === "completed" ? (
            <CheckCircle2 className="h-3.5 w-3.5" />
          ) : status === "failed" ? (
            <AlertTriangle className="h-3.5 w-3.5" />
          ) : (
            <Clock className="h-3.5 w-3.5" />
          )}
          {sc.label}
        </div>
      </div>

      {/* Main layout */}
      <div className="grid sm:grid-cols-3 gap-4">

        {/* Left: agents + progress */}
        <div className="sm:col-span-2 space-y-3">

          {/* Progress bar card */}
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border flex items-center justify-between">
              <h2 className="text-xs font-semibold">Pipeline Progress</h2>
              <span className="font-mono text-xs font-bold text-primary">{progress}%</span>
            </div>
            <div className="p-4 space-y-3">
              <div className="h-1.5 bg-muted rounded-full overflow-hidden">
                <div
                  className={`h-full rounded-full transition-all duration-700 ease-out ${
                    status === "completed" ? "bg-green-500" :
                    status === "failed"    ? "bg-destructive" :
                    "bg-primary progress-pulse"
                  }`}
                  style={{ width: `${progress}%` }}
                />
              </div>

              {/* Agent rows */}
              <div className="space-y-0 divide-y divide-border">
                {AGENTS.map(({ key, label, color, textColor }) => {
                  const done      = !!workers[key];
                  const isRunning = status === "running" && !done;
                  const w         = workers[key];
                  const conf      = w?.confidence !== undefined
                    ? `${Math.round(w.confidence * 100)}%`
                    : null;

                  return (
                    <div key={key} className="flex items-center gap-3 py-2.5">
                      <div className={`h-1.5 w-1.5 rounded-full flex-shrink-0 ${
                        done      ? color :
                        isRunning ? `${color} animate-pulse` :
                        "bg-muted-foreground/20"
                      }`} />

                      <div className="flex-1 min-w-0">
                        <span className={`text-xs font-medium ${done ? "text-foreground" : "text-muted-foreground"}`}>
                          {label}
                        </span>
                        {w?.summary && (
                          <p className="text-[10px] text-muted-foreground truncate mt-0.5">{w.summary}</p>
                        )}
                      </div>

                      <div className="flex items-center gap-2 flex-shrink-0">
                        {conf && <span className={`font-mono text-[10px] font-bold ${textColor}`}>{conf}</span>}
                        {done ? (
                          <CheckCircle2 className={`h-3.5 w-3.5 ${textColor}`} />
                        ) : isRunning ? (
                          <Loader2 className="h-3.5 w-3.5 text-primary animate-spin" />
                        ) : (
                          <div className="h-3.5 w-3.5 rounded-full border border-border/60" />
                        )}
                      </div>
                    </div>
                  );
                })}
              </div>
            </div>
          </div>

          {/* Recommendation */}
          {job?.payload?.recommendation && (
            <div className="card-flush overflow-hidden">
              <div className="px-4 py-3 flex items-center justify-between">
                <div>
                  <p className="label-xs mb-1">Final Recommendation</p>
                  <p className="text-xs text-muted-foreground">Based on {doneCount} agent analyses</p>
                </div>
                <span className={`text-sm font-bold px-3 py-1 rounded badge-base ${
                  job.payload.recommendation === "GO" || job.payload.recommendation === "Go"
                    ? "status-go"
                    : job.payload.recommendation === "No-Go" || job.payload.recommendation === "NO-GO"
                    ? "status-no-go"
                    : "status-investigate"
                }`}>
                  {job.payload.recommendation}
                </span>
              </div>
            </div>
          )}

          {/* Actions when complete */}
          {status === "completed" && (
            <div className="flex gap-2">
              <Link href={`/results?id=${jobId}`} className="flex-1">
                <Button size="sm" className="w-full h-8 text-xs bg-primary text-primary-foreground gap-1.5">
                  <ArrowRight className="h-3 w-3" />
                  View Full Results
                </Button>
              </Link>
              <Link href={`/reports?id=${jobId}`}>
                <Button variant="outline" size="sm" className="h-8 text-xs gap-1.5">
                  <FileText className="h-3 w-3" />
                  Report
                </Button>
              </Link>
            </div>
          )}
        </div>

        {/* Right: job metadata */}
        <div className="space-y-3">
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border">
              <h2 className="text-xs font-semibold">Job Details</h2>
            </div>
            <div className="divide-y divide-border">
              {[
                { label: "Status",   value: status, mono: false, capitalize: true },
                { label: "Agents",   value: `${doneCount}/4`, mono: true },
                { label: "Progress", value: `${progress}%`, mono: true },
                { label: "Duration", value: duration, mono: true },
              ].map(({ label, value, mono, capitalize }) => (
                <div key={label} className="flex items-center justify-between px-4 py-2.5">
                  <span className="label-xs">{label}</span>
                  <span className={`text-xs font-bold ${mono ? "font-mono" : ""} ${capitalize ? "capitalize" : ""} text-foreground`}>
                    {value}
                  </span>
                </div>
              ))}
            </div>
          </div>

          {/* Raw JSON */}
          <div className="card-flush overflow-hidden">
            <details>
              <summary className="flex items-center gap-2 px-4 py-2.5 cursor-pointer hover:bg-muted/30 transition-colors border-b border-border">
                <Terminal className="h-3 w-3 text-muted-foreground" />
                <span className="text-xs font-semibold">Raw Payload</span>
              </summary>
              <div className="p-3">
                <pre className="text-[10px] font-mono text-muted-foreground bg-muted/30 rounded p-3 overflow-x-auto max-h-48 overflow-y-auto whitespace-pre-wrap">
                  {JSON.stringify(job?.payload, null, 2)}
                </pre>
              </div>
            </details>
          </div>
        </div>
      </div>
    </div>
  );
}

export default function ProtectedJobStatusPage() {
  return (
    <ProtectedRoute>
      <JobStatusPage />
    </ProtectedRoute>
  );
}
