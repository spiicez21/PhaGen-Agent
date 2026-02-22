"use client";

import { useEffect, useState } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { Button } from "@/components/ui/button";
import {
  ArrowRight, BarChart3, BookOpen, FileText, Microscope, Scale,
  AlertTriangle, CheckCircle2, Info, Loader2, Download, FileJson,
  Sparkles, Shield, Activity
} from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

const exportToJSON = (data: any, filename: string) => {
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a"); link.href = url; link.download = filename; link.click();
  URL.revokeObjectURL(url);
};

const exportToCSV = (data: any, filename: string) => {
  const rows: string[][] = [["Worker", "Summary", "Confidence", "Confidence Band"]];
  Object.entries(data.workers || {}).forEach(([name, w]: [string, any]) => {
    rows.push([name, w?.summary || "", w?.confidence?.toString() || "", w?.confidence_band || ""]);
  });
  const csv = rows.map(r => r.map(c => `"${c}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a"); link.href = url; link.download = filename; link.click();
  URL.revokeObjectURL(url);
};

const WORKER_META: Record<string, { icon: any; color: string; textColor: string }> = {
  clinical:   { icon: Microscope, color: "bg-blue-500/10",   textColor: "text-blue-500"   },
  literature: { icon: BookOpen,   color: "bg-purple-500/10", textColor: "text-purple-500" },
  patent:     { icon: Scale,      color: "bg-amber-500/10",  textColor: "text-amber-500"  },
  market:     { icon: BarChart3,  color: "bg-teal-500/10",   textColor: "text-teal-500"   },
};

function ResultsPage() {
  const searchParams = useSearchParams();
  const router       = useRouter();
  const urlJobId     = searchParams.get("id");

  const [loading,   setLoading]   = useState(true);
  const [error,     setError]     = useState<string | null>(null);
  const [jobData,   setJobData]   = useState<any>(null);
  const [exporting, setExporting] = useState(false);
  const [jobId,     setJobId]     = useState<string | null>(urlJobId);

  useEffect(() => {
    if (!urlJobId) {
      const recentJobs = jobStore.getAll();
      let target = recentJobs.find(j => j.status.toLowerCase() === "completed");
      if (!target) {
        const ongoing = recentJobs.find(j => ["running", "pending"].includes(j.status.toLowerCase()));
        if (ongoing) { router.push(`/job?id=${ongoing.job_id}`); return; }
      }
      if (!target && recentJobs.length > 0) target = recentJobs[0];
      if (target) { setJobId(target.job_id); router.replace(`/results?id=${target.job_id}`); return; }
      setError("No jobs found. Please start a new analysis first.");
      setLoading(false);
      return;
    }
    setJobId(urlJobId);
  }, [urlJobId, router]);

  useEffect(() => {
    if (!jobId) return;
    const fetchJob = async () => {
      try {
        const job = await api.getJob(jobId);
        const status = job.status.toLowerCase();
        if (status === "completed" && job.payload && Object.keys(job.payload).length > 0) {
          setJobData(job.payload); setError(null);
        } else if (status === "completed") {
          setError("Job completed but results data is empty.");
        } else if (status === "failed") {
          const msg = typeof job.payload === "object" && job.payload && "error" in job.payload
            ? String((job.payload as any).error) : "Unknown error";
          setError("Job failed: " + msg);
        } else {
          setError(`Job is ${status}. Check the pipeline page.`);
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch job");
      } finally {
        setLoading(false);
      }
    };
    fetchJob();
  }, [jobId]);

  const handleExportPDF = async () => {
    if (!jobId) return;
    setExporting(true);
    try {
      const blob = await api.downloadReport(jobId);
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `${jobData?.molecule || "report"}-${jobId}.pdf`;
      link.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      alert("Failed to download PDF: " + (err instanceof Error ? err.message : "Unknown error"));
    } finally {
      setExporting(false);
    }
  };

  /* ── Loading ── */
  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-[50vh]">
        <div className="flex items-center gap-3 text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin text-primary" />
          <span className="text-xs">Loading results...</span>
        </div>
      </div>
    );
  }

  /* ── Error ── */
  if (error) {
    return (
      <div className="max-w-sm mx-auto mt-20">
        <div className="card-flush p-6 text-center space-y-4">
          <AlertTriangle className="h-8 w-8 text-destructive mx-auto" />
          <div>
            <p className="text-sm font-semibold">Results Unavailable</p>
            <p className="text-xs text-muted-foreground mt-1">{error}</p>
          </div>
          <div className="flex gap-2 justify-center">
            <Link href="/molecule"><Button size="sm" className="text-xs h-7 bg-primary text-primary-foreground">New Analysis</Button></Link>
            <Link href="/history"><Button variant="outline" size="sm" className="text-xs h-7">History</Button></Link>
          </div>
        </div>
      </div>
    );
  }

  const workers       = jobData?.workers || {};
  const validation    = jobData?.validation;
  const reportVersion = jobData?.report_version ?? 1;
  const structure     = jobData?.structure;
  const quality       = jobData?.quality;
  const moleculeName  = jobData?.molecule || "Unknown Molecule";
  const recommendation = jobData?.recommendation || "INVESTIGATE";
  const innovationStory = jobData?.innovation_story || "No innovation story available.";

  const recClass =
    recommendation === "GO" || recommendation === "Go" ? "status-go" :
    recommendation === "No-Go" || recommendation === "NO-GO" ? "status-no-go" :
    "status-investigate";

  return (
    <div className="space-y-4">

      {/* Page header */}
      <div className="flex items-start justify-between gap-4">
        <div>
          <div className="flex items-center gap-3">
            <h1 className="text-base font-semibold text-foreground">{moleculeName}</h1>
            <span className={`badge-base text-xs px-2 py-0.5 rounded ${recClass}`}>{recommendation}</span>
          </div>
          <p className="font-mono text-[10px] text-muted-foreground mt-0.5">
            JOB/{jobId?.slice(0, 8).toUpperCase()} &middot; REPORT v{reportVersion}
          </p>
        </div>

        {/* Export buttons */}
        <div className="flex items-center gap-1.5">
          <Button
            size="sm"
            onClick={handleExportPDF}
            disabled={exporting}
            className="h-7 px-3 text-xs bg-primary text-primary-foreground gap-1.5"
          >
            {exporting ? <Loader2 className="h-3 w-3 animate-spin" /> : <FileText className="h-3 w-3" />}
            PDF
          </Button>
          <Button size="sm" variant="outline" className="h-7 px-3 text-xs gap-1.5"
            onClick={() => jobData && exportToJSON(jobData, `${moleculeName}-${jobId}.json`)}>
            <FileJson className="h-3 w-3" /> JSON
          </Button>
          <Button size="sm" variant="outline" className="h-7 px-3 text-xs gap-1.5"
            onClick={() => jobData && exportToCSV(jobData, `${moleculeName}-${jobId}.csv`)}>
            <Download className="h-3 w-3" /> CSV
          </Button>
        </div>
      </div>

      {/* Main grid */}
      <div className="grid lg:grid-cols-3 gap-4">

        {/* Left: story + workers */}
        <div className="lg:col-span-2 space-y-4">

          {/* Innovation story */}
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border flex items-center gap-2">
              <Sparkles className="h-3.5 w-3.5 text-primary" />
              <h2 className="text-xs font-semibold">Innovation Story</h2>
            </div>
            <div className="px-4 py-3">
              <p className="text-xs text-muted-foreground leading-relaxed">{innovationStory}</p>
            </div>
          </div>

          {/* Worker summaries */}
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border">
              <h2 className="text-xs font-semibold">Agent Summaries</h2>
            </div>
            <div className="divide-y divide-border">
              {Object.entries(workers).map(([key, worker]: [string, any]) => {
                const meta = WORKER_META[key] || { icon: Activity, color: "bg-muted/50", textColor: "text-muted-foreground" };
                const Icon = meta.icon;
                const conf = worker?.confidence;
                const confPct = conf !== undefined ? Math.round(conf * 100) : null;
                const confColor =
                  conf >= 0.8 ? "text-green-500" :
                  conf >= 0.6 ? "text-amber-500" :
                  "text-muted-foreground";

                return (
                  <div key={key} className="px-4 py-3 flex items-start gap-3 hover:bg-muted/20 transition-colors">
                    <div className={`h-7 w-7 rounded flex items-center justify-center flex-shrink-0 mt-0.5 ${meta.color}`}>
                      <Icon className={`h-3.5 w-3.5 ${meta.textColor}`} />
                    </div>
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2 mb-1">
                        <span className="text-xs font-semibold capitalize text-foreground">{key}</span>
                        {confPct !== null && (
                          <span className={`font-mono text-[10px] font-bold ${confColor}`}>{confPct}%</span>
                        )}
                        <span className="text-[10px] text-muted-foreground/50 ml-auto">
                          {worker?.evidence?.length || 0} items
                        </span>
                      </div>
                      <p className="text-[11px] text-muted-foreground leading-relaxed line-clamp-2">
                        {worker?.summary || "No summary available"}
                      </p>
                    </div>
                  </div>
                );
              })}
            </div>
          </div>
        </div>

        {/* Right: sidebar */}
        <div className="space-y-4">

          {/* Structure */}
          {structure && (
            <div className="card-flush overflow-hidden">
              <div className="px-4 py-2.5 border-b border-border">
                <h2 className="text-xs font-semibold">Chemical Structure</h2>
              </div>
              <div className="p-3">
                {structure.svg ? (
                  <div
                    className="w-full aspect-square flex items-center justify-center bg-white dark:bg-white/95 rounded p-2 overflow-hidden [&>svg]:w-full [&>svg]:h-full [&>svg]:object-contain"
                    dangerouslySetInnerHTML={{ __html: structure.svg }}
                  />
                ) : (
                  <div className="aspect-square flex items-center justify-center bg-muted/30 rounded text-xs text-muted-foreground">
                    Preview Unavailable
                  </div>
                )}
                <div className="mt-2 font-mono text-[10px] break-all bg-muted/30 p-2 rounded text-muted-foreground select-all leading-relaxed">
                  {structure.smiles}
                </div>
                {structure.source_type && (
                  <p className="mt-1.5 text-[10px] text-muted-foreground/50 text-center uppercase tracking-wider">
                    {structure.source_type}
                  </p>
                )}
              </div>
            </div>
          )}

          {/* Retrieval quality */}
          {quality && (
            <div className="card-flush overflow-hidden">
              <div className="px-4 py-2.5 border-b border-border flex items-center justify-between">
                <h2 className="text-xs font-semibold">Retrieval Health</h2>
                <div className={`flex items-center gap-1 text-[11px] font-semibold ${
                  quality.status === "pass" ? "text-green-500" : "text-amber-500"
                }`}>
                  {quality.status === "pass"
                    ? <CheckCircle2 className="h-3 w-3" />
                    : <AlertTriangle className="h-3 w-3" />}
                  {quality.status === "pass" ? "Pass" : "Attention"}
                </div>
              </div>

              {Object.keys(quality.alerts || {}).length > 0 && (
                <div className="p-3 space-y-2 border-b border-border">
                  {Object.entries(quality.alerts).map(([worker, alerts]: [string, any]) => (
                    <div key={worker} className="rounded bg-amber-500/6 border border-amber-500/15 p-2.5">
                      <p className="text-[10px] font-bold text-amber-500 uppercase tracking-wider mb-1">{worker}</p>
                      <ul className="text-[10px] text-muted-foreground space-y-0.5">
                        {Array.isArray(alerts) && alerts.map((a: string, i: number) => <li key={i}>{a}</li>)}
                      </ul>
                    </div>
                  ))}
                </div>
              )}

              <div className="divide-y divide-border">
                {Object.entries(quality.metrics || {}).map(([worker, m]: [string, any]) => (
                  <div key={worker} className="px-4 py-2.5 flex items-center justify-between">
                    <span className="text-xs capitalize font-medium text-foreground">{worker}</span>
                    <div className="text-right">
                      <p className="font-mono text-xs font-bold text-foreground">{m?.evidence_count || 0} <span className="text-muted-foreground font-normal">items</span></p>
                      <p className="font-mono text-[10px] text-muted-foreground">
                        {m?.coverage_ratio ? (m.coverage_ratio * 100).toFixed(0) : 0}% cov
                      </p>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Evidence tabs */}
      <EvidenceTabs workers={workers} reportJobId={jobId || undefined} reportVersion={reportVersion} />

      {/* Claim traceability */}
      {validation && (
        <div className="card-flush overflow-hidden">
          <div className="px-4 py-2.5 border-b border-border flex items-center justify-between">
            <div className="flex items-center gap-2">
              <Shield className="h-3.5 w-3.5 text-primary" />
              <h2 className="text-xs font-semibold">Claim Traceability</h2>
            </div>
            <span className="font-mono text-[11px] text-muted-foreground">
              {validation.claims_linked}/{validation.claims_total} linked
            </span>
          </div>
          <div className="divide-y divide-border">
            {validation.claim_links.map((claim: any) => (
              <div key={claim.claim_id} className="px-4 py-3 space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="label-xs capitalize">{claim.worker}</span>
                  <div className={`flex items-center gap-1 text-[11px] font-semibold ${
                    claim.status === "linked" ? "text-green-500" : "text-amber-500"
                  }`}>
                    {claim.status === "linked"
                      ? <CheckCircle2 className="h-3 w-3" />
                      : <AlertTriangle className="h-3 w-3" />}
                    {claim.status === "linked" ? "Linked" : "Review"}
                  </div>
                </div>
                <p className="text-[11px] text-muted-foreground leading-relaxed">{claim.claim_text}</p>
                {claim.evidence_ids?.length > 0 && (
                  <p className="font-mono text-[10px] text-muted-foreground/50">
                    Evidence: {claim.evidence_ids.join(", ")}
                  </p>
                )}
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

export default function ProtectedResultsPage() {
  return (
    <ProtectedRoute>
      <ResultsPage />
    </ProtectedRoute>
  );
}
