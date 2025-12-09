"use client";

import { useEffect, useState } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ArrowRight, BarChart3, BookOpen, FileText, Microscope, Scale, AlertTriangle, CheckCircle2, Info, Loader2, Download, FileJson } from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

// Export utilities
const exportToJSON = (data: any, filename: string) => {
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
};

const exportToCSV = (data: any, filename: string) => {
  // Convert workers data to CSV format
  const rows: string[][] = [];

  // Header
  rows.push(['Worker', 'Summary', 'Confidence', 'Confidence Band']);

  // Data rows
  Object.entries(data.workers || {}).forEach(([workerName, workerData]: [string, any]) => {
    rows.push([
      workerName,
      workerData?.summary || '',
      workerData?.confidence?.toString() || '',
      workerData?.confidence_band || ''
    ]);
  });

  const csvContent = rows.map(row => row.map(cell => `"${cell}"`).join(',')).join('\n');
  const blob = new Blob([csvContent], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
};

function ResultsPage() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const urlJobId = searchParams.get("id");

  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [jobData, setJobData] = useState<any>(null);
  const [exporting, setExporting] = useState(false);
  const [jobId, setJobId] = useState<string | null>(urlJobId);

  // Auto-detect completed job if no ID provided
  useEffect(() => {
    const initJob = async () => {
      if (!urlJobId) {
        const recentJobs = jobStore.getAll();
        console.log('[Results] Auto-detection: Found', recentJobs.length, 'jobs in history');

        // Priority 1: Find most recent completed job
        let targetJob = recentJobs.find(j => j.status.toLowerCase() === "completed");

        // Priority 2: If no completed jobs, check for any ongoing jobs and redirect
        if (!targetJob) {
          const ongoingJob = recentJobs.find(j => {
            const status = j.status.toLowerCase();
            return status === "running" || status === "pending";
          });
          if (ongoingJob) {
            console.log('[Results] No completed jobs, but found ongoing job. Redirecting to job status:', ongoingJob.job_id);
            router.push(`/job?id=${ongoingJob.job_id}`);
            return;
          }
        }

        // Priority 3: Fall back to most recent job regardless of status
        if (!targetJob && recentJobs.length > 0) {
          targetJob = recentJobs[0]; // Most recent job
          console.log('[Results] Using most recent job from history:', targetJob.job_id, 'with status:', targetJob.status);
        }

        if (targetJob) {
          console.log('[Results] Auto-detected job:', targetJob.job_id);
          setJobId(targetJob.job_id);
          router.replace(`/results?id=${targetJob.job_id}`);
          return;
        }

        // No jobs at all in history
        console.log('[Results] No jobs found in history');
        setError("No jobs found. Please start a new analysis first.");
        setLoading(false);
        return;
      }

      setJobId(urlJobId);
    };

    initJob();
  }, [urlJobId, router]);

  useEffect(() => {
    if (!jobId) return;

    const fetchJob = async () => {
      try {
        const job = await api.getJob(jobId);
        console.log('[Results] Fetched job:', job.job_id, 'Status:', job.status, 'Has payload:', !!job.payload);

        const status = job.status.toLowerCase();

        if (status === "completed") {
          // Check if payload exists and is not empty
          if (job.payload) {
            const payloadKeys = typeof job.payload === 'object' ? Object.keys(job.payload) : [];
            console.log('[Results] Payload keys:', payloadKeys);

            if (payloadKeys.length > 0) {
              console.log('[Results] Setting job data');
              setJobData(job.payload);
              setError(null);
            } else {
              console.log('[Results] Payload is empty');
              setError("Job completed but results data is empty. Please try refreshing or contact support.");
            }
          } else {
            console.log('[Results] No payload found');
            setError("Job completed but results are not available. Please try again or contact support.");
          }
        } else if (status === "failed") {
          const errorMsg = typeof job.payload === 'object' && job.payload && 'error' in job.payload
            ? String((job.payload as any).error)
            : "Unknown error";
          setError("Job failed: " + errorMsg);
        } else if (status === "running") {
          setError("Job is still running. Please wait for it to complete or check the job status page.");
        } else if (status === "pending") {
          setError("Job is pending and hasn't started yet. Please check the job status page.");
        } else {
          setError("Job status: " + job.status + ". Please check the job status page for details.");
        }
      } catch (err) {
        console.error('[Results] Error fetching job:', err);
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
      const link = document.createElement('a');
      link.href = url;
      link.download = `${jobData?.molecule || 'report'}-${jobId}.pdf`;
      link.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      alert('Failed to download PDF: ' + (err instanceof Error ? err.message : 'Unknown error'));
    } finally {
      setExporting(false);
    }
  };

  const handleExportJSON = () => {
    if (!jobData) return;
    const filename = `${jobData?.molecule || 'results'}-${jobId}.json`;
    exportToJSON(jobData, filename);
  };

  const handleExportCSV = () => {
    if (!jobData) return;
    const filename = `${jobData?.molecule || 'results'}-${jobId}.csv`;
    exportToCSV(jobData, filename);
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center min-h-[60vh]">
        <div className="text-center space-y-4">
          <Loader2 className="h-12 w-12 animate-spin mx-auto text-primary" />
          <p className="text-muted-foreground">Loading results...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center min-h-[60vh]">
        <Card className="max-w-md">
          <CardHeader className="space-y-3">
            <div className="h-1 w-8 bg-destructive rounded-full" />
            <CardTitle className="flex items-center gap-2 text-destructive">
              <AlertTriangle className="h-5 w-5" />
              Error Loading Results
            </CardTitle>
            <CardDescription>Unable to fetch job results</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="p-3 rounded-md bg-destructive/10 border border-destructive/20">
              <p className="text-sm text-destructive">{error}</p>
            </div>
            <div className="flex gap-3">
              <Link href="/molecule" className="flex-1">
                <Button className="w-full" size="sm">Start New Analysis</Button>
              </Link>
              <Link href="/history" className="flex-1">
                <Button variant="outline" className="w-full" size="sm">View History</Button>
              </Link>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  const workers = jobData?.workers || {};
  const validation = jobData?.validation;
  const reportVersion = jobData?.report_version ?? 1;
  const structure = jobData?.structure;
  const quality = jobData?.quality;
  const moleculeName = jobData?.molecule || "Unknown Molecule";
  const recommendation = jobData?.recommendation || "INVESTIGATE";
  const innovationStory = jobData?.innovation_story || "No innovation story available";

  const qualityStatusTone: Record<string, string> = {
    pass: "bg-emerald-100 text-emerald-800",
    needs_attention: "bg-amber-100 text-amber-800",
    investigate: "bg-rose-100 text-rose-800"
  };

  return (
    <div className="space-y-8 py-8">
      <div className="grid gap-8 md:grid-cols-3">
        <div className="md:col-span-2 space-y-8">
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <div className="space-y-3">
                  <div className="h-1 w-8 bg-primary rounded-full" />
                  <CardTitle className="text-2xl">{moleculeName} Overview</CardTitle>
                  <CardDescription className="leading-relaxed">Innovation Story & Recommendation</CardDescription>
                </div>
                <Badge variant={recommendation === "GO" ? "default" : "destructive"} className="text-lg px-4 py-1">
                  {recommendation}
                </Badge>
              </div>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="prose prose-invert max-w-none">
                <p className="text-lg leading-relaxed text-muted-foreground">
                  {innovationStory}
                </p>
              </div>

              <div className="flex flex-wrap gap-4 pt-6 mt-6 border-t">
                <Button onClick={handleExportPDF} disabled={exporting}>
                  {exporting ? (
                    <><Loader2 className="mr-2 h-4 w-4 animate-spin" /> Downloading...</>
                  ) : (
                    <><FileText className="mr-2 h-4 w-4" /> Download PDF</>
                  )}
                </Button>
                <Button variant="secondary" onClick={handleExportJSON}>
                  <FileJson className="mr-2 h-4 w-4" /> Export JSON
                </Button>
                <Button variant="outline" onClick={handleExportCSV}>
                  <Download className="mr-2 h-4 w-4" /> Export CSV
                </Button>
                <Link href="/evidence/clinical">
                  <Button variant="ghost">
                    View Evidence <ArrowRight className="ml-2 h-4 w-4" />
                  </Button>
                </Link>
              </div>
            </CardContent>
          </Card>

          <div className="grid gap-4 md:grid-cols-2">
            {Object.entries(workers).map(([key, worker]: [string, any]) => (
              <Card key={key} className="bg-secondary/5 border-secondary/20">
                <CardHeader className="pb-2">
                  <CardTitle className="text-base capitalize flex items-center gap-2">
                    {key === 'clinical' && <Microscope className="h-4 w-4" />}
                    {key === 'literature' && <BookOpen className="h-4 w-4" />}
                    {key === 'patent' && <Scale className="h-4 w-4" />}
                    {key === 'market' && <BarChart3 className="h-4 w-4" />}
                    {key} Analysis
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="text-sm text-muted-foreground line-clamp-3">
                    {worker?.summary || 'No summary available'}
                  </p>
                  <Link href={`/evidence/${key}`} className="text-xs text-primary hover:underline mt-2 inline-block">
                    View details
                  </Link>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>

        <div className="space-y-8">
          {structure && (
            <Card>
              <CardHeader>
                <CardTitle>Structure</CardTitle>
              </CardHeader>
              <CardContent className="flex flex-col items-center justify-center bg-muted/20 rounded-md p-6 space-y-4">
                {structure.svg ? (
                  <div
                    className="w-full h-64 flex items-center justify-center bg-white rounded-lg p-4 shadow-sm overflow-hidden [&>svg]:w-full [&>svg]:h-full [&>svg]:object-contain"
                    dangerouslySetInnerHTML={{ __html: structure.svg }}
                  />
                ) : (
                  <div className="text-center text-muted-foreground text-sm py-8">
                    Structure Preview Unavailable
                  </div>
                )}
                <div className="text-center w-full">
                  <div className="text-xs font-mono break-all bg-background/50 p-3 rounded border border-border/50 text-muted-foreground select-all shadow-inner">
                    {structure.smiles}
                  </div>
                  {structure.source_type && (
                    <p className="text-[10px] text-muted-foreground mt-2 font-medium tracking-wide">
                      SOURCE: {structure.source_type.toUpperCase()}
                    </p>
                  )}
                </div>
              </CardContent>
            </Card>
          )}
        </div>
      </div>

      {quality && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-1">
                <CardTitle>Retrieval Health</CardTitle>
                <CardDescription>Quality guardrails and evidence coverage</CardDescription>
              </div>
              <Badge variant={quality.status === "pass" ? "outline" : "destructive"} className={quality.status === "pass" ? "text-green-500 border-green-500" : ""}>
                {quality.status === "pass" ? "Pass" : "Needs Attention"}
              </Badge>
            </div>
          </CardHeader>
          <CardContent className="space-y-6">
            {Object.keys(quality.alerts).length > 0 ? (
              <div className="space-y-3">
                {Object.entries(quality.alerts).map(([worker, alerts]: [string, any]) => (
                  <div key={worker} className="rounded-md border border-amber-500/20 bg-amber-500/10 p-4">
                    <div className="flex items-center gap-2 text-sm font-semibold text-amber-500 mb-2">
                      <AlertTriangle className="h-4 w-4" />
                      <span className="capitalize">{worker} Alerts</span>
                    </div>
                    <ul className="list-disc list-inside text-sm text-muted-foreground">
                      {Array.isArray(alerts) && alerts.map((alert: string, index: number) => (
                        <li key={index}>{alert}</li>
                      ))}
                    </ul>
                  </div>
                ))}
              </div>
            ) : (
              <div className="rounded-md border border-green-500/20 bg-green-500/10 p-4 flex items-center gap-2 text-green-500">
                <CheckCircle2 className="h-5 w-5" />
                <span className="text-sm font-medium">All workers met minimum evidence thresholds.</span>
              </div>
            )}

            <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
              {Object.entries(quality.metrics).map(([worker, metrics]: [string, any]) => (
                <div key={worker} className="rounded-lg border bg-card p-4 space-y-3">
                  <div className="flex items-center justify-between">
                    <Badge variant="secondary" className="capitalize">{worker}</Badge>
                    <span className="text-xs text-muted-foreground">{metrics?.evidence_count || 0} items</span>
                  </div>
                  <div className="space-y-1">
                    <div className="flex justify-between text-sm">
                      <span className="text-muted-foreground">Coverage</span>
                      <span className="font-medium">{metrics?.coverage_ratio ? (metrics.coverage_ratio * 100).toFixed(0) : 0}%</span>
                    </div>
                    <div className="flex justify-between text-sm">
                      <span className="text-muted-foreground">Precision</span>
                      <span className="font-medium">{metrics?.precision_proxy ? (metrics.precision_proxy * 100).toFixed(0) : 0}%</span>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
      )}

      <div className="pt-8">
        <EvidenceTabs workers={workers} reportJobId={jobId || undefined} reportVersion={reportVersion} />
      </div>

      {validation && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-1">
                <CardTitle>Claim Traceability</CardTitle>
                <CardDescription>Validation pass for innovation story claims</CardDescription>
              </div>
              <Badge variant="outline">
                {validation.claims_linked}/{validation.claims_total} claims linked
              </Badge>
            </div>
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              {validation.claim_links.map((claim: any) => (
                <div key={claim.claim_id} className="rounded-lg border bg-card/50 p-4 space-y-2">
                  <div className="flex items-center justify-between">
                    <Badge variant="secondary" className="capitalize">{claim.worker}</Badge>
                    <Badge variant={claim.status === "linked" ? "outline" : "destructive"} className={claim.status === "linked" ? "text-green-500 border-green-500" : ""}>
                      {claim.status === "linked" ? "Linked" : "Needs Review"}
                    </Badge>
                  </div>
                  <p className="text-sm">{claim.claim_text}</p>
                  <div className="flex items-center gap-2 text-xs text-muted-foreground">
                    <Info className="h-3 w-3" />
                    Evidence IDs: {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
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
