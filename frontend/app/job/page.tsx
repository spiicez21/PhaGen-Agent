"use client";

import { useEffect, useState } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { JobTimeline } from "../components/JobTimeline";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ArrowRight, FileText, Terminal, Loader2, AlertTriangle, CheckCircle2 } from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

function JobStatusPage() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const urlJobId = searchParams.get("id");
  
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [job, setJob] = useState<any>(null);
  const [jobId, setJobId] = useState<string | null>(urlJobId);

  useEffect(() => {
    const initJob = async () => {
      // If no job ID in URL, try to find the most recent running/pending job
      if (!urlJobId) {
        const recentJobs = jobStore.getAll();
        console.log('[Job Status] Auto-detection: Found', recentJobs.length, 'jobs in history');
        
        // Priority 1: Look for running jobs
        let targetJob = recentJobs.find(j => j.status.toLowerCase() === "running");
        
        // Priority 2: Look for pending jobs
        if (!targetJob) {
          targetJob = recentJobs.find(j => j.status.toLowerCase() === "pending");
        }
        
        // Priority 3: Get the most recent job regardless of status
        if (!targetJob && recentJobs.length > 0) {
          targetJob = recentJobs[0]; // Most recent job (already sorted by time)
          console.log('[Job Status] No ongoing jobs, using most recent job:', targetJob.job_id, 'with status:', targetJob.status);
        }
        
        if (targetJob) {
          const status = targetJob.status.toLowerCase();
          // Check if it's completed - redirect to results instead
          if (status === "completed") {
            console.log('[Job Status] Most recent job is completed, redirecting to results:', targetJob.job_id);
            router.push(`/results?id=${targetJob.job_id}`);
            return;
          }
          // Check if it's failed - still show it for review
          else if (status === "failed") {
            console.log('[Job Status] Most recent job failed, showing details:', targetJob.job_id);
            setJobId(targetJob.job_id);
            router.replace(`/job?id=${targetJob.job_id}`);
            return;
          }
          // Ongoing job (running/pending)
          else {
            console.log('[Job Status] Auto-detected ongoing job:', targetJob.job_id);
            setJobId(targetJob.job_id);
            router.replace(`/job?id=${targetJob.job_id}`);
            return;
          }
        }
        
        // No jobs at all in history
        console.log('[Job Status] No jobs found in history');
        setError("No jobs found. Please start a new analysis or check your history.");
        setLoading(false);
        return;
      }
      
      setJobId(urlJobId);
    };

    initJob();
  }, [urlJobId, router]);

  useEffect(() => {
    if (!jobId) return;

    let pollTimeout: NodeJS.Timeout;
    let redirectTimeout: NodeJS.Timeout;
    
    const pollJob = async () => {
      try {
        const jobData = await api.getJob(jobId);
        console.log('[Job Status] Poll result:', jobData.status, 'for job', jobId);
        setJob(jobData);
        
        const status = jobData.status.toLowerCase();
        
        // Update job status in store
        jobStore.update(jobId, {
          status: jobData.status,
          updated_at: jobData.updated_at,
          recommendation: jobData.payload?.recommendation,
        });
        
        // If completed, redirect to results
        if (status === "completed") {
          console.log('[Job Status] Job completed, redirecting to results...');
          redirectTimeout = setTimeout(() => {
            router.push(`/results?id=${jobId}`);
          }, 1500);
        }
        // If failed, stop polling
        else if (status === "failed") {
          console.log('[Job Status] Job failed');
          setError("Job failed: " + (jobData.payload?.error || "Unknown error"));
        }
        // Continue polling if pending or running
        else if (status === "pending" || status === "running") {
          pollTimeout = setTimeout(() => pollJob(), 3000); // Poll every 3 seconds
        }
      } catch (err) {
        console.error('[Job Status] Polling error:', err);
        setError(err instanceof Error ? err.message : "Failed to fetch job");
      } finally {
        setLoading(false);
      }
    };

    pollJob();
    
    // Cleanup function
    return () => {
      if (pollTimeout) clearTimeout(pollTimeout);
      if (redirectTimeout) clearTimeout(redirectTimeout);
    };
  }, [jobId, router]);

  if (loading && !job) {
    return (
      <div className="flex items-center justify-center min-h-[60vh]">
        <div className="text-center space-y-4">
          <Loader2 className="h-12 w-12 animate-spin mx-auto text-primary" />
          <p className="text-muted-foreground">Loading job status...</p>
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
              Job Error
            </CardTitle>
            <CardDescription>Unable to load job information</CardDescription>
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
  const statusBadgeVariant = job?.status === "running" ? "default" : job?.status === "completed" ? "default" : "secondary";
  const moleculeName = job?.payload?.molecule || job?.molecule || "Unknown";

  const workers = job?.payload?.workers || {};
  const hasWorkerData = Object.keys(workers).length > 0;
  
  // Calculate progress metrics
  const totalWorkers = 4; // clinical, literature, market, patent
  const completedWorkers = hasWorkerData ? Object.keys(workers).length : 0;
  const progressPercentage = job?.status === "completed" ? 100 : 
                            job?.status === "running" ? Math.min(Math.round((completedWorkers / totalWorkers) * 100) || 50, 90) : 
                            job?.status === "pending" ? 5 : 0;

  return (
    <div className="space-y-6 py-8">
      {/* Header Card with Job Overview */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div className="space-y-3">
              <div className="h-1 w-8 bg-primary rounded-full" />
              <div>
                <CardTitle className="text-3xl">Job {jobId?.slice(0, 8)}</CardTitle>
                <CardDescription className="text-base mt-2">
                  {moleculeName} • Analysis Pipeline
                </CardDescription>
              </div>
            </div>
            <Badge 
              variant={statusBadgeVariant} 
              className="flex items-center gap-2 text-base px-4 py-2"
            >
              {job?.status === "running" && <Loader2 className="h-4 w-4 animate-spin" />}
              {job?.status === "completed" && <CheckCircle2 className="h-4 w-4" />}
              {job?.status === "pending" && <div className="h-2 w-2 rounded-full bg-current animate-pulse" />}
              {job?.status?.toUpperCase()}
            </Badge>
          </div>
        </CardHeader>
      </Card>

      {/* Main Content Grid */}
      <div className="grid gap-6 lg:grid-cols-3">
        {/* Left Column - Job Details */}
        <div className="lg:col-span-2 space-y-6">
          {/* Job Information Card */}
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Job Information</CardTitle>
            </CardHeader>
            <CardContent>
              <dl className="grid grid-cols-2 gap-6">
                <div className="space-y-1">
                  <dt className="text-sm text-muted-foreground">Job ID</dt>
                  <dd className="font-mono text-sm bg-muted px-2 py-1 rounded">{jobId}</dd>
                </div>
                <div className="space-y-1">
                  <dt className="text-sm text-muted-foreground">Molecule</dt>
                  <dd className="font-semibold text-lg">{moleculeName}</dd>
                </div>
                <div className="space-y-1 col-span-2">
                  <dt className="text-sm text-muted-foreground">SMILES Structure</dt>
                  <dd className="font-mono text-xs bg-muted px-3 py-2 rounded break-all">
                    {job?.smiles || job?.payload?.smiles || "Not provided"}
                  </dd>
                </div>
                <div className="space-y-1">
                  <dt className="text-sm text-muted-foreground">Status</dt>
                  <dd className="font-medium capitalize text-lg">{job?.status || "Unknown"}</dd>
                </div>
                <div className="space-y-1">
                  <dt className="text-sm text-muted-foreground">Created At</dt>
                  <dd className="text-sm">
                    {job?.created_at ? new Date(job.created_at).toLocaleString() : "-"}
                  </dd>
                </div>
                <div className="space-y-1">
                  <dt className="text-sm text-muted-foreground">Last Updated</dt>
                  <dd className="text-sm">
                    {job?.updated_at ? new Date(job.updated_at).toLocaleString() : "-"}
                  </dd>
                </div>
              </dl>
            </CardContent>
          </Card>

          {/* Worker Analysis Card - Only show if completed and has worker data */}
          {job?.status === "completed" && hasWorkerData && (
            <Card>
              <CardHeader>
                <CardTitle className="text-lg">Agent Analysis Summary</CardTitle>
                <CardDescription>Results from specialized research agents</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="space-y-4">
                  {Object.entries(workers).map(([key, worker]: [string, any]) => (
                    <div key={key} className="p-4 rounded-lg border bg-muted/30">
                      <div className="flex items-center justify-between mb-2">
                        <h4 className="font-semibold capitalize">{key} Agent</h4>
                        {worker?.confidence && (
                          <Badge variant="outline" className="font-mono">
                            {worker.confidence}% confidence
                          </Badge>
                        )}
                      </div>
                      {worker?.summary && (
                        <p className="text-sm text-muted-foreground line-clamp-2">
                          {worker.summary}
                        </p>
                      )}
                    </div>
                  ))}
                </div>
              </CardContent>
            </Card>
          )}

          {/* Action Buttons */}
          {job?.status === "completed" && (
            <Card className="border-primary/20 bg-primary/5">
              <CardContent className="pt-6">
                <div className="flex gap-3">
                  <Link href={`/results?id=${jobId}`} className="flex-1">
                    <Button className="w-full" size="lg">
                      <CheckCircle2 className="mr-2 h-5 w-5" />
                      View Complete Results
                    </Button>
                  </Link>
                  <Link href={`/reports?id=${jobId}`}>
                    <Button variant="outline" size="lg">
                      <FileText className="mr-2 h-5 w-5" />
                      Generate Report
                    </Button>
                  </Link>
                </div>
              </CardContent>
            </Card>
          )}
        </div>

        {/* Right Column - Progress & Timeline */}
        <div className="space-y-6">
          {/* Progress Card */}
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Analysis Progress</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              {job?.status === "pending" && (
                <>
                  <div className="flex items-center gap-3 p-3 rounded-lg bg-yellow-50 dark:bg-yellow-950/20 border border-yellow-200 dark:border-yellow-900">
                    <div className="h-3 w-3 rounded-full bg-yellow-500 animate-pulse" />
                    <div className="flex-1">
                      <p className="text-sm font-medium">Queued</p>
                      <p className="text-xs text-muted-foreground">Waiting for available resources</p>
                    </div>
                  </div>
                  <div className="space-y-3 pt-2">
                    <div className="flex items-center justify-between text-xs">
                      <span className="text-muted-foreground">Initializing agents...</span>
                      <span className="font-medium">{progressPercentage}%</span>
                    </div>
                    <div className="h-2 bg-muted rounded-full overflow-hidden">
                      <div 
                        className="h-full bg-yellow-500 animate-pulse transition-all duration-500" 
                        style={{ width: `${progressPercentage}%` }}
                      />
                    </div>
                  </div>
                </>
              )}

              {job?.status === "running" && (
                <>
                  <div className="flex items-center gap-3 p-3 rounded-lg bg-blue-50 dark:bg-blue-950/20 border border-blue-200 dark:border-blue-900">
                    <Loader2 className="h-5 w-5 animate-spin text-blue-500" />
                    <div className="flex-1">
                      <p className="text-sm font-medium">Processing</p>
                      <p className="text-xs text-muted-foreground">
                        {completedWorkers > 0 ? `${completedWorkers} of ${totalWorkers} agents complete` : 'Agents are analyzing the molecule'}
                      </p>
                    </div>
                  </div>
                  <div className="space-y-2">
                    {['clinical', 'literature', 'market', 'patent'].map((agent) => (
                      <div key={agent} className="flex items-center justify-between text-xs p-2 rounded bg-muted/50">
                        <span className="text-muted-foreground capitalize">{agent}</span>
                        {workers[agent] ? (
                          <CheckCircle2 className="h-3 w-3 text-green-500" />
                        ) : (
                          <Loader2 className="h-3 w-3 animate-spin text-blue-500" />
                        )}
                      </div>
                    ))}
                  </div>
                  <div className="space-y-3 pt-2">
                    <div className="flex items-center justify-between text-xs">
                      <span className="text-muted-foreground">Overall Progress</span>
                      <span className="font-medium">{progressPercentage}%</span>
                    </div>
                    <div className="h-2 bg-muted rounded-full overflow-hidden">
                      <div 
                        className="h-full bg-blue-500 animate-pulse transition-all duration-500" 
                        style={{ width: `${progressPercentage}%` }}
                      />
                    </div>
                  </div>
                </>
              )}

              {job?.status === "completed" && (
                <>
                  <div className="flex items-center gap-3 p-3 rounded-lg bg-green-50 dark:bg-green-950/20 border border-green-200 dark:border-green-900">
                    <CheckCircle2 className="h-5 w-5 text-green-500" />
                    <div className="flex-1">
                      <p className="text-sm font-medium">Completed</p>
                      <p className="text-xs text-muted-foreground">
                        {hasWorkerData ? `All ${Object.keys(workers).length} agents finished` : 'All agents finished successfully'}
                      </p>
                    </div>
                  </div>
                  <div className="space-y-2">
                    {['clinical', 'literature', 'market', 'patent'].map((agent) => (
                      <div key={agent} className="flex items-center justify-between text-xs p-2 rounded bg-muted/50">
                        <span className="text-muted-foreground capitalize">{agent}</span>
                        <CheckCircle2 className="h-3 w-3 text-green-500" />
                      </div>
                    ))}
                  </div>
                  <div className="space-y-3 pt-2">
                    <div className="flex items-center justify-between text-xs">
                      <span className="text-muted-foreground">Overall Progress</span>
                      <span className="font-medium text-green-600">100%</span>
                    </div>
                    <div className="h-2 bg-muted rounded-full overflow-hidden">
                      <div className="h-full bg-green-500 w-full transition-all duration-500" />
                    </div>
                  </div>
                  {job?.payload?.recommendation && (
                    <div className="pt-4 border-t">
                      <div className="text-center space-y-2">
                        <p className="text-xs text-muted-foreground">Final Recommendation</p>
                        <Badge 
                          variant={job.payload.recommendation === "GO" ? "default" : "destructive"}
                          className="text-lg px-4 py-1"
                        >
                          {job.payload.recommendation}
                        </Badge>
                      </div>
                    </div>
                  )}
                </>
              )}

              {job?.status === "failed" && (
                <div className="flex items-center gap-3 p-3 rounded-lg bg-red-50 dark:bg-red-950/20 border border-red-200 dark:border-red-900">
                  <AlertTriangle className="h-5 w-5 text-red-500" />
                  <div>
                    <p className="text-sm font-medium">Failed</p>
                    <p className="text-xs text-muted-foreground">Job encountered an error</p>
                  </div>
                </div>
              )}

              {/* Default/Unknown Status */}
              {!job?.status && (
                <>
                  <div className="flex items-center gap-3 p-3 rounded-lg bg-muted/50 border">
                    <Loader2 className="h-5 w-5 animate-spin text-muted-foreground" />
                    <div className="flex-1">
                      <p className="text-sm font-medium">Loading Job Data</p>
                      <p className="text-xs text-muted-foreground">Fetching job information...</p>
                    </div>
                  </div>
                  <div className="space-y-2">
                    {['Clinical', 'Literature', 'Market', 'Patent'].map((agent) => (
                      <div key={agent} className="flex items-center justify-between text-xs p-2 rounded bg-muted/30">
                        <span className="text-muted-foreground">{agent} Agent</span>
                        <div className="h-2 w-2 rounded-full bg-muted-foreground/30" />
                      </div>
                    ))}
                  </div>
                  <div className="space-y-3 pt-2">
                    <div className="flex items-center justify-between text-xs">
                      <span className="text-muted-foreground">Overall Progress</span>
                      <span className="font-medium">0%</span>
                    </div>
                    <div className="h-2 bg-muted rounded-full overflow-hidden">
                      <div className="h-full bg-muted-foreground/30 w-0" />
                    </div>
                  </div>
                </>
              )}
            </CardContent>
          </Card>

          {/* Quick Actions */}
          <Card>
            <CardHeader>
              <CardTitle className="text-lg">Quick Actions</CardTitle>
            </CardHeader>
            <CardContent className="space-y-2">
              {job?.status === "completed" && (
                <Link href={`/results?id=${jobId}`}>
                  <Button className="w-full justify-start" size="sm">
                    <CheckCircle2 className="mr-2 h-4 w-4" />
                    View Results
                  </Button>
                </Link>
              )}
              <Link href="/molecule">
                <Button variant="outline" className="w-full justify-start" size="sm">
                  <ArrowRight className="mr-2 h-4 w-4" />
                  Start New Analysis
                </Button>
              </Link>
              <Link href="/history">
                <Button variant="outline" className="w-full justify-start" size="sm">
                  <FileText className="mr-2 h-4 w-4" />
                  View All Jobs
                </Button>
              </Link>
            </CardContent>
          </Card>
        </div>
      </div>

      {/* Job Metrics & Statistics */}
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2">
            <Terminal className="h-5 w-5 text-muted-foreground" />
            <CardTitle>Job Metrics & Statistics</CardTitle>
          </div>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <div className="p-4 rounded-lg bg-muted/30 border">
              <p className="text-xs text-muted-foreground mb-1">Status</p>
              <p className="text-2xl font-bold capitalize">{job?.status || "Unknown"}</p>
            </div>
            <div className="p-4 rounded-lg bg-muted/30 border">
              <p className="text-xs text-muted-foreground mb-1">Agents</p>
              <p className="text-2xl font-bold">
                {hasWorkerData ? Object.keys(workers).length : 0}/{totalWorkers}
              </p>
            </div>
            <div className="p-4 rounded-lg bg-muted/30 border">
              <p className="text-xs text-muted-foreground mb-1">Progress</p>
              <p className="text-2xl font-bold">{progressPercentage}%</p>
            </div>
            <div className="p-4 rounded-lg bg-muted/30 border">
              <p className="text-xs text-muted-foreground mb-1">Duration</p>
              <p className="text-2xl font-bold">
                {job?.created_at && job?.updated_at 
                  ? `${Math.round((new Date(job.updated_at).getTime() - new Date(job.created_at).getTime()) / 1000)}s`
                  : "-"}
              </p>
            </div>
          </div>

          {job?.payload?.recommendation && (
            <div className="mt-4 p-4 rounded-lg bg-primary/5 border border-primary/20">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium">Recommendation</p>
                  <p className="text-xs text-muted-foreground mt-1">
                    Based on {hasWorkerData ? Object.keys(workers).length : 0} agent analyses
                  </p>
                </div>
                <Badge 
                  variant={job.payload.recommendation === "GO" ? "default" : "destructive"}
                  className="text-xl px-6 py-2"
                >
                  {job.payload.recommendation}
                </Badge>
              </div>
            </div>
          )}

          {job?.payload?.market_score !== undefined && (
            <div className="mt-4 p-4 rounded-lg bg-muted/30 border">
              <p className="text-sm font-medium mb-3">Market Score</p>
              <div className="flex items-center gap-3">
                <div className="flex-1 h-3 bg-muted rounded-full overflow-hidden">
                  <div 
                    className="h-full bg-primary transition-all duration-500"
                    style={{ width: `${job.payload.market_score}%` }}
                  />
                </div>
                <span className="text-lg font-bold">{job.payload.market_score}</span>
              </div>
            </div>
          )}

          {/* Raw JSON Toggle */}
          <details className="mt-4">
            <summary className="cursor-pointer text-sm text-muted-foreground hover:text-foreground flex items-center gap-2">
              <Terminal className="h-4 w-4" />
              View Raw JSON Data
            </summary>
            <div className="mt-2 rounded-md bg-black/50 border border-border p-4 font-mono text-xs max-h-[300px] overflow-y-auto">
              <pre className="text-muted-foreground">{JSON.stringify(job, null, 2)}</pre>
            </div>
          </details>
        </CardContent>
      </Card>
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
