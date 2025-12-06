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
  const jobId = searchParams.get("id");
  
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [job, setJob] = useState<any>(null);

  useEffect(() => {
    if (!jobId) {
      setError("No job ID provided");
      setLoading(false);
      return;
    }

    const pollJob = async () => {
      try {
        const jobData = await api.getJob(jobId);
        setJob(jobData);
        
        // Update job status in store
        jobStore.update(jobId, {
          status: jobData.status,
          updated_at: jobData.updated_at,
          recommendation: jobData.payload?.recommendation,
        });
        
        // If completed, redirect to results
        if (jobData.status === "completed") {
          setTimeout(() => {
            router.push(`/results?id=${jobId}`);
          }, 1000);
        }
        // If failed, stop polling
        else if (jobData.status === "failed") {
          setError("Job failed: " + (jobData.payload?.error || "Unknown error"));
        }
        // Continue polling if pending or running
        else if (jobData.status === "pending" || jobData.status === "running") {
          setTimeout(() => pollJob(), 3000); // Poll every 3 seconds
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch job");
      } finally {
        setLoading(false);
      }
    };

    pollJob();
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
          <CardHeader>
            <CardTitle className="flex items-center gap-2 text-destructive">
              <AlertTriangle className="h-5 w-5" />
              Job Error
            </CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-muted-foreground">{error}</p>
            <div className="mt-4 space-x-2">
              <Link href="/molecule">
                <Button>Start New Analysis</Button>
              </Link>
              <Link href="/history">
                <Button variant="outline">View History</Button>
              </Link>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }
  const statusBadgeVariant = job?.status === "running" ? "default" : job?.status === "completed" ? "default" : "secondary";
  const moleculeName = job?.payload?.molecule || "Unknown";

  return (
    <div className="space-y-8 py-8">
      <div className="grid gap-8 md:grid-cols-2">
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-3">
                <div className="h-1 w-8 bg-primary rounded-full" />
                <CardTitle>Job {jobId?.slice(0, 8)}</CardTitle>
                <CardDescription className="leading-relaxed">Analysis Status</CardDescription>
              </div>
              <Badge variant={statusBadgeVariant} className="flex items-center gap-2">
                {job?.status === "running" && <Loader2 className="h-3 w-3 animate-spin" />}
                {job?.status === "completed" && <CheckCircle2 className="h-3 w-3" />}
                {job?.status?.toUpperCase()}
              </Badge>
            </div>
          </CardHeader>
          <CardContent className="space-y-6">
            <dl className="grid grid-cols-2 gap-4 text-sm">
              <div>
                <dt className="text-muted-foreground">Molecule</dt>
                <dd className="font-medium">{moleculeName}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Created</dt>
                <dd className="font-medium">{job?.created_at ? new Date(job.created_at).toLocaleTimeString() : "-"}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Status</dt>
                <dd className="font-medium capitalize">{job?.status || "Unknown"}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Updated</dt>
                <dd className="font-medium">{job?.updated_at ? new Date(job.updated_at).toLocaleTimeString() : "-"}</dd>
              </div>
            </dl>
            
            {job?.status === "completed" && (
              <div className="flex flex-col gap-3 pt-4">
                <Link href={`/results?id=${jobId}`}>
                  <Button className="w-full">
                    <CheckCircle2 className="mr-2 h-4 w-4" />
                    View Results
                  </Button>
                </Link>
              </div>
            )}
            
            {(job?.status === "pending" || job?.status === "running") && (
              <div className="flex items-center justify-center py-8">
                <div className="text-center space-y-2">
                  <Loader2 className="h-8 w-8 animate-spin mx-auto text-primary" />
                  <p className="text-sm text-muted-foreground">
                    {job?.status === "pending" ? "Queued for processing..." : "Analyzing molecule..."}
                  </p>
                </div>
              </div>
            )}
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Progress</CardTitle>
          </CardHeader>
          <CardContent className="pt-6">
            <div className="space-y-4">
              {job?.status === "pending" && (
                <div className="flex items-center gap-3">
                  <div className="h-2 w-2 rounded-full bg-yellow-500 animate-pulse" />
                  <span className="text-sm">Waiting in queue...</span>
                </div>
              )}
              {job?.status === "running" && (
                <div className="flex items-center gap-3">
                  <Loader2 className="h-4 w-4 animate-spin text-primary" />
                  <span className="text-sm">Running analysis agents...</span>
                </div>
              )}
              {job?.status === "completed" && (
                <div className="flex items-center gap-3">
                  <CheckCircle2 className="h-4 w-4 text-green-500" />
                  <span className="text-sm">Analysis completed successfully</span>
                </div>
              )}
            </div>
          </CardContent>
        </Card>
      </div>

      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <Terminal className="h-5 w-5 text-muted-foreground" />
              <CardTitle>Job Information</CardTitle>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <div className="rounded-md bg-black/50 border border-border p-4 font-mono text-xs h-[200px] overflow-y-auto">
            <pre className="text-muted-foreground">{JSON.stringify(job, null, 2)}</pre>
          </div>
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
