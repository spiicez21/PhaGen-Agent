"use client";

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { FileText, ExternalLink, Clock, CheckCircle2, AlertCircle, Loader2, Trash2 } from "lucide-react";
import Link from "next/link";
import { jobStore, StoredJob } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

function HistoryPage() {
  const [jobs, setJobs] = useState<StoredJob[]>([]);

  // Load jobs only on client side to avoid hydration mismatch
  useEffect(() => {
    const storedJobs = jobStore.getAll();
    // Use setTimeout to defer setState and avoid warning
    const timer = setTimeout(() => setJobs(storedJobs), 0);
    return () => clearTimeout(timer);
  }, []);

  const handleClearHistory = () => {
    if (confirm("Are you sure you want to clear all history?")) {
      jobStore.clear();
      setJobs([]);
    }
  };

  const formatDate = (dateString: string) => {
    try {
      return new Date(dateString).toLocaleString();
    } catch {
      return dateString;
    }
  };

  const getStatusBadge = (status: string) => {
    switch (status) {
      case "completed":
        return <Badge variant="default" className="bg-green-600">Completed</Badge>;
      case "running":
        return <Badge variant="default" className="bg-blue-600">Running</Badge>;
      case "pending":
        return <Badge variant="secondary">Pending</Badge>;
      case "failed":
        return <Badge variant="destructive">Failed</Badge>;
      default:
        return <Badge variant="outline">{status}</Badge>;
    }
  };

  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div className="space-y-3">
              <div className="h-1 w-8 bg-primary rounded-full" />
              <CardTitle className="text-2xl">Job History</CardTitle>
              <CardDescription className="leading-relaxed">
                Access your past analysis runs and review results.{jobs.length > 0 && ` (${jobs.length} total)`}
              </CardDescription>
            </div>
            {jobs.length > 0 && (
              <Button variant="outline" size="sm" onClick={handleClearHistory}>
                <Trash2 className="h-4 w-4 mr-2" />
                Clear History
              </Button>
            )}
          </div>
        </CardHeader>
        <CardContent>
          {jobs.length === 0 ? (
            <div className="text-center py-12 text-muted-foreground">
              <AlertCircle className="h-12 w-12 mx-auto mb-4 opacity-50" />
              <p className="text-lg font-medium mb-2">No jobs yet</p>
              <p className="text-sm mb-4">Start by submitting a new molecule analysis</p>
              <Link href="/molecule">
                <Button>Start New Analysis</Button>
              </Link>
            </div>
          ) : (
            <div className="rounded-md border">
              <table className="w-full text-sm text-left">
                <thead className="bg-muted/50 text-muted-foreground font-medium">
                  <tr>
                    <th className="h-12 px-4 align-middle">Molecule / SMILES</th>
                    <th className="h-12 px-4 align-middle">Created</th>
                    <th className="h-12 px-4 align-middle">Recommendation</th>
                    <th className="h-12 px-4 align-middle">Status</th>
                    <th className="h-12 px-4 align-middle text-right">Actions</th>
                  </tr>
                </thead>
                <tbody className="divide-y">
                  {jobs.map((job) => (
                    <tr key={job.job_id} className="hover:bg-muted/50 transition-colors">
                      <td className="p-4 font-medium">
                        <div className="space-y-1">
                          <div>{job.molecule || "Unknown"}</div>
                          <div className="text-xs text-muted-foreground font-mono truncate max-w-[300px]">
                            {job.smiles}
                          </div>
                        </div>
                      </td>
                      <td className="p-4 text-muted-foreground">
                        <div className="flex items-center gap-2">
                          <Clock className="h-3 w-3" />
                          {formatDate(job.created_at)}
                        </div>
                      </td>
                      <td className="p-4">
                        {job.recommendation ? (
                          <Badge variant={job.recommendation === "GO" ? "default" : "destructive"}>
                            {job.recommendation}
                          </Badge>
                        ) : (
                          <span className="text-muted-foreground text-xs">-</span>
                        )}
                      </td>
                      <td className="p-4">
                        <div className="flex items-center gap-2">
                          {job.status === "completed" ? (
                            <CheckCircle2 className="h-4 w-4 text-green-500" />
                          ) : job.status === "running" ? (
                            <Loader2 className="h-4 w-4 text-blue-500 animate-spin" />
                          ) : (
                            <AlertCircle className="h-4 w-4 text-amber-500" />
                          )}
                          {getStatusBadge(job.status)}
                        </div>
                      </td>
                      <td className="p-4 text-right">
                        <div className="flex items-center justify-end gap-2">
                          {job.status === "completed" ? (
                            <>
                              <Button variant="ghost" size="sm" asChild>
                                <Link href={`/results?id=${job.job_id}`}>
                                  <ExternalLink className="h-4 w-4 mr-2" /> View
                                </Link>
                              </Button>
                              <Button variant="outline" size="sm" asChild>
                                <a href={`http://localhost:8000/api/jobs/${job.job_id}/report.pdf`} target="_blank" rel="noopener noreferrer">
                                  <FileText className="h-4 w-4 mr-2" /> PDF
                                </a>
                              </Button>
                            </>
                          ) : (
                            <Button variant="ghost" size="sm" asChild>
                              <Link href={`/job?id=${job.job_id}`}>
                                <Clock className="h-4 w-4 mr-2" /> Status
                              </Link>
                            </Button>
                          )}
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}

export default function ProtectedHistoryPage() {
  return (
    <ProtectedRoute>
      <HistoryPage />
    </ProtectedRoute>
  );
}
