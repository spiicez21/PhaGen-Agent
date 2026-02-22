"use client";

import { useState, useEffect } from "react";
import { Button } from "@/components/ui/button";
import {
  FileText, Clock, CheckCircle2, Loader2,
  Trash2, AlertTriangle, ArrowRight, Search, Plus
} from "lucide-react";
import Link from "next/link";
import { jobStore, StoredJob } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

function HistoryPage() {
  const [jobs, setJobs] = useState<StoredJob[]>([]);

  useEffect(() => {
    const t = setTimeout(() => setJobs(jobStore.getAll()), 0);
    return () => clearTimeout(t);
  }, []);

  const handleClear = () => {
    if (confirm("Clear all analysis history?")) {
      jobStore.clear();
      setJobs([]);
    }
  };

  const formatDate = (d: string) => {
    try {
      return new Date(d).toLocaleDateString(undefined, {
        year: "numeric", month: "short", day: "numeric",
        hour: "2-digit", minute: "2-digit",
      });
    } catch { return d; }
  };

  const getStatusConfig = (status: string) => {
    switch (status) {
      case "completed": return { color: "text-green-500",   icon: CheckCircle2,    label: "Completed" };
      case "running":   return { color: "text-blue-500",    icon: Loader2,         label: "Running"   };
      case "failed":    return { color: "text-destructive", icon: AlertTriangle,   label: "Failed"    };
      default:          return { color: "text-muted-foreground", icon: Clock,      label: "Pending"   };
    }
  };

  const getRecClass = (rec?: string) => {
    if (!rec) return "";
    if (rec === "GO" || rec === "Go") return "status-go";
    if (rec === "No-Go" || rec === "NO-GO") return "status-no-go";
    return "status-investigate";
  };

  const completed = jobs.filter(j => j.status === "completed").length;

  return (
    <div className="space-y-4">

      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-base font-semibold text-foreground">Analysis History</h1>
          <p className="text-[11px] text-muted-foreground mt-0.5">
            <span className="font-mono font-bold text-foreground">{jobs.length}</span> total runs —{" "}
            <span className="font-mono font-bold text-green-500">{completed}</span> completed
          </p>
        </div>
        <div className="flex items-center gap-2">
          {jobs.length > 0 && (
            <button
              onClick={handleClear}
              className="flex items-center gap-1.5 text-[11px] text-muted-foreground hover:text-destructive transition-colors"
            >
              <Trash2 className="h-3 w-3" />
              Clear
            </button>
          )}
          <Link href="/molecule">
            <Button size="sm" className="h-7 px-3 text-xs bg-primary text-primary-foreground gap-1.5">
              <Plus className="h-3 w-3" />
              New
            </Button>
          </Link>
        </div>
      </div>

      {/* Content */}
      {jobs.length === 0 ? (
        <div className="card-flush flex flex-col items-center justify-center py-20 text-center gap-4">
          <Search className="h-10 w-10 text-muted-foreground/20" />
          <div>
            <p className="text-sm font-medium text-muted-foreground">No analyses yet</p>
            <p className="text-[11px] text-muted-foreground/60 mt-1">
              Submit your first molecule to begin building history
            </p>
          </div>
          <Link href="/molecule">
            <Button size="sm" className="h-7 text-xs bg-primary text-primary-foreground gap-1.5">
              <Search className="h-3 w-3" /> Start Analysis
            </Button>
          </Link>
        </div>
      ) : (
        <div className="card-flush overflow-hidden">
          <table className="data-table">
            <thead>
              <tr>
                <th style={{ width: 12 }}></th>
                <th>Molecule</th>
                <th>SMILES</th>
                <th>Status</th>
                <th>Decision</th>
                <th>Job ID</th>
                <th>Date</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody>
              {jobs.map((job) => {
                const sc   = getStatusConfig(job.status);
                const Icon = sc.icon;
                const href = job.status === "completed"
                  ? `/results?id=${job.job_id}`
                  : `/job?id=${job.job_id}`;

                return (
                  <tr
                    key={job.job_id}
                    onClick={() => window.location.href = href}
                    className="group"
                  >
                    {/* Status dot */}
                    <td>
                      <span className={`inline-block h-1.5 w-1.5 rounded-full ${
                        job.status === "completed" ? "bg-green-500" :
                        job.status === "running"   ? "bg-blue-500 animate-pulse" :
                        job.status === "failed"    ? "bg-destructive" :
                        "bg-muted-foreground/30"
                      }`} />
                    </td>

                    {/* Molecule */}
                    <td>
                      <span className="font-medium text-foreground text-xs">
                        {job.molecule || "Unknown"}
                      </span>
                    </td>

                    {/* SMILES */}
                    <td>
                      <span className="font-mono text-[10px] text-muted-foreground block max-w-[180px] truncate">
                        {job.smiles}
                      </span>
                    </td>

                    {/* Status */}
                    <td>
                      <div className="flex items-center gap-1.5">
                        <Icon className={`h-3 w-3 ${sc.color} ${job.status === "running" ? "animate-spin" : ""}`} />
                        <span className={`text-xs font-medium ${sc.color}`}>{sc.label}</span>
                      </div>
                    </td>

                    {/* Decision */}
                    <td>
                      {job.recommendation ? (
                        <span className={`badge-base ${getRecClass(job.recommendation)}`}>
                          {job.recommendation}
                        </span>
                      ) : (
                        <span className="text-[10px] text-muted-foreground/40">—</span>
                      )}
                    </td>

                    {/* Job ID */}
                    <td>
                      <span className="font-mono text-[10px] text-muted-foreground">
                        {job.job_id.slice(0, 8)}
                      </span>
                    </td>

                    {/* Date */}
                    <td>
                      <span className="text-[10px] text-muted-foreground whitespace-nowrap">
                        {formatDate(job.created_at)}
                      </span>
                    </td>

                    {/* Actions */}
                    <td>
                      <div className="flex items-center gap-1" onClick={(e) => e.stopPropagation()}>
                        {job.status === "completed" && (
                          <a
                            href={`/api/jobs/${job.job_id}/report.pdf`}
                            target="_blank"
                            rel="noopener noreferrer"
                            className="p-1 rounded hover:bg-muted/60 transition-colors text-muted-foreground hover:text-primary"
                            title="Download PDF"
                          >
                            <FileText className="h-3 w-3" />
                          </a>
                        )}
                        <span className="p-1 text-muted-foreground/30 group-hover:text-primary transition-colors">
                          <ArrowRight className="h-3 w-3" />
                        </span>
                      </div>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>

          {/* Table footer */}
          <div className="px-4 py-2 border-t border-border bg-muted/20 flex items-center justify-between">
            <p className="text-[10px] text-muted-foreground font-mono">
              {jobs.length} record{jobs.length !== 1 ? "s" : ""}
            </p>
            <p className="text-[10px] text-muted-foreground">
              Click any row to view details
            </p>
          </div>
        </div>
      )}
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
