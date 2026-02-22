"use client";

import { useState, useEffect } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import {
  ArrowRight, Brain, FileText, Layers, Search,
  Zap, CheckCircle2, Loader2, AlertTriangle, Clock, Plus
} from "lucide-react";
import { jobStore, StoredJob } from "@/lib/store";

const FEATURES = [
  { icon: Brain,    title: "ReAct Agents",     desc: "Reason + Act loop for multi-step analysis" },
  { icon: Layers,   title: "4 Domain Workers", desc: "Clinical, literature, patent & market" },
  { icon: Zap,      title: "Live Tool-Use",    desc: "PubMed, ClinicalTrials.gov, Semantic Scholar" },
  { icon: FileText, title: "PDF Reports",      desc: "Citation-backed with claim traceability" },
];

function StatusDot({ status }: { status: string }) {
  const s = status.toLowerCase();
  if (s === "completed") return <span className="inline-block h-1.5 w-1.5 rounded-full bg-green-500" />;
  if (s === "running")   return <span className="inline-block h-1.5 w-1.5 rounded-full bg-blue-500 animate-pulse" />;
  if (s === "failed")    return <span className="inline-block h-1.5 w-1.5 rounded-full bg-red-500" />;
  return <span className="inline-block h-1.5 w-1.5 rounded-full bg-muted-foreground/40" />;
}

function RecBadge({ rec }: { rec?: string }) {
  if (!rec) return null;
  const cls =
    rec === "GO" || rec === "Go" ? "status-go" :
    rec === "No-Go" || rec === "NO-GO" ? "status-no-go" :
    "status-investigate";
  return <span className={`badge-base ${cls}`}>{rec}</span>;
}

function StatusIcon({ status }: { status: string }) {
  const s = status.toLowerCase();
  if (s === "completed") return <CheckCircle2 className="h-3 w-3 text-green-500" />;
  if (s === "running")   return <Loader2 className="h-3 w-3 text-blue-500 animate-spin" />;
  if (s === "failed")    return <AlertTriangle className="h-3 w-3 text-red-500" />;
  return <Clock className="h-3 w-3 text-muted-foreground" />;
}

export default function DashboardPage() {
  const [jobs, setJobs] = useState<StoredJob[]>([]);

  useEffect(() => {
    const t = setTimeout(() => setJobs(jobStore.getAll()), 0);
    return () => clearTimeout(t);
  }, []);

  const completed = jobs.filter(j => j.status === "completed").length;
  const running   = jobs.filter(j => j.status === "running").length;
  const goCount   = jobs.filter(j => j.recommendation === "GO" || j.recommendation === "Go").length;

  return (
    <div className="space-y-5">

      {/* Page title row */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-base font-semibold text-foreground">R&amp;D Dashboard</h1>
          <p className="text-[11px] text-muted-foreground mt-0.5">Molecule repurposing intelligence platform</p>
        </div>
        <Link href="/molecule">
          <Button size="sm" className="h-7 px-3 text-xs bg-primary text-primary-foreground gap-1.5">
            <Plus className="h-3 w-3" />
            New Analysis
          </Button>
        </Link>
      </div>

      {/* Stats row */}
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
        {[
          { label: "Total Runs",   value: jobs.length,  mono: true },
          { label: "Completed",    value: completed,     mono: true },
          { label: "Active",       value: running,       mono: true },
          { label: "GO Decisions", value: goCount,       mono: true },
        ].map(({ label, value }) => (
          <div key={label} className="card-flush px-4 py-3">
            <p className="label-xs mb-1.5">{label}</p>
            <p className="font-mono text-xl font-bold text-foreground">{value}</p>
          </div>
        ))}
      </div>

      {/* Main content: recent jobs + capabilities */}
      <div className="grid lg:grid-cols-3 gap-4">

        {/* Recent analyses table */}
        <div className="lg:col-span-2 card-flush overflow-hidden">
          <div className="flex items-center justify-between px-4 py-2.5 border-b border-border">
            <h2 className="text-xs font-semibold text-foreground">Recent Analyses</h2>
            <Link href="/history" className="text-[11px] text-primary hover:underline flex items-center gap-1">
              View all <ArrowRight className="h-2.5 w-2.5" />
            </Link>
          </div>

          {jobs.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-14 text-center gap-3">
              <Search className="h-8 w-8 text-muted-foreground/30" />
              <div>
                <p className="text-xs font-medium text-muted-foreground">No analyses yet</p>
                <p className="text-[11px] text-muted-foreground/60 mt-0.5">Submit a molecule to begin scouting</p>
              </div>
              <Link href="/molecule">
                <Button size="sm" className="h-7 text-xs bg-primary text-primary-foreground">
                  <Search className="mr-1.5 h-3 w-3" /> Start Analysis
                </Button>
              </Link>
            </div>
          ) : (
            <table className="data-table">
              <thead>
                <tr>
                  <th>Molecule</th>
                  <th>SMILES</th>
                  <th>Status</th>
                  <th>Decision</th>
                  <th>Date</th>
                  <th></th>
                </tr>
              </thead>
              <tbody>
                {jobs.slice(0, 8).map((job) => {
                  const href = job.status === "completed"
                    ? `/results?id=${job.job_id}`
                    : `/job?id=${job.job_id}`;
                  return (
                    <tr key={job.job_id} onClick={() => window.location.href = href}>
                      <td>
                        <div className="flex items-center gap-2">
                          <StatusDot status={job.status} />
                          <span className="font-medium text-foreground text-xs">{job.molecule || "Unknown"}</span>
                        </div>
                      </td>
                      <td>
                        <span className="font-mono text-[11px] text-muted-foreground max-w-[120px] block truncate">
                          {job.smiles}
                        </span>
                      </td>
                      <td>
                        <div className="flex items-center gap-1.5">
                          <StatusIcon status={job.status} />
                          <span className="text-xs capitalize text-muted-foreground">{job.status}</span>
                        </div>
                      </td>
                      <td><RecBadge rec={job.recommendation} /></td>
                      <td>
                        <span className="font-mono text-[10px] text-muted-foreground">
                          {new Date(job.created_at).toLocaleDateString(undefined, {
                            month: "short", day: "numeric"
                          })}
                        </span>
                      </td>
                      <td>
                        <ArrowRight className="h-3 w-3 text-muted-foreground/40" />
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          )}
        </div>

        {/* Capabilities panel */}
        <div className="card-flush overflow-hidden">
          <div className="px-4 py-2.5 border-b border-border">
            <h2 className="text-xs font-semibold text-foreground">Agent Capabilities</h2>
          </div>
          <div className="divide-y divide-border">
            {FEATURES.map(({ icon: Icon, title, desc }) => (
              <div key={title} className="flex items-start gap-3 px-4 py-3 hover:bg-muted/30 transition-colors">
                <div className="h-6 w-6 rounded bg-primary/10 flex items-center justify-center flex-shrink-0 mt-0.5">
                  <Icon className="h-3.5 w-3.5 text-primary" />
                </div>
                <div>
                  <p className="text-xs font-semibold text-foreground">{title}</p>
                  <p className="text-[11px] text-muted-foreground mt-0.5 leading-relaxed">{desc}</p>
                </div>
              </div>
            ))}
          </div>

          {/* Quick actions */}
          <div className="px-4 py-3 border-t border-border space-y-2">
            <p className="label-xs mb-2">Quick Actions</p>
            <Link href="/molecule" className="flex items-center justify-between py-1.5 px-2 rounded hover:bg-muted/40 transition-colors group">
              <div className="flex items-center gap-2">
                <Search className="h-3 w-3 text-primary" />
                <span className="text-xs font-medium">Run Molecule Analysis</span>
              </div>
              <ArrowRight className="h-3 w-3 text-muted-foreground/40 group-hover:text-primary transition-colors" />
            </Link>
            <Link href="/history" className="flex items-center justify-between py-1.5 px-2 rounded hover:bg-muted/40 transition-colors group">
              <div className="flex items-center gap-2">
                <Clock className="h-3 w-3 text-muted-foreground" />
                <span className="text-xs font-medium">Browse History</span>
              </div>
              <ArrowRight className="h-3 w-3 text-muted-foreground/40 group-hover:text-primary transition-colors" />
            </Link>
          </div>
        </div>

      </div>
    </div>
  );
}
