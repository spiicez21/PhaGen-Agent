"use client";

import { useMemo } from "react";
import type { WorkerKind, WorkerResultPayload } from "../types";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { FileText, ExternalLink, CheckCircle2, AlertCircle, HelpCircle, Microscope, BookOpen, Scale, BarChart3, Activity } from "lucide-react";

const PANEL_META: Record<WorkerKind, { label: string; icon: any; description: string; empty: string }> = {
  clinical: {
    label: "Clinical",
    icon: Microscope,
    description: "Signals from trials, registries, and observational programs.",
    empty: "No clinical snippets yet. Run the crawler + index builder, then re-run this molecule."
  },
  literature: {
    label: "Literature",
    icon: BookOpen,
    description: "Mechanism-of-action insights across preclinical and clinical papers.",
    empty: "Literature worker has no evidence for this run."
  },
  patent: {
    label: "Patent & IP",
    icon: Scale,
    description: "Freedom-to-operate, priority dates, and regulatory guardrails.",
    empty: "Patent worker metadata is empty. Add USPTO/FDA corpora to populate it."
  },
  market: {
    label: "Market",
    icon: BarChart3,
    description: "Commercial sizing, competitor snapshots, and unmet-need commentary.",
    empty: "Market worker requires incidence / prevalence data to render cards."
  }
};

interface EvidenceTabsProps {
  workers?: Record<string, WorkerResultPayload>;
  reportJobId?: string;
  reportVersion?: number;
}

const formatMetadataValue = (value: string): string => {
  const trimmed = value?.trim?.() ?? "";
  if (!trimmed) return "--";
  if ((trimmed.startsWith("[") && trimmed.endsWith("]")) || (trimmed.startsWith("{") && trimmed.endsWith("}"))) {
    try {
      const parsed = JSON.parse(trimmed);
      if (Array.isArray(parsed)) return `${parsed.length} entries`;
      if (typeof parsed === "object" && parsed !== null) return `${Object.keys(parsed).length} fields`;
    } catch {
      return value;
    }
  }
  return value;
};

const getConfidenceDisplay = (band?: WorkerResultPayload["confidence_band"], score?: number) => {
  if (band === "high" || (score !== undefined && score >= 0.8)) {
    return { label: "High", color: "text-emerald-500", bg: "bg-emerald-500/10", icon: CheckCircle2 };
  }
  if (band === "medium" || (score !== undefined && score >= 0.6)) {
    return { label: "Moderate", color: "text-amber-500", bg: "bg-amber-500/10", icon: AlertCircle };
  }
  return { label: "Low", color: "text-muted-foreground", bg: "bg-muted/50", icon: HelpCircle };
};

export const EvidenceTabs = ({ workers = {}, reportJobId, reportVersion }: EvidenceTabsProps) => {
  const workerKeys = useMemo(() => Object.keys(PANEL_META) as WorkerKind[], []);
  const available = useMemo(() => workerKeys.filter((key) => Boolean(workers[key])), [workerKeys, workers]);
  const defaultTab = available[0] ?? "clinical";

  return (
    <div className="space-y-4">
      <div className="flex items-center gap-2">
        <Activity className="h-4 w-4 text-primary" />
        <h2 className="text-lg font-semibold tracking-tight">Evidence Dashboard</h2>
      </div>

      <Tabs defaultValue={defaultTab} className="w-full">
        <TabsList className="grid w-full grid-cols-4 bg-muted/30 p-1 rounded-xl h-auto">
          {workerKeys.map((key) => {
            const meta = PANEL_META[key];
            const Icon = meta.icon;
            const hasData = !!workers[key];
            return (
              <TabsTrigger
                key={key}
                value={key}
                disabled={!hasData}
                className="flex items-center gap-1.5 text-xs py-2 rounded-lg data-[state=active]:bg-background data-[state=active]:shadow-sm"
              >
                <Icon className="h-3.5 w-3.5" />
                {meta.label}
                {hasData && (
                  <span className="ml-1 text-[10px] bg-primary/10 text-primary px-1.5 py-0.5 rounded-full font-medium">
                    {workers[key]?.evidence?.length || 0}
                  </span>
                )}
              </TabsTrigger>
            );
          })}
        </TabsList>

        {workerKeys.map((key) => {
          const worker = workers[key];
          const meta = PANEL_META[key];
          if (!worker) return null;

          const conf = getConfidenceDisplay(worker.confidence_band, worker.confidence);
          const ConfIcon = conf.icon;

          return (
            <TabsContent key={key} value={key} className="mt-4 space-y-4">
              <div className="grid gap-4 lg:grid-cols-3">
                {/* Summary sidebar */}
                <div className="space-y-3">
                  <Card className="border-border/50 bg-card/80">
                    <CardContent className="pt-5 pb-4 space-y-3">
                      <p className="text-[10px] uppercase tracking-wider text-muted-foreground font-medium">Summary</p>
                      <p className="text-xs text-muted-foreground leading-relaxed">{worker.summary}</p>

                      <div className={`flex items-center gap-1.5 text-xs font-medium px-2.5 py-1.5 rounded-lg w-fit ${conf.bg} ${conf.color}`}>
                        <ConfIcon className="h-3 w-3" />
                        {conf.label} confidence
                      </div>

                      {reportJobId && (
                        <Button variant="outline" size="sm" className="w-full text-xs" asChild>
                          <a href={`/api/jobs/${reportJobId}/report.pdf`} target="_blank" rel="noreferrer">
                            <FileText className="mr-1.5 h-3 w-3" />
                            PDF {reportVersion ? `v${reportVersion}` : ""}
                          </a>
                        </Button>
                      )}
                    </CardContent>
                  </Card>

                  {Object.keys(worker.metadata ?? {}).length > 0 && (
                    <Card className="border-border/50 bg-card/60">
                      <CardContent className="pt-4 pb-3">
                        <p className="text-[10px] uppercase tracking-wider text-muted-foreground font-medium mb-2">Metadata</p>
                        <dl className="space-y-1.5 text-xs">
                          {Object.entries(worker.metadata).map(([k, v]) => (
                            <div key={k} className="flex justify-between gap-2 py-1 border-b border-border/20 last:border-0">
                              <dt className="text-muted-foreground capitalize truncate">{k.replaceAll("_", " ")}</dt>
                              <dd className="text-right truncate max-w-[120px] font-mono text-[11px]" title={v}>{formatMetadataValue(v)}</dd>
                            </div>
                          ))}
                        </dl>
                      </CardContent>
                    </Card>
                  )}
                </div>

                {/* Evidence items */}
                <div className="lg:col-span-2">
                  <Card className="border-border/50 bg-card/60 h-full">
                    <CardContent className="pt-5 pb-0">
                      <div className="flex items-center justify-between mb-4">
                        <p className="text-[10px] uppercase tracking-wider text-muted-foreground font-medium">
                          Evidence Items
                        </p>
                        <span className="text-[11px] text-muted-foreground/60">
                          {worker.evidence?.length ?? 0} retrieved
                        </span>
                      </div>

                      {worker.evidence?.length ? (
                        <div className="divide-y divide-border/30 -mx-6">
                          {worker.evidence.map((item, idx) => (
                            <div key={idx} className="px-6 py-4 hover:bg-muted/20 transition-colors">
                              <div className="flex items-center justify-between mb-2">
                                <span className="text-[10px] font-medium uppercase tracking-wider text-primary/70 bg-primary/5 px-2 py-0.5 rounded capitalize">
                                  {item.type}
                                </span>
                                <span className={`text-[10px] font-bold px-2 py-0.5 rounded-full ${
                                  item.confidence >= 0.8 ? "bg-emerald-500/10 text-emerald-500"
                                  : item.confidence >= 0.6 ? "bg-amber-500/10 text-amber-500"
                                  : "bg-muted text-muted-foreground"
                                }`}>
                                  {Math.round(item.confidence * 100)}%
                                </span>
                              </div>
                              <p className="text-xs leading-relaxed text-foreground/80 mb-2">{item.text}</p>
                              <div className="flex items-center justify-between">
                                {item.evidence_id && (
                                  <span className="text-[10px] text-muted-foreground/50 font-mono">ID: {item.evidence_id}</span>
                                )}
                                {item.url && (
                                  <a href={item.url} target="_blank" rel="noreferrer" className="text-[10px] flex items-center gap-1 text-primary hover:underline">
                                    Source <ExternalLink className="h-2.5 w-2.5" />
                                  </a>
                                )}
                              </div>
                            </div>
                          ))}
                        </div>
                      ) : (
                        <div className="py-12 text-center text-xs text-muted-foreground/60">
                          {meta.empty}
                        </div>
                      )}
                    </CardContent>
                  </Card>
                </div>
              </div>
            </TabsContent>
          );
        })}
      </Tabs>
    </div>
  );
};
