"use client";

import { useMemo } from "react";
import type { WorkerKind, WorkerResultPayload } from "../types";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ScrollArea } from "@/components/ui/scroll-area";
import { FileText, ExternalLink, AlertCircle, CheckCircle2, HelpCircle } from "lucide-react";

const PANEL_META: Record<WorkerKind, { label: string; description: string; empty: string }> = {
  clinical: {
    label: "Clinical",
    description: "Signals from trials, registries, and observational programs.",
    empty: "No clinical snippets yet. Run the crawler + index builder, then re-run this molecule."
  },
  literature: {
    label: "Literature",
    description: "Mechanism-of-action insights across preclinical and clinical papers.",
    empty: "Literature worker has no evidence for this run."
  },
  patent: {
    label: "Patent & Regulatory",
    description: "Freedom-to-operate, priority dates, and regulatory guardrails.",
    empty: "Patent worker metadata is empty. Add USPTO/FDA corpora to populate it."
  },
  market: {
    label: "Market",
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

const getConfidenceBadge = (band?: WorkerResultPayload["confidence_band"], score?: number) => {
  let label = "Signal needs validation";
  let variant: "default" | "secondary" | "destructive" | "outline" = "outline";
  let icon = <HelpCircle className="w-3 h-3 mr-1" />;

  if (band === "high" || (score !== undefined && score >= 0.8)) {
    label = "High confidence";
    variant = "default"; // Using default (primary color) for high confidence
    icon = <CheckCircle2 className="w-3 h-3 mr-1" />;
  } else if (band === "medium" || (score !== undefined && score >= 0.6)) {
    label = "Moderate confidence";
    variant = "secondary";
    icon = <AlertCircle className="w-3 h-3 mr-1" />;
  }

  return (
    <Badge variant={variant} className="flex items-center w-fit">
      {icon} {label}
    </Badge>
  );
};

export const EvidenceTabs = ({ workers = {}, reportJobId, reportVersion }: EvidenceTabsProps) => {
  const workerKeys = useMemo(() => (Object.keys(PANEL_META) as WorkerKind[]), []);
  const available = useMemo(
    () => workerKeys.filter((key) => Boolean(workers[key])),
    [workerKeys, workers]
  );
  
  const defaultTab = available[0] ?? "clinical";

  return (
    <div className="space-y-6">
      <div className="space-y-1">
        <h2 className="text-2xl font-semibold tracking-tight">Evidence Dashboard</h2>
        <p className="text-muted-foreground">
          Detailed breakdown of evidence collected by each specialized worker.
        </p>
      </div>

      <Tabs defaultValue={defaultTab} className="w-full">
        <TabsList className="grid w-full grid-cols-4">
          {workerKeys.map((key) => (
            <TabsTrigger key={key} value={key} disabled={!workers[key]}>
              {PANEL_META[key].label}
            </TabsTrigger>
          ))}
        </TabsList>
        
        {workerKeys.map((key) => {
          const worker = workers[key];
          const meta = PANEL_META[key];
          
          if (!worker) return null;

          return (
            <TabsContent key={key} value={key} className="space-y-6 mt-6">
              <div className="grid gap-6 md:grid-cols-3">
                {/* Summary Column */}
                <div className="md:col-span-1 space-y-6">
                  <Card>
                    <CardHeader>
                      <CardTitle>Summary</CardTitle>
                      <CardDescription>{meta.description}</CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-4">
                      <div className="text-sm text-muted-foreground">
                        {worker.summary}
                      </div>
                      
                      {getConfidenceBadge(worker.confidence_band, worker.confidence)}

                      {reportJobId && (
                        <Button variant="outline" className="w-full" asChild>
                          <a href={`/api/jobs/${reportJobId}/report.pdf`} target="_blank" rel="noreferrer">
                            <FileText className="mr-2 h-4 w-4" />
                            Download PDF {reportVersion ? `(V${reportVersion})` : ''}
                          </a>
                        </Button>
                      )}
                    </CardContent>
                  </Card>

                  {Object.keys(worker.metadata ?? {}).length > 0 && (
                    <Card>
                      <CardHeader>
                        <CardTitle className="text-sm">Metadata</CardTitle>
                      </CardHeader>
                      <CardContent>
                        <dl className="space-y-2 text-sm">
                          {Object.entries(worker.metadata).map(([k, v]) => (
                            <div key={k} className="flex justify-between border-b pb-1 last:border-0 last:pb-0">
                              <dt className="font-medium text-muted-foreground capitalize">{k.replaceAll("_", " ")}</dt>
                              <dd className="text-right truncate max-w-[150px]" title={v}>{formatMetadataValue(v)}</dd>
                            </div>
                          ))}
                        </dl>
                      </CardContent>
                    </Card>
                  )}
                </div>

                {/* Evidence List Column */}
                <div className="md:col-span-2">
                  <Card className="h-full flex flex-col">
                    <CardHeader>
                      <CardTitle>Evidence Items</CardTitle>
                      <CardDescription>
                        {worker.evidence?.length ?? 0} items retrieved
                      </CardDescription>
                    </CardHeader>
                    <CardContent className="flex-1 p-0">
                      {worker.evidence?.length ? (
                        <div className="divide-y">
                          {worker.evidence.map((item, idx) => (
                            <div key={idx} className="p-6 hover:bg-muted/50 transition-colors">
                              <div className="flex items-center justify-between mb-2">
                                <Badge variant="secondary" className="capitalize">
                                  {item.type}
                                </Badge>
                                <span className="text-xs text-muted-foreground font-mono">
                                  {Math.round(item.confidence * 100)}% conf
                                </span>
                              </div>
                              <p className="text-sm leading-relaxed mb-3">
                                {item.text}
                              </p>
                              <div className="flex items-center justify-between">
                                {item.evidence_id && (
                                  <span className="text-xs text-muted-foreground font-mono">
                                    ID: {item.evidence_id}
                                  </span>
                                )}
                                {item.url && (
                                  <a 
                                    href={item.url} 
                                    target="_blank" 
                                    rel="noreferrer" 
                                    className="text-xs flex items-center text-primary hover:underline"
                                  >
                                    Source <ExternalLink className="ml-1 h-3 w-3" />
                                  </a>
                                )}
                              </div>
                            </div>
                          ))}
                        </div>
                      ) : (
                        <div className="p-8 text-center text-muted-foreground">
                          <p>{meta.empty}</p>
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
