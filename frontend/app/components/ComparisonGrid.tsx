"use client";

import Link from "next/link";
import type { ComparisonSlot, WorkerResultPayload } from "../types";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ArrowRight, FileText, Activity, CheckCircle2, AlertCircle, HelpCircle } from "lucide-react";

interface ComparisonGridProps {
  slots: ComparisonSlot[];
}

const getConfidenceBadge = (band: WorkerResultPayload["confidence_band"]) => {
  const config = {
    high: { label: "High Confidence", icon: CheckCircle2, variant: "default" as const },
    medium: { label: "Medium Confidence", icon: AlertCircle, variant: "secondary" as const },
    low: { label: "Low Confidence", icon: HelpCircle, variant: "outline" as const }
  };
  
  const { label, icon: Icon, variant } = config[band] || config.low;
  
  return (
    <Badge variant={variant} className="flex items-center gap-1 w-fit text-xs">
      <Icon className="h-3 w-3" /> {label}
    </Badge>
  );
};

const formatWorkerLabel = (key: string): string => key.charAt(0).toUpperCase() + key.slice(1);

export function ComparisonGrid({ slots }: ComparisonGridProps) {
  if (!slots.length) {
    return (
      <div className="flex flex-col items-center justify-center py-12 text-center space-y-4 border rounded-lg bg-muted/10 border-dashed">
        <div className="p-4 rounded-full bg-muted">
          <Activity className="h-8 w-8 text-muted-foreground" />
        </div>
        <div className="space-y-2">
          <h3 className="text-lg font-medium">No molecules selected yet</h3>
          <p className="text-muted-foreground max-w-sm mx-auto">
            Enter at least two job IDs to compare their evidence side by side.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="grid gap-6 md:grid-cols-2 xl:grid-cols-3">
      {slots.map((slot) => (
        <Card key={slot.jobId} className="flex flex-col h-full">
          <CardHeader className="space-y-4">
            <div className="flex items-start justify-between">
              <div className="space-y-1">
                <div className="text-xs text-muted-foreground font-mono">
                  {slot.lastUpdated}
                </div>
                <CardTitle className="text-xl">{slot.molecule}</CardTitle>
              </div>
              <Badge variant={slot.payload.recommendation === "GO" ? "default" : "destructive"}>
                {slot.payload.recommendation}
              </Badge>
            </div>
            
            <div className="flex flex-wrap gap-2">
              <Badge variant="outline" className="font-mono text-xs">
                Job: {slot.jobId}
              </Badge>
              {(slot.reportVersion ?? slot.payload.report_version) && (
                <Badge variant="secondary" className="text-xs">
                  V{slot.reportVersion ?? slot.payload.report_version}
                </Badge>
              )}
              {slot.payload.validation && (
                <Badge 
                  variant={slot.payload.validation.status === "pass" ? "outline" : "destructive"}
                  className={slot.payload.validation.status === "pass" ? "text-green-600 border-green-600" : ""}
                >
                  {slot.payload.validation.claims_linked}/{slot.payload.validation.claims_total} claims
                </Badge>
              )}
            </div>
          </CardHeader>
          
          <CardContent className="space-y-6 flex-1">
            <div className="space-y-2">
              <h4 className="text-sm font-medium text-muted-foreground">Innovation Story</h4>
              <p className="text-sm line-clamp-4 text-muted-foreground">
                {slot.payload.innovation_story}
              </p>
            </div>

            <div className="grid grid-cols-2 gap-4 py-4 border-y">
              <div className="space-y-1">
                <span className="text-xs text-muted-foreground uppercase tracking-wider">Market Score</span>
                <p className="text-2xl font-bold">{slot.payload.market_score}</p>
              </div>
              <div className="space-y-1">
                <span className="text-xs text-muted-foreground uppercase tracking-wider">Evidence Panels</span>
                <p className="text-2xl font-bold">{Object.keys(slot.payload.workers).length}</p>
              </div>
            </div>

            <div className="space-y-4">
              <h4 className="text-sm font-medium">Worker Analysis</h4>
              {Object.entries(slot.payload.workers).map(([workerKey, worker]) => (
                <div key={`${slot.jobId}-${workerKey}`} className="space-y-2 p-3 rounded-lg bg-muted/50">
                  <div className="flex items-center justify-between">
                    <span className="text-sm font-medium">{formatWorkerLabel(workerKey)}</span>
                    {getConfidenceBadge(worker.confidence_band)}
                  </div>
                  <p className="text-xs text-muted-foreground line-clamp-2">
                    {worker.summary}
                  </p>
                  {worker.evidence?.[0] && (
                    <div className="pt-2 mt-2 border-t border-border/50">
                      <p className="text-[10px] text-muted-foreground">
                        <span className="font-semibold text-primary">Top Cite:</span> {worker.evidence[0].text.substring(0, 100)}...
                      </p>
                    </div>
                  )}
                </div>
              ))}
            </div>

            <div className="flex gap-2 pt-4 mt-auto">
              <Button variant="outline" size="sm" className="flex-1" asChild>
                <Link href="/reports">
                  <FileText className="mr-2 h-3 w-3" /> Report
                </Link>
              </Button>
              <Button variant="outline" size="sm" className="flex-1" asChild>
                <Link href={`/job?jobId=${slot.jobId}`}>
                  Timeline <ArrowRight className="ml-2 h-3 w-3" />
                </Link>
              </Button>
            </div>
          </CardContent>
        </Card>
      ))}
    </div>
  );
}
