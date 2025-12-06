"use client";

import type { FC } from "react";
import type { TimelineStep } from "../types";
import { CheckCircle2, Circle, Clock, XCircle } from "lucide-react";
import { cn } from "@/lib/utils";

const StatusIcon = ({ status }: { status: TimelineStep["status"] }) => {
  switch (status) {
    case "completed":
      return <CheckCircle2 className="h-5 w-5 text-green-500" />;
    case "running":
      return <Clock className="h-5 w-5 text-amber-500 animate-pulse" />;
    case "failed":
      return <XCircle className="h-5 w-5 text-destructive" />;
    default:
      return <Circle className="h-5 w-5 text-muted-foreground" />;
  }
};

export const JobTimeline: FC<{ steps: TimelineStep[] }> = ({ steps }) => {
  return (
    <div className="space-y-6">
      <div>
        <h3 className="text-lg font-semibold">Agent Pipeline</h3>
        <p className="text-sm text-muted-foreground">Real-time execution status</p>
      </div>
      <div className="relative space-y-0">
        {steps.map((step, index) => (
          <div key={step.label} className="flex gap-4 pb-8 last:pb-0 relative">
            {index !== steps.length - 1 && (
              <div className="absolute left-2.5 top-6 bottom-0 w-px bg-border" />
            )}
            <div className="relative z-10 bg-background rounded-full">
              <StatusIcon status={step.status} />
            </div>
            <div className="flex-1 pt-0.5">
              <div className="flex items-center justify-between">
                <p className={cn("font-medium", step.status === "pending" && "text-muted-foreground")}>
                  {step.label}
                </p>
                <span className={cn(
                  "text-xs px-2 py-0.5 rounded-full capitalize",
                  step.status === "completed" && "bg-green-500/10 text-green-500",
                  step.status === "running" && "bg-amber-500/10 text-amber-500",
                  step.status === "failed" && "bg-destructive/10 text-destructive",
                  step.status === "pending" && "bg-muted text-muted-foreground"
                )}>
                  {step.status}
                </span>
              </div>
              {step.description && (
                <p className="text-sm text-muted-foreground mt-1">{step.description}</p>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};
