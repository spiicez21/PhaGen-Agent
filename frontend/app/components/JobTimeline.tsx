"use client";

import type { FC } from "react";

export type TimelineStep = {
  label: string;
  status: "pending" | "running" | "completed" | "failed";
};

const statusColor: Record<TimelineStep["status"], string> = {
  pending: "text-slate-400",
  running: "text-amber-400",
  completed: "text-emerald-400",
  failed: "text-rose-500"
};

export const JobTimeline: FC<{ steps: TimelineStep[] }> = ({ steps }) => {
  return (
    <div className="card space-y-4">
      <h3 className="text-lg font-semibold">Job timeline</h3>
      <ol className="space-y-2">
        {steps.map((step) => (
          <li key={step.label} className={`flex items-center gap-3 ${statusColor[step.status]}`}>
            <span className="w-2 h-2 rounded-full bg-current" />
            <span>{step.label}</span>
          </li>
        ))}
      </ol>
    </div>
  );
};
