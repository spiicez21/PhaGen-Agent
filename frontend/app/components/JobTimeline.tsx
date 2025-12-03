"use client";

import type { FC } from "react";
import type { TimelineStep } from "../types";

const statusColor: Record<TimelineStep["status"], string> = {
  pending: "border-slate-600/50 text-slate-300",
  running: "border-amber-400/70 text-amber-300",
  completed: "border-emerald-400/70 text-emerald-300",
  failed: "border-rose-500/70 text-rose-400"
};

export const JobTimeline: FC<{ steps: TimelineStep[] }> = ({ steps }) => {
  return (
    <div className="space-y-4">
      <div>
        <p className="eyebrow">Agent pipeline</p>
        <h3 className="text-lg font-semibold">Job timeline</h3>
      </div>
      <ol className="timeline">
        {steps.map((step) => (
          <li key={step.label} className={`timeline__item ${statusColor[step.status]}`}>
            <span className="timeline__bullet" />
            <div>
              <p className="font-medium">{step.label}</p>
              {step.description && <p className="timeline__description">{step.description}</p>}
            </div>
          </li>
        ))}
      </ol>
    </div>
  );
};
