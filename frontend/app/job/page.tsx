import Link from "next/link";
import { JobTimeline } from "../components/JobTimeline";
import { DEMO_JOB, JOB_LOGS, JOB_TIMELINE } from "../sample-data";

export default function JobStatusPage() {
  return (
    <div className="section-stack">
      <section className="grid-two">
        <div className="section-card space-y-3">
          <p className="eyebrow">Job {DEMO_JOB.id}</p>
          <h1 className="text-2xl font-semibold">Agent pipeline status</h1>
          <dl className="metadata-grid">
            <div>
              <dt>Status</dt>
              <dd>{DEMO_JOB.status}</dd>
            </div>
            <div>
              <dt>Depth</dt>
              <dd>{DEMO_JOB.depth}</dd>
            </div>
            <div>
              <dt>Started (UTC)</dt>
              <dd>{new Date(DEMO_JOB.startedAt).toUTCString()}</dd>
            </div>
            <div>
              <dt>ETA</dt>
              <dd>{DEMO_JOB.etaMinutes} min</dd>
            </div>
          </dl>
          <div className="flex gap-3">
            <Link href="/results" className="btn-primary">
              View provisional results
            </Link>
            <Link href="/reports" className="btn-secondary">
              Prepare report
            </Link>
          </div>
        </div>
        <div className="section-card">
          <JobTimeline steps={JOB_TIMELINE} />
        </div>
      </section>

      <section className="section-card space-y-3">
        <div className="flex items-center justify-between">
          <h2 className="text-xl font-semibold">Live logs</h2>
          <span className="chip">Streaming</span>
        </div>
        <div className="table-shell">
          <table>
            <tbody>
              {JOB_LOGS.map((log) => (
                <tr key={log}>
                  <td>{log}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
