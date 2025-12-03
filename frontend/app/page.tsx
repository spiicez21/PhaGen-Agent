import Link from "next/link";
import { HISTORY_RUNS, MARKET_METRICS, SAMPLE_PAYLOAD } from "./sample-data";

export default function LandingPage() {
  const storyHighlights = Object.entries(SAMPLE_PAYLOAD.workers)
    .slice(0, 3)
    .map(([label, worker]) => ({
      label,
      summary: worker.summary
    }));

  const recent = HISTORY_RUNS.slice(0, 3);

  return (
    <div className="section-stack">
      <section className="glass-card">
        <p className="eyebrow">PhaGen Agentic</p>
        <div className="grid-two items-center">
          <div className="space-y-4">
            <h1 className="text-3xl font-semibold">
              Monotone intelligence stack for molecule repurposing decisions
            </h1>
            <p className="subtle-text">
              A calm, compliance-ready cockpit that orchestrates clinical, literature, patent, and market workers into a single innovation story.
            </p>
            <div className="flex flex-wrap gap-3">
              <Link href="/molecule" className="btn-primary">
                Start new analysis
              </Link>
              <Link href="/history" className="btn-secondary">
                Review saved runs
              </Link>
            </div>
          </div>
          <div className="section-card space-y-3">
            <p className="eyebrow">Latest master summary</p>
            <p className="text-lg text-white">{SAMPLE_PAYLOAD.innovation_story}</p>
            <ul className="space-y-2 text-sm text-white/80">
              {storyHighlights.map((item) => (
                <li key={item.label}>
                  <span className="font-semibold capitalize">{item.label} · </span>
                  {item.summary}
                </li>
              ))}
            </ul>
          </div>
        </div>
      </section>

      <section>
        <h2 className="section-title">System signals</h2>
        <div className="metrics-grid">
          {MARKET_METRICS.map((metric) => (
            <div key={metric.label} className="metric-card space-y-2">
              <p className="eyebrow">{metric.label}</p>
              <p className="text-2xl font-semibold">{metric.value}</p>
              <p className="subtle-text">{metric.description}</p>
            </div>
          ))}
        </div>
      </section>

      <section className="grid-two">
        <div className="section-card space-y-4">
          <p className="eyebrow">Navigator</p>
          <h3 className="text-xl font-semibold">Evidence dashboards</h3>
          <p className="subtle-text">
            Review dedicated panels for clinical, literature, patent, and market intelligence. Monotone palette keeps focus on the signal.
          </p>
          <div className="flex flex-wrap gap-2">
            <Link className="chip" href="/evidence/clinical">
              Clinical
            </Link>
            <Link className="chip" href="/evidence/literature">
              Literature
            </Link>
            <Link className="chip" href="/evidence/patent">
              Patent
            </Link>
            <Link className="chip" href="/evidence/market">
              Market
            </Link>
          </div>
        </div>
        <div className="section-card space-y-4">
          <p className="eyebrow">Workspace</p>
          <h3 className="text-xl font-semibold">Reporting & history</h3>
          <p className="subtle-text">
            Export PDF/JSON deliverables in the report viewer, or reopen any past run straight from the history table.
          </p>
          <div className="flex gap-3">
            <Link href="/reports" className="btn-secondary">
              Open report layout
            </Link>
            <Link href="/history" className="btn-secondary">
              View history
            </Link>
          </div>
        </div>
      </section>

      <section className="section-card space-y-4">
        <div className="flex flex-wrap items-center justify-between gap-4">
          <div>
            <p className="eyebrow">Recent molecules</p>
            <h3 className="text-xl font-semibold">Latest saved runs</h3>
          </div>
          <Link href="/history" className="btn-secondary">
            See all
          </Link>
        </div>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>Molecule</th>
                <th>Date</th>
                <th>Recommendation</th>
                <th>Status</th>
              </tr>
            </thead>
            <tbody>
              {recent.map((run) => (
                <tr key={run.molecule}>
                  <td>{run.molecule}</td>
                  <td>{run.date}</td>
                  <td>{run.recommendation}</td>
                  <td>{run.status}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
