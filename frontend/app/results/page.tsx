import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { DEMO_JOB, MARKET_METRICS, SAMPLE_PAYLOAD } from "../sample-data";

export default function ResultsPage() {
  const workers = SAMPLE_PAYLOAD.workers;

  return (
    <div className="section-stack">
      <section className="grid-two">
        <div className="section-card space-y-4">
          <p className="eyebrow">Innovation story</p>
          <h1 className="text-2xl font-semibold">Pirfenidone overview</h1>
          <p className="subtle-text">Recommendation: {SAMPLE_PAYLOAD.recommendation}</p>
          <p>{SAMPLE_PAYLOAD.innovation_story}</p>
          <div className="flex gap-3">
            <Link className="btn-primary" href="/reports">
              Download-ready report
            </Link>
            <Link className="btn-secondary" href="/evidence/clinical">
              Dive into evidence
            </Link>
          </div>
        </div>
        <div className="section-card">
          <h2 className="section-title">Key figures</h2>
          <div className="metrics-grid">
            {MARKET_METRICS.map((metric) => (
              <div key={metric.label} className="metric-card space-y-2">
                <p className="eyebrow">{metric.label}</p>
                <p className="text-2xl font-semibold">{metric.value}</p>
                <p className="subtle-text">{metric.description}</p>
              </div>
            ))}
          </div>
        </div>
      </section>

      <EvidenceTabs workers={workers} reportJobId={DEMO_JOB.id} />
    </div>
  );
}
