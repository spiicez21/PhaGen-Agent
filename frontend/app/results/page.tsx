import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { DEMO_JOB, MARKET_METRICS, SAMPLE_PAYLOAD } from "../sample-data";

export default function ResultsPage() {
  const workers = SAMPLE_PAYLOAD.workers;
  const validation = SAMPLE_PAYLOAD.validation;
  const reportVersion = SAMPLE_PAYLOAD.report_version ?? 1;

  return (
    <div className="section-stack">
      <section className="grid-two">
        <div className="section-card space-y-4">
          <p className="eyebrow">Innovation story</p>
          <h1 className="text-2xl font-semibold">Pirfenidone overview</h1>
          <p className="subtle-text">
            Recommendation: {SAMPLE_PAYLOAD.recommendation} · Report V{reportVersion}
          </p>
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

      <EvidenceTabs workers={workers} reportJobId={DEMO_JOB.id} reportVersion={reportVersion} />

      {validation && (
        <section className="section-card space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <p className="eyebrow">Validation pass</p>
              <h2 className="text-xl font-semibold">Claim traceability</h2>
            </div>
            <span className={`validation-chip validation-chip--${validation.status}`}>
              {validation.claims_linked}/{validation.claims_total} claims linked
            </span>
          </div>
          <p className="subtle-text">
            Every sentence in the innovation story carries an evidence ID so reviewers can trace claims back to primary sources.
          </p>
          <div className="validation-claims">
            {validation.claim_links.map((claim) => (
              <article key={claim.claim_id} className="validation-claim space-y-2">
                <header className="flex flex-wrap items-center gap-2">
                  <span className="badge">{claim.worker}</span>
                  <span className="confidence-pill">
                    {claim.status === "linked" ? "Linked" : "Needs review"}
                  </span>
                </header>
                <p>{claim.claim_text}</p>
                <p className="subtle-text">
                  Evidence IDs · {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                </p>
              </article>
            ))}
          </div>
        </section>
      )}
    </div>
  );
}
