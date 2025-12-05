import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { DEMO_JOB, MARKET_METRICS, SAMPLE_PAYLOAD } from "../sample-data";

export default function ResultsPage() {
  const workers = SAMPLE_PAYLOAD.workers;
  const validation = SAMPLE_PAYLOAD.validation;
  const reportVersion = SAMPLE_PAYLOAD.report_version ?? 1;
  const structure = SAMPLE_PAYLOAD.structure;
  const quality = SAMPLE_PAYLOAD.quality;

  const qualityStatusTone: Record<string, string> = {
    pass: "bg-emerald-100 text-emerald-800",
    needs_attention: "bg-amber-100 text-amber-800",
    investigate: "bg-rose-100 text-rose-800"
  };

  const apiBudgetTone: Record<string, string> = {
    ok: "bg-emerald-100 text-emerald-800",
    warning: "bg-amber-100 text-amber-800",
    exceeded: "bg-rose-100 text-rose-800"
  };

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

      {structure && (
        <section className="section-card space-y-4">
          <div>
            <p className="eyebrow">Molecular structure</p>
            <h2 className="text-xl font-semibold">SMILES preview</h2>
            <p className="subtle-text">{structure.smiles}</p>
            {structure.source_type && structure.source_reference && (
              <p className="text-sm text-neutral-600">
                Source: {structure.source_type.toUpperCase()} · {structure.source_reference}
              </p>
            )}
            {structure.generated_at && (
              <p className="text-xs text-neutral-500">Generated {structure.generated_at}</p>
            )}
            {structure.license && (
              <p className="text-xs text-neutral-500">License: {structure.license}</p>
            )}
            {structure.image_id && (
              <p className="text-xs text-neutral-500">Image ID: {structure.image_id}</p>
            )}
          </div>
          {structure.svg ? (
            <div
              className="rounded-2xl border border-neutral-200 bg-white p-4"
              dangerouslySetInnerHTML={{ __html: structure.svg }}
            />
          ) : (
            <p className="text-sm text-red-700">{structure.error ?? "Structure unavailable."}</p>
          )}
          {structure.path && (
            <p className="subtle-text">SVG asset: {structure.path}</p>
          )}
          {structure.metadata_path && (
            <p className="subtle-text">Provenance: {structure.metadata_path}</p>
          )}
        </section>
      )}

      {quality && (
        <section className="section-card space-y-4">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <div>
              <p className="eyebrow">Quality guardrails</p>
              <h2 className="text-xl font-semibold">Retrieval health</h2>
            </div>
            <span
              className={`px-3 py-1 rounded-full text-sm font-semibold ${qualityStatusTone[quality.status] ?? "bg-neutral-200 text-neutral-800"}`}
            >
              {quality.status === "pass"
                ? "Pass"
                : quality.status === "needs_attention"
                  ? "Needs attention"
                  : "Investigate"}
            </span>
          </div>
          <p className="subtle-text">
            Metrics and alerts reflect how much evidence each worker gathered. Use them to catch thin coverage or anomalous retrieval runs before sharing results.
          </p>

          {Object.keys(quality.alerts).length ? (
            <div className="space-y-3">
              {Object.entries(quality.alerts).map(([worker, alerts]) => (
                <article key={worker} className="rounded-2xl border border-amber-200 bg-amber-50 p-3">
                  <header className="flex items-center gap-2 text-sm font-semibold text-amber-900">
                    <span className="badge">{worker}</span>
                    <span>Alerts</span>
                  </header>
                  <ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-amber-900">
                    {alerts.map((alert, index) => (
                      <li key={`${worker}-alert-${index}`}>{alert}</li>
                    ))}
                  </ul>
                </article>
              ))}
            </div>
          ) : (
            <p className="rounded-2xl border border-emerald-100 bg-emerald-50 p-3 text-sm text-emerald-900">
              All workers met the minimum evidence thresholds. No alerts at this time.
            </p>
          )}

          <div className="grid-two">
            {Object.entries(quality.metrics).map(([worker, metrics]) => (
              <article key={worker} className="space-y-2 rounded-2xl border border-neutral-200 p-4">
                <header className="flex items-center justify-between">
                  <span className="badge">{worker}</span>
                  <span className="text-sm text-neutral-500">
                    {metrics.evidence_count} evidence · {metrics.unique_sources} sources
                  </span>
                </header>
                <dl className="grid grid-cols-2 gap-2 text-sm">
                  <div>
                    <dt className="subtle-text">Coverage</dt>
                    <dd className="font-semibold">{(metrics.coverage_ratio * 100).toFixed(0)}%</dd>
                  </div>
                  <div>
                    <dt className="subtle-text">Precision proxy</dt>
                    <dd className="font-semibold">{(metrics.precision_proxy * 100).toFixed(0)}%</dd>
                  </div>
                  <div>
                    <dt className="subtle-text">Passages pulled</dt>
                    <dd className="font-semibold">{metrics.final_passages}/{metrics.retriever_top_k}</dd>
                  </div>
                  <div>
                    <dt className="subtle-text">High-confidence cites</dt>
                    <dd className="font-semibold">{metrics.high_conf_evidence}</dd>
                  </div>
                </dl>
              </article>
            ))}
          </div>

          {quality.story_checks && (
            <div className="space-y-3 rounded-2xl border border-neutral-200 p-4">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div>
                  <p className="eyebrow">Story guardrail</p>
                  <h3 className="text-lg font-semibold">Hallucination scan</h3>
                </div>
                <span
                  className={`px-3 py-1 rounded-full text-sm font-semibold ${qualityStatusTone[quality.story_checks.status] ?? "bg-neutral-200 text-neutral-800"}`}
                >
                  {quality.story_checks.status === "pass"
                    ? "Pass"
                    : quality.story_checks.status === "needs_attention"
                      ? "Needs attention"
                      : "Investigate"}
                </span>
              </div>
              <p className="subtle-text">
                Every innovation-story sentence is compared against its cited evidence. Low lexical overlap or missing citations are flagged so reviewers can pause before sharing results.
              </p>
              {quality.story_checks.flagged_claims.length ? (
                <ul className="space-y-3">
                  {quality.story_checks.flagged_claims.map((claim) => (
                    <li key={claim.claim_id} className="rounded-2xl border border-rose-200 bg-rose-50 p-3">
                      <p className="text-sm font-semibold text-rose-900">{claim.claim_text}</p>
                      <p className="text-sm text-rose-900">{claim.reason}</p>
                      <p className="text-xs text-rose-800">
                        Support score {(claim.support_score * 100).toFixed(0)}% · Evidence IDs {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                      </p>
                    </li>
                  ))}
                </ul>
              ) : (
                <p className="rounded-2xl border border-emerald-100 bg-emerald-50 p-3 text-sm text-emerald-900">
                  All innovation-story claims have citations with sufficient lexical overlap. No hallucination risk detected.
                </p>
              )}
            </div>
          )}

          {quality.api_budgets && (
            <div className="space-y-3 rounded-2xl border border-neutral-200 p-4">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div>
                  <p className="eyebrow">API budgets</p>
                  <h3 className="text-lg font-semibold">NCBI & OpenFDA quotas</h3>
                </div>
                <span className="text-sm text-neutral-500">
                  Live view of minute/day usage vs documented limits
                </span>
              </div>
              <div className="grid-two">
                {Object.entries(quality.api_budgets).map(([provider, budget]) => (
                  <article key={provider} className="space-y-2 rounded-2xl border border-neutral-200 p-4">
                    <header className="flex items-center justify-between">
                      <span className="badge">{provider}</span>
                      <span
                        className={`px-2 py-0.5 rounded-full text-xs font-semibold ${apiBudgetTone[budget.status] ?? "bg-neutral-200 text-neutral-800"}`}
                      >
                        {budget.status === "ok"
                          ? "Within budget"
                          : budget.status === "warning"
                            ? "Near limit"
                            : "Exceeded"}
                      </span>
                    </header>
                    <p className="text-sm text-neutral-600">{budget.label}</p>
                    <dl className="grid grid-cols-2 gap-2 text-sm">
                      <div>
                        <dt className="subtle-text">Minute window</dt>
                        <dd className="font-semibold">
                          {budget.minute_used}/{budget.per_minute_limit}
                        </dd>
                        <p className="text-xs text-neutral-500">
                          Resets in {budget.reset_in_seconds.minute}s
                        </p>
                      </div>
                      <div>
                        <dt className="subtle-text">Day window</dt>
                        <dd className="font-semibold">
                          {budget.day_used}/{budget.per_day_limit}
                        </dd>
                        <p className="text-xs text-neutral-500">
                          Resets in {Math.round(budget.reset_in_seconds.day / 3600)}h
                        </p>
                      </div>
                    </dl>
                    <p className="text-xs text-neutral-500">
                      Last call {budget.last_call_at ?? "n/a"} · Lifetime {budget.lifetime_calls} requests
                    </p>
                  </article>
                ))}
              </div>
            </div>
          )}
        </section>
      )}

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
