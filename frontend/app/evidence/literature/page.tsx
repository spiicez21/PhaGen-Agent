import { LITERATURE_EVIDENCE } from "../../sample-data";

export default function LiteratureEvidencePage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Evidence dashboard</p>
        <h1 className="text-2xl font-semibold">Literature intelligence</h1>
        <div className="filters-bar">
          <span className="chip">Type: Clinical + Mechanistic</span>
          <span className="chip">Strength ≥ 0.6</span>
        </div>
        <div className="evidence-stack">
          {LITERATURE_EVIDENCE.map((entry) => (
            <article key={entry.title} className="evidence-card">
              <div className="flex items-center justify-between">
                <div>
                  <p className="eyebrow">{entry.doi}</p>
                  <h2 className="text-xl font-semibold">{entry.title}</h2>
                </div>
                <span className="confidence-pill">{Math.round(entry.strength * 100)}% strength</span>
              </div>
              <p className="subtle-text">Mechanism: {entry.mechanism}</p>
              <p>{entry.snippet}</p>
              <a className="evidence-link" href={`https://doi.org/${entry.doi}`} target="_blank" rel="noreferrer">
                View DOI ↗
              </a>
            </article>
          ))}
        </div>
      </section>
    </div>
  );
}
