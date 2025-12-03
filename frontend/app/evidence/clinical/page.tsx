import { CLINICAL_TRIALS } from "../../sample-data";

export default function ClinicalEvidencePage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <div>
          <p className="eyebrow">Evidence dashboard</p>
          <h1 className="text-2xl font-semibold">Clinical trials</h1>
          <p className="subtle-text">Filtered view of trial and registry evidence powering the recommendation.</p>
        </div>
        <div className="filters-bar">
          <span className="chip">Phase: All</span>
          <span className="chip">Status: Recruiting + Completed</span>
          <span className="chip">Condition: PF-ILD</span>
        </div>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>NCT ID</th>
                <th>Phase</th>
                <th>Status</th>
                <th>Condition</th>
                <th>Outcome</th>
              </tr>
            </thead>
            <tbody>
              {CLINICAL_TRIALS.map((trial) => (
                <tr key={trial.nctId}>
                  <td>{trial.nctId}</td>
                  <td>{trial.phase}</td>
                  <td>{trial.status}</td>
                  <td>{trial.condition}</td>
                  <td>{trial.outcome}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <section className="grid-two">
        {CLINICAL_TRIALS.map((trial) => (
          <article key={`${trial.nctId}-card`} className="section-card space-y-2">
            <p className="eyebrow">{trial.nctId}</p>
            <h2 className="text-xl font-semibold">{trial.condition}</h2>
            <p className="subtle-text">
              Phase {trial.phase} · {trial.status}
            </p>
            <p>{trial.population} | Sample size {trial.sampleSize}</p>
            <p className="subtle-text">Outcome: {trial.outcome}</p>
            <a className="evidence-link" href={trial.source} target="_blank" rel="noreferrer">
              Source ↗
            </a>
          </article>
        ))}
      </section>
    </div>
  );
}
