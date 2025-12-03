import { MARKET_COMPETITORS, MARKET_METRICS } from "../../sample-data";

export default function MarketEvidencePage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Evidence dashboard</p>
        <h1 className="text-2xl font-semibold">Market viability</h1>
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

      <section className="section-card space-y-3">
        <div className="flex items-center justify-between">
          <h2 className="text-xl font-semibold">Competitor landscape</h2>
          <span className="chip">Updated weekly</span>
        </div>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>Drug</th>
                <th>Company</th>
                <th>Status</th>
                <th>Share</th>
              </tr>
            </thead>
            <tbody>
              {MARKET_COMPETITORS.map((row) => (
                <tr key={row.drug}>
                  <td>{row.drug}</td>
                  <td>{row.company}</td>
                  <td>{row.status}</td>
                  <td>{row.share}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="chart-placeholder">Incidence and prevalence charts render here.</div>
      </section>
    </div>
  );
}
