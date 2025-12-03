import { HISTORY_RUNS } from "../sample-data";

export default function HistoryPage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-3">
        <p className="eyebrow">Workspace</p>
        <h1 className="text-2xl font-semibold">Saved reports</h1>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>Molecule</th>
                <th>Date</th>
                <th>Recommendation</th>
                <th>Status</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody>
              {HISTORY_RUNS.map((run) => (
                <tr key={`${run.molecule}-${run.date}`}>
                  <td>{run.molecule}</td>
                  <td>{run.date}</td>
                  <td>{run.recommendation}</td>
                  <td>{run.status}</td>
                  <td>
                    <div className="table-actions">
                      <button className="btn-secondary" type="button">
                        Open
                      </button>
                      <button className="btn-secondary" type="button">
                        PDF
                      </button>
                    </div>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
