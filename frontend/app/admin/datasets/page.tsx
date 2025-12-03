import { INDEX_VERSIONS } from "../../sample-data";

export default function DatasetManagerPage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-3">
        <p className="eyebrow">System intelligence</p>
        <h1 className="text-2xl font-semibold">Dataset / index manager</h1>
        <div className="flex gap-3">
          <button className="btn-primary" type="button">
            Rebuild active index
          </button>
          <button className="btn-secondary" type="button">
            Purge stale embeddings
          </button>
        </div>
      </section>

      <section className="section-card space-y-3">
        <h2 className="text-xl font-semibold">Index versions</h2>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>Version</th>
                <th>Date</th>
                <th>Size</th>
                <th>Status</th>
              </tr>
            </thead>
            <tbody>
              {INDEX_VERSIONS.map((version) => (
                <tr key={version.version}>
                  <td>{version.version}</td>
                  <td>{version.date}</td>
                  <td>{version.size}</td>
                  <td>{version.status}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
