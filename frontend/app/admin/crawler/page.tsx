import { CRAWLER_QUEUE, ROBOTS_STATUS } from "../../sample-data";

export default function CrawlerStatusPage() {
  const stats = {
    active: 4,
    completed: 320,
    failed: 2
  };

  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">System intelligence</p>
        <h1 className="text-2xl font-semibold">Crawler status</h1>
        <div className="metrics-grid">
          <div className="metric-card"><p className="eyebrow">Active</p><p className="text-2xl font-semibold">{stats.active}</p></div>
          <div className="metric-card"><p className="eyebrow">Completed</p><p className="text-2xl font-semibold">{stats.completed}</p></div>
          <div className="metric-card"><p className="eyebrow">Failed</p><p className="text-2xl font-semibold">{stats.failed}</p></div>
        </div>
      </section>

      <section className="section-card space-y-3">
        <div className="flex items-center justify-between">
          <h2 className="text-xl font-semibold">Queue</h2>
          <button className="btn-secondary" type="button">
            Add URL
          </button>
        </div>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>URL</th>
                <th>Status</th>
                <th>Retries</th>
              </tr>
            </thead>
            <tbody>
              {CRAWLER_QUEUE.map((item) => (
                <tr key={item.url}>
                  <td className="truncate max-w-xs">{item.url}</td>
                  <td>{item.status}</td>
                  <td>{item.retries}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <section className="section-card space-y-3">
        <h2 className="text-xl font-semibold">Robots.txt summary</h2>
        <div className="table-shell">
          <table>
            <thead>
              <tr>
                <th>Domain</th>
                <th>Access</th>
                <th>Note</th>
              </tr>
            </thead>
            <tbody>
              {ROBOTS_STATUS.map((row) => (
                <tr key={row.domain}>
                  <td>{row.domain}</td>
                  <td>{row.access}</td>
                  <td>{row.note}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>
    </div>
  );
}
