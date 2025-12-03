import { PATENT_SUMMARY, REGULATORY_NOTES } from "../../sample-data";

export default function PatentEvidencePage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Evidence dashboard</p>
        <h1 className="text-2xl font-semibold">Patent & regulatory signals</h1>
        <div className="grid-two">
          <div className="section-card space-y-2">
            <p className="eyebrow">Patent summary</p>
            <h2 className="text-xl font-semibold">{PATENT_SUMMARY.assignee}</h2>
            <p className="subtle-text">Priority date: {PATENT_SUMMARY.priorityDate}</p>
            <p className="subtle-text">Risk level: {PATENT_SUMMARY.risk}</p>
            <ul className="space-y-2">
              {PATENT_SUMMARY.blockingClaims.map((claim) => (
                <li key={claim}>{claim}</li>
              ))}
            </ul>
          </div>
          <div className="section-card space-y-2">
            <p className="eyebrow">Regulatory notes</p>
            <ul className="space-y-3">
              {REGULATORY_NOTES.map((note) => (
                <li key={note.label}>
                  <p className="text-sm uppercase tracking-[0.2em]">{note.label}</p>
                  <p>{note.detail}</p>
                </li>
              ))}
            </ul>
          </div>
        </div>
      </section>
    </div>
  );
}
