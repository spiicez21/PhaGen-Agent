import { REPORT_SECTIONS } from "../sample-data";

export default function ReportsPage() {
  return (
    <div className="section-stack">
      <section className="section-card space-y-4">
        <p className="eyebrow">Report workspace</p>
        <h1 className="text-2xl font-semibold">Full deliverable layout</h1>
        <p className="subtle-text">
          Review, edit, and export the structured report. PDF + JSON export hooks connect to backend automation.
        </p>
        <div className="flex gap-3">
          <button className="btn-primary" type="button">
            Download PDF
          </button>
          <button className="btn-secondary" type="button">
            Export JSON
          </button>
        </div>
      </section>

      <section className="section-card space-y-4">
        <p className="eyebrow">Table of contents</p>
        <ol className="space-y-2 list-decimal pl-6">
          {REPORT_SECTIONS.map((section) => (
            <li key={section.title}>{section.title}</li>
          ))}
        </ol>
      </section>

      {REPORT_SECTIONS.map((section) => (
        <section key={section.title} className="section-card space-y-2">
          <h2 className="text-xl font-semibold">{section.title}</h2>
          <p>{section.body}</p>
        </section>
      ))}
    </div>
  );
}
