"use client";

import { useCallback, useState } from "react";
import { JobTimeline } from "../components/JobTimeline";
import { JOB_TIMELINE } from "../sample-data";

const API_BASE = process.env.NEXT_PUBLIC_API_URL ?? "http://localhost:8000";

type JobResponse = {
  job_id: string;
  status: string;
};

export default function MoleculePage() {
  const [molecule, setMolecule] = useState("Pirfenidone");
  const [synonyms, setSynonyms] = useState("Esbriet\nPFD");
  const [smiles, setSmiles] = useState("");
  const [inchikey, setInchikey] = useState("KUFJONJOBWVSNK-UHFFFAOYSA-N");
  const [depth, setDepth] = useState<"quick" | "full">("full");
  const [note, setNote] = useState("Focus on PF-ILD cohort.");
  const [jobId, setJobId] = useState<string | null>(null);
  const [statusMessage, setStatusMessage] = useState<string>("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const submit = useCallback(
    async (event: React.FormEvent<HTMLFormElement>) => {
      event.preventDefault();
      setLoading(true);
      setError(null);
      setStatusMessage("Queueing job...");
      const synonymList = synonyms
        .split(/[\n,]/)
        .map((value) => value.trim())
        .filter(Boolean);

      try {
        const response = await fetch(`${API_BASE}/api/jobs`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ molecule, synonyms: synonymList, smiles, inchikey, depth, note })
        });

        if (!response.ok) {
          throw new Error("Unable to enqueue job");
        }

        const payload: JobResponse = await response.json();
        setJobId(payload.job_id);
        setStatusMessage(`Job ${payload.job_id} queued with status ${payload.status}`);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Unknown error");
      } finally {
        setLoading(false);
      }
    },
    [depth, inchikey, molecule, note, synonyms, smiles]
  );

  return (
    <div className="section-stack">
      <div className="grid-two">
        <form className="section-card space-y-4" onSubmit={submit}>
          <div>
            <p className="eyebrow">Molecule search</p>
            <h1 className="text-2xl font-semibold">Submit a new analysis</h1>
            <p className="subtle-text">
              Provide identifiers or upload structured context -- the platform will normalise inputs before orchestration.
            </p>
          </div>

          <div>
            <label className="eyebrow" htmlFor="molecule">
              Molecule / asset
            </label>
            <input
              id="molecule"
              className="input"
              value={molecule}
              onChange={(event) => setMolecule(event.target.value)}
              required
            />
          </div>

          <div className="grid-two">
            <div>
              <label className="eyebrow" htmlFor="smiles">
                SMILES
              </label>
              <input
                id="smiles"
                className="input"
                value={smiles}
                onChange={(event) => setSmiles(event.target.value)}
                placeholder="O=C1NC(=O)N(C)C=C1C"
              />
            </div>
            <div>
              <label className="eyebrow" htmlFor="inchikey">
                InChI key
              </label>
              <input
                id="inchikey"
                className="input"
                value={inchikey}
                onChange={(event) => setInchikey(event.target.value)}
                placeholder="XXXXXXXXXXXXXX"
              />
            </div>
          </div>

          <div>
            <label className="eyebrow" htmlFor="synonyms">
              Synonyms / aliases
            </label>
            <textarea
              id="synonyms"
              className="textarea"
              value={synonyms}
              onChange={(event) => setSynonyms(event.target.value)}
            />
          </div>

          <div className="grid-two">
            {(["quick", "full"] as const).map((value) => (
              <label key={value} className={`chip ${depth === value ? "tab--active" : ""}`}>
                <input
                  type="radio"
                  name="depth"
                  className="mr-2"
                  checked={depth === value}
                  onChange={() => setDepth(value)}
                />
                {value === "quick" ? "Quick scan" : "Full evaluation"}
              </label>
            ))}
          </div>

          <div>
            <label className="eyebrow" htmlFor="note">
              Analyst note
            </label>
            <textarea
              id="note"
              className="textarea"
              value={note}
              onChange={(event) => setNote(event.target.value)}
              placeholder="Call out safety watch-outs, preferred datasets, etc."
            />
          </div>

          <div className="space-y-2">
            <button className="btn-primary w-full justify-center" type="submit" disabled={loading}>
              {loading ? "Submitting..." : "Run analysis"}
            </button>
            <label className="btn-secondary w-full cursor-pointer justify-center text-sm">
              Upload context (.json/.txt)
              <input type="file" className="hidden" accept=".json,.txt,.csv" />
            </label>
          </div>

          {statusMessage && <p className="subtle-text">{statusMessage}</p>}
          {jobId && <p className="eyebrow">Tracking job #{jobId}</p>}
          {error && <p className="text-rose-300">{error}</p>}
        </form>

        <div className="section-card space-y-4">
          <p className="eyebrow">Agent pipeline</p>
          <JobTimeline steps={JOB_TIMELINE} />
          <div className="chart-placeholder">LLM orchestration, retrieval, and QA status feed will render here.</div>
        </div>
      </div>
    </div>
  );
}
