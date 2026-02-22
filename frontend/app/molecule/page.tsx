"use client";

import { useCallback, useState } from "react";
import { useRouter } from "next/navigation";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import {
  AlertCircle, ChevronDown, ChevronUp,
  FlaskConical, Loader2, Sparkles, Info
} from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

const EXAMPLE_SMILES = [
  { name: "Pirfenidone",  smiles: "O=C1C=CC=CN1c1ccc(C)cc1" },
  { name: "Aspirin",      smiles: "CC(=O)Oc1ccccc1C(=O)O" },
  { name: "Metformin",    smiles: "CN(C)C(=N)NC(=N)N" },
];

function MoleculePage() {
  const router = useRouter();
  const [molecule, setMolecule]         = useState("Pirfenidone");
  const [smiles, setSmiles]             = useState("");
  const [query, setQuery]               = useState("");
  const [loading, setLoading]           = useState(false);
  const [error, setError]               = useState<string | null>(null);
  const [showAdvanced, setShowAdvanced] = useState(false);

  const fill = (name: string, s: string) => {
    setMolecule(name);
    setSmiles(s);
  };

  const submit = useCallback(
    async (event: React.FormEvent<HTMLFormElement>) => {
      event.preventDefault();
      setLoading(true);
      setError(null);
      try {
        const job = await api.createJob({
          smiles: smiles || "",
          molecule: molecule || undefined,
          query: query || undefined,
        });
        jobStore.add({
          job_id: job.job_id,
          smiles,
          molecule: molecule || undefined,
          status: job.status,
          created_at: job.created_at,
        });
        router.push(`/job?id=${job.job_id}`);
      } catch (err) {
        setError(err instanceof Error ? err.message : "An unknown error occurred");
        setLoading(false);
      }
    },
    [molecule, smiles, query, router]
  );

  return (
    <div className="max-w-4xl mx-auto space-y-4">

      {/* Page header */}
      <div>
        <h1 className="text-base font-semibold text-foreground">New Analysis</h1>
        <p className="text-[11px] text-muted-foreground mt-0.5">
          Configure and launch a multi-agent molecule scouting job
        </p>
      </div>

      <div className="grid lg:grid-cols-5 gap-4">

        {/* Form */}
        <div className="lg:col-span-3 card-flush overflow-hidden">
          <div className="px-4 py-2.5 border-b border-border flex items-center gap-2">
            <FlaskConical className="h-3.5 w-3.5 text-primary" />
            <h2 className="text-xs font-semibold">Molecule Input</h2>
          </div>

          <form onSubmit={submit}>
            <div className="p-4 space-y-4">

              {/* SMILES */}
              <div className="space-y-1.5">
                <Label htmlFor="smiles" className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wide">
                  SMILES String <span className="text-destructive">*</span>
                </Label>
                <Input
                  id="smiles"
                  placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  required
                  className="font-mono text-xs h-9 bg-muted/30 border-border"
                />
                <p className="text-[10px] text-muted-foreground">
                  SMILES notation for the molecule structure
                </p>
              </div>

              {/* Advanced toggle */}
              <div>
                <button
                  type="button"
                  className="flex items-center gap-1.5 text-[11px] text-muted-foreground hover:text-foreground transition-colors py-1"
                  onClick={() => setShowAdvanced(!showAdvanced)}
                >
                  {showAdvanced ? <ChevronUp className="h-3 w-3" /> : <ChevronDown className="h-3 w-3" />}
                  <span className="font-medium">Optional Parameters</span>
                </button>

                {showAdvanced && (
                  <div className="mt-3 space-y-3 pl-4 border-l-2 border-border animate-in fade-in slide-in-from-top-1 duration-150">
                    <div className="space-y-1.5">
                      <Label htmlFor="molecule" className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wide">
                        Molecule Name
                      </Label>
                      <Input
                        id="molecule"
                        placeholder="e.g. Aspirin"
                        value={molecule}
                        onChange={(e) => setMolecule(e.target.value)}
                        className="h-9 bg-muted/30 border-border text-xs"
                      />
                    </div>
                    <div className="space-y-1.5">
                      <Label htmlFor="query" className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wide">
                        Custom Query
                      </Label>
                      <Textarea
                        id="query"
                        placeholder="Specific focus areas or questions for the analysis..."
                        value={query}
                        onChange={(e) => setQuery(e.target.value)}
                        className="min-h-16 bg-muted/30 border-border text-xs resize-none"
                      />
                    </div>
                  </div>
                )}
              </div>

              {/* Error */}
              {error && (
                <div className="flex items-start gap-2 p-3 rounded bg-destructive/8 border border-destructive/20 text-destructive text-xs">
                  <AlertCircle className="h-3.5 w-3.5 mt-0.5 shrink-0" />
                  <span>{error}</span>
                </div>
              )}
            </div>

            {/* Footer */}
            <div className="flex items-center justify-between px-4 py-3 border-t border-border bg-muted/20">
              <button
                type="button"
                onClick={() => router.back()}
                className="text-[11px] text-muted-foreground hover:text-foreground transition-colors"
              >
                Cancel
              </button>
              <Button
                type="submit"
                disabled={loading || !smiles.trim()}
                size="sm"
                className="h-7 px-4 text-xs bg-primary text-primary-foreground disabled:opacity-50"
              >
                {loading ? (
                  <>
                    <Loader2 className="mr-1.5 h-3 w-3 animate-spin" />
                    Submitting...
                  </>
                ) : (
                  <>
                    <Sparkles className="mr-1.5 h-3 w-3" />
                    Run Analysis
                  </>
                )}
              </Button>
            </div>
          </form>
        </div>

        {/* Right panel: examples + agent info */}
        <div className="lg:col-span-2 space-y-4">

          {/* Example molecules */}
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border">
              <h3 className="text-xs font-semibold">Example Molecules</h3>
            </div>
            <div className="divide-y divide-border">
              {EXAMPLE_SMILES.map(({ name, smiles: s }) => (
                <button
                  key={name}
                  type="button"
                  onClick={() => fill(name, s)}
                  className="w-full text-left px-4 py-2.5 hover:bg-muted/40 transition-colors group"
                >
                  <p className="text-xs font-semibold text-foreground group-hover:text-primary transition-colors">{name}</p>
                  <p className="font-mono text-[10px] text-muted-foreground truncate mt-0.5">{s}</p>
                </button>
              ))}
            </div>
          </div>

          {/* Agent pipeline info */}
          <div className="card-flush overflow-hidden">
            <div className="px-4 py-2.5 border-b border-border flex items-center gap-2">
              <Info className="h-3 w-3 text-muted-foreground" />
              <h3 className="text-xs font-semibold">Agent Pipeline</h3>
            </div>
            <div className="divide-y divide-border">
              {[
                { name: "Clinical",    color: "bg-blue-500",   desc: "ClinicalTrials.gov search" },
                { name: "Literature",  color: "bg-purple-500", desc: "PubMed + Semantic Scholar" },
                { name: "Patent",      color: "bg-amber-500",  desc: "Patent database search" },
                { name: "Market",      color: "bg-teal-500",   desc: "Market & competitive analysis" },
              ].map(({ name, color, desc }) => (
                <div key={name} className="flex items-center gap-3 px-4 py-2.5">
                  <div className={`h-1.5 w-1.5 rounded-full ${color} flex-shrink-0`} />
                  <div>
                    <span className="text-xs font-medium text-foreground">{name}</span>
                    <span className="text-[10px] text-muted-foreground ml-2">{desc}</span>
                  </div>
                </div>
              ))}
            </div>
            <div className="px-4 py-2.5 bg-muted/20 border-t border-border">
              <p className="text-[10px] text-muted-foreground leading-relaxed">
                All 4 agents run in parallel. Results are synthesized into a GO / INVESTIGATE / NO-GO recommendation.
              </p>
            </div>
          </div>

        </div>
      </div>
    </div>
  );
}

export default function ProtectedMoleculePage() {
  return (
    <ProtectedRoute>
      <MoleculePage />
    </ProtectedRoute>
  );
}
