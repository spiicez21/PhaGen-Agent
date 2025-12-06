"use client";

import { useCallback, useState } from "react";
import { useRouter } from "next/navigation";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from "@/components/ui/card";
import { AlertCircle, Beaker, CheckCircle2, ChevronDown, ChevronUp, Loader2, Upload } from "lucide-react";
import { cn } from "@/lib/utils";

const API_BASE = process.env.NEXT_PUBLIC_API_URL ?? "http://localhost:8000";

type JobResponse = {
  job_id: string;
  status: string;
};

export default function MoleculePage() {
  const router = useRouter();
  const [molecule, setMolecule] = useState("Pirfenidone");
  const [synonyms, setSynonyms] = useState("Esbriet\nPFD");
  const [smiles, setSmiles] = useState("");
  const [inchikey, setInchikey] = useState("KUFJONJOBWVSNK-UHFFFAOYSA-N");
  const [depth, setDepth] = useState<"quick" | "full">("full");
  const [note, setNote] = useState("Focus on PF-ILD cohort.");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showAdvanced, setShowAdvanced] = useState(false);

  const submit = useCallback(
    async (event: React.FormEvent<HTMLFormElement>) => {
      event.preventDefault();
      setLoading(true);
      setError(null);
      
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
        router.push(`/job?id=${payload.job_id}`);
      } catch (err) {
        setError(err instanceof Error ? err.message : "An unknown error occurred");
        setLoading(false);
      }
    },
    [molecule, synonyms, smiles, inchikey, depth, note, router]
  );

  return (
    <div className="flex justify-center items-center min-h-[calc(100vh-100px)]">
      <Card className="w-full max-w-2xl shadow-lg border-muted/40">
        <CardHeader>
          <CardTitle className="text-2xl flex items-center gap-2">
            <Beaker className="h-6 w-6 text-primary" />
            New Molecule Analysis
          </CardTitle>
          <CardDescription>
            Configure the agents to scout for repurposing opportunities.
          </CardDescription>
        </CardHeader>
        <form onSubmit={submit}>
          <CardContent className="space-y-6">
            <div className="space-y-2">
              <Label htmlFor="molecule">Molecule Name</Label>
              <Input
                id="molecule"
                placeholder="e.g. Aspirin"
                value={molecule}
                onChange={(e) => setMolecule(e.target.value)}
                required
                className="text-lg"
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="synonyms">Synonyms (one per line)</Label>
              <Textarea
                id="synonyms"
                placeholder="Enter synonyms..."
                value={synonyms}
                onChange={(e) => setSynonyms(e.target.value)}
                className="min-h-20"
              />
            </div>

            <div className="space-y-2">
              <Label>Analysis Depth</Label>
              <div className="grid grid-cols-2 gap-4">
                <div
                  className={cn(
                    "cursor-pointer rounded-lg border-2 p-4 hover:bg-accent transition-all",
                    depth === "quick" ? "border-primary bg-primary/5" : "border-muted"
                  )}
                  onClick={() => setDepth("quick")}
                >
                  <div className="font-semibold">Quick Scan</div>
                  <div className="text-xs text-muted-foreground">Rapid feasibility check</div>
                </div>
                <div
                  className={cn(
                    "cursor-pointer rounded-lg border-2 p-4 hover:bg-accent transition-all",
                    depth === "full" ? "border-primary bg-primary/5" : "border-muted"
                  )}
                  onClick={() => setDepth("full")}
                >
                  <div className="font-semibold">Full Evaluation</div>
                  <div className="text-xs text-muted-foreground">Deep multi-agent report</div>
                </div>
              </div>
            </div>

            <div className="border-t pt-4">
              <button
                type="button"
                className="flex items-center text-sm text-muted-foreground hover:text-foreground transition-colors"
                onClick={() => setShowAdvanced(!showAdvanced)}
              >
                {showAdvanced ? <ChevronUp className="h-4 w-4 mr-1" /> : <ChevronDown className="h-4 w-4 mr-1" />}
                Advanced Inputs
              </button>
              
              {showAdvanced && (
                <div className="mt-4 space-y-4 animate-in fade-in slide-in-from-top-2">
                  <div className="space-y-2">
                    <Label htmlFor="smiles">SMILES String</Label>
                    <Input
                      id="smiles"
                      placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
                      value={smiles}
                      onChange={(e) => setSmiles(e.target.value)}
                      className="font-mono text-xs"
                    />
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="inchikey">InChI Key</Label>
                    <Input
                      id="inchikey"
                      placeholder="e.g. BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
                      value={inchikey}
                      onChange={(e) => setInchikey(e.target.value)}
                      className="font-mono text-xs"
                    />
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="note">Research Note</Label>
                    <Textarea
                      id="note"
                      placeholder="Add context for the agents..."
                      value={note}
                      onChange={(e) => setNote(e.target.value)}
                    />
                  </div>
                  <div className="pt-2">
                     <Button variant="outline" type="button" className="w-full border-dashed">
                        <Upload className="mr-2 h-4 w-4" /> Upload Structure File (.json/.txt)
                     </Button>
                  </div>
                </div>
              )}
            </div>

            {error && (
              <div className="p-3 rounded-md bg-destructive/10 text-destructive text-sm flex items-center">
                <AlertCircle className="h-4 w-4 mr-2" />
                {error}
              </div>
            )}
          </CardContent>
          <CardFooter className="flex justify-end gap-4">
            <Button variant="ghost" type="button" onClick={() => router.back()}>Cancel</Button>
            <Button type="submit" disabled={loading} className="min-w-[140px]">
              {loading ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Analyzing...
                </>
              ) : (
                <>
                  Run Analysis
                  <CheckCircle2 className="ml-2 h-4 w-4" />
                </>
              )}
            </Button>
          </CardFooter>
        </form>
      </Card>
    </div>
  );
}
