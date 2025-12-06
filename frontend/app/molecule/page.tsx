"use client";

import { useCallback, useState } from "react";
import { useRouter } from "next/navigation";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from "@/components/ui/card";
import { AlertCircle, Beaker, CheckCircle2, ChevronDown, ChevronUp, Loader2 } from "lucide-react";
import { api } from "@/lib/api";
import { jobStore } from "@/lib/store";
import ProtectedRoute from "@/components/ProtectedRoute";

function MoleculePage() {
  const router = useRouter();
  const [molecule, setMolecule] = useState("Pirfenidone");
  const [smiles, setSmiles] = useState("");
  const [query, setQuery] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showAdvanced, setShowAdvanced] = useState(false);

  const submit = useCallback(
    async (event: React.FormEvent<HTMLFormElement>) => {
      event.preventDefault();
      setLoading(true);
      setError(null);

      try {
        const job = await api.createJob({
          smiles: smiles || '',
          molecule: molecule || undefined,
          query: query || undefined,
        });
        
        // Save to local storage for history
        jobStore.add({
          job_id: job.job_id,
          smiles: smiles,
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
    <div className="flex justify-center items-center min-h-[calc(100vh-100px)]">
      <Card className="w-full max-w-2xl">
        <CardHeader className="space-y-3">
          <div className="h-1 w-8 bg-primary rounded-full" />
          <CardTitle className="text-2xl flex items-center gap-2">
            <Beaker className="h-6 w-6 text-primary" />
            New Molecule Analysis
          </CardTitle>
          <CardDescription className="leading-relaxed">
            Configure the agents to scout for repurposing opportunities.
          </CardDescription>
        </CardHeader>
        <form onSubmit={submit}>
          <CardContent className="space-y-6">
            <div className="space-y-2">
              <Label htmlFor="smiles">SMILES String (Required)</Label>
              <Input
                id="smiles"
                placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                required
                className="font-mono text-sm"
              />
              <p className="text-xs text-muted-foreground">Enter the molecular structure as a SMILES string</p>
            </div>

            <div className="border-t pt-4">
              <button
                type="button"
                className="flex items-center text-sm text-muted-foreground hover:text-foreground transition-colors"
                onClick={() => setShowAdvanced(!showAdvanced)}
              >
                {showAdvanced ? <ChevronUp className="h-4 w-4 mr-1" /> : <ChevronDown className="h-4 w-4 mr-1" />}
                Optional Fields
              </button>
              
              {showAdvanced && (
                <div className="mt-4 space-y-4 animate-in fade-in slide-in-from-top-2">
                  <div className="space-y-2">
                    <Label htmlFor="molecule">Molecule Name</Label>
                    <Input
                      id="molecule"
                      placeholder="e.g. Aspirin"
                      value={molecule}
                      onChange={(e) => setMolecule(e.target.value)}
                      className="text-lg"
                    />
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="query">Custom Query</Label>
                    <Textarea
                      id="query"
                      placeholder="Enter specific questions or focus areas for the analysis..."
                      value={query}
                      onChange={(e) => setQuery(e.target.value)}
                      className="min-h-20"
                    />
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
          <CardFooter className="flex justify-end gap-3 pt-6">
            <Button variant="ghost" type="button" onClick={() => router.back()} size="sm">Cancel</Button>
            <Button type="submit" disabled={loading} className="min-w-[140px]" size="sm">
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

export default function ProtectedMoleculePage() {
  return (
    <ProtectedRoute>
      <MoleculePage />
    </ProtectedRoute>
  );
}
