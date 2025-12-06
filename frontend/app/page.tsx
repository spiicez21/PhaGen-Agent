import Link from "next/link";
import { HISTORY_RUNS } from "./sample-data";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { ArrowRight, Beaker, FileText, History, Search } from "lucide-react";

export default function LandingPage() {
  const recent = HISTORY_RUNS.slice(0, 3);

  return (
    <div className="flex flex-col items-center justify-center space-y-12 py-12 md:py-24 lg:py-32">
      <div className="container flex flex-col items-center text-center gap-6 max-w-3xl">
        <div className="space-y-2">
          <h1 className="text-4xl font-bold tracking-tighter sm:text-5xl md:text-6xl lg:text-7xl bg-clip-text text-transparent bg-linear-to-r from-white to-gray-400">
            PhaGen Agentic
          </h1>
          <p className="mx-auto max-w-[700px] text-muted-foreground md:text-xl">
            Molecule Scouting in Minutes.
            <br />
            Orchestrate clinical, literature, patent, and market agents to accelerate repurposing decisions.
          </p>
        </div>
        <div className="flex flex-col gap-4 min-[400px]:flex-row">
          <Link href="/molecule">
            <Button size="lg" className="h-12 px-8">
              <Search className="mr-2 h-4 w-4" />
              Start New Analysis
            </Button>
          </Link>
          <Link href="/molecule?demo=true">
            <Button variant="outline" size="lg" className="h-12 px-8">
              <Beaker className="mr-2 h-4 w-4" />
              Try Demo Molecule
            </Button>
          </Link>
        </div>
      </div>

      <div className="container max-w-5xl space-y-8">
        <div className="flex items-center justify-between">
          <h2 className="text-2xl font-bold tracking-tight flex items-center gap-2">
            <History className="h-6 w-6" />
            Recent Analyses
          </h2>
          <Link href="/history">
            <Button variant="ghost" className="text-muted-foreground">View all</Button>
          </Link>
        </div>

        <div className="grid gap-4 md:grid-cols-3">
          {recent.map((run) => (
            <Card key={run.id} className="bg-card/50 backdrop-blur transition-all hover:bg-card/80 hover:shadow-md border-muted">
              <CardHeader className="pb-2">
                <div className="flex justify-between items-start">
                  <CardTitle className="text-lg font-medium">{run.molecule}</CardTitle>
                  <Badge variant={run.recommendation === "GO" ? "default" : run.recommendation === "NO-GO" ? "destructive" : "secondary"}>
                    {run.recommendation}
                  </Badge>
                </div>
                <CardDescription>{run.date}</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="flex justify-between items-center mt-4">
                  <span className="text-sm text-muted-foreground">Status: {run.status}</span>
                  <Link href={`/results?id=${run.id}`}>
                    <Button variant="ghost" size="sm" className="gap-1">
                      Open <ArrowRight className="h-3 w-3" />
                    </Button>
                  </Link>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>
      
      <div className="container max-w-5xl grid gap-8 md:grid-cols-3 pt-12">
         <div className="flex flex-col items-center text-center space-y-2 p-6 rounded-lg border bg-card/30">
            <div className="p-3 rounded-full bg-primary/10 text-primary">
               <FileText className="h-6 w-6" />
            </div>
            <h3 className="text-xl font-bold">Comprehensive Reports</h3>
            <p className="text-sm text-muted-foreground">Generate full PDF reports with citations, evidence tables, and market analysis.</p>
         </div>
         <div className="flex flex-col items-center text-center space-y-2 p-6 rounded-lg border bg-card/30">
            <div className="p-3 rounded-full bg-primary/10 text-primary">
               <Beaker className="h-6 w-6" />
            </div>
            <h3 className="text-xl font-bold">Multi-Agent AI</h3>
            <p className="text-sm text-muted-foreground">Specialized agents for clinical trials, literature, patents, and market data.</p>
         </div>
         <div className="flex flex-col items-center text-center space-y-2 p-6 rounded-lg border bg-card/30">
            <div className="p-3 rounded-full bg-primary/10 text-primary">
               <Search className="h-6 w-6" />
            </div>
            <h3 className="text-xl font-bold">Deep Search</h3>
            <p className="text-sm text-muted-foreground">RAG-powered retrieval across millions of documents and trial records.</p>
         </div>
      </div>
    </div>
  );
}
