"use client";

import { useState, useEffect } from "react";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { ArrowRight, Beaker, FileText, Search } from "lucide-react";
import { jobStore, StoredJob } from "@/lib/store";

export default function LandingPage() {
  const [recent, setRecent] = useState<StoredJob[]>([]);

  useEffect(() => {
    // Load recent jobs client-side only
    const timer = setTimeout(() => setRecent(jobStore.getAll().slice(0, 3)), 0);
    return () => clearTimeout(timer);
  }, []);

  return (
    <div className="flex flex-col space-y-16 py-16 md:py-24">
      <div className="container flex flex-col items-center text-center gap-8 max-w-4xl mx-auto">
        <div className="space-y-4">
          <div className="inline-block rounded-lg bg-muted px-3 py-1 text-sm font-medium mb-4">
            Pharmaceutical Intelligence
          </div>
          <h1 className="text-5xl font-bold tracking-tight sm:text-6xl md:text-7xl">
            Molecule Scouting
            <span className="block text-muted-foreground mt-2">in Minutes</span>
          </h1>
          <p className="mx-auto max-w-[600px] text-muted-foreground text-lg leading-relaxed">
            Orchestrate clinical, literature, patent, and market agents to accelerate repurposing decisions with confidence.
          </p>
        </div>
        <div className="flex flex-col sm:flex-row gap-4">
          <Link href="/molecule">
            <Button size="lg" className="h-12 px-8 font-medium">
              <Search className="mr-2 h-5 w-5" />
              Start Analysis
            </Button>
          </Link>
          <Link href="/molecule?demo=true">
            <Button variant="outline" size="lg" className="h-12 px-8 font-medium">
              <Beaker className="mr-2 h-5 w-5" />
              Try Demo
            </Button>
          </Link>
        </div>
      </div>

      <div className="container max-w-6xl space-y-6">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="h-1 w-8 bg-primary rounded-full" />
            <h2 className="text-2xl font-semibold tracking-tight">Recent Analyses</h2>
          </div>
          <Link href="/history">
            <Button variant="ghost" size="sm" className="text-muted-foreground">
              View All <ArrowRight className="ml-2 h-4 w-4" />
            </Button>
          </Link>
        </div>

        {recent.length === 0 ? (
          <Card className="border-dashed">
            <CardContent className="flex flex-col items-center justify-center py-12 text-center">
              <FileText className="h-12 w-12 text-muted-foreground mb-4" />
              <p className="text-muted-foreground text-sm">No recent analyses yet</p>
              <p className="text-muted-foreground text-xs mt-1">Start your first analysis to see results here</p>
            </CardContent>
          </Card>
        ) : (
          <div className="grid gap-4 md:grid-cols-3">
            {recent.map((job) => (
              <Card key={job.job_id}>
                <CardHeader className="pb-3">
                  <div className="flex justify-between items-start gap-2">
                    <CardTitle className="text-base font-semibold">{job.molecule || "Unknown"}</CardTitle>
                    {job.recommendation && (
                      <Badge variant={job.recommendation === "GO" ? "default" : "destructive"} className="shrink-0">
                        {job.recommendation}
                      </Badge>
                    )}
                  </div>
                  <CardDescription className="text-xs font-mono truncate">{job.smiles}</CardDescription>
                </CardHeader>
                <CardContent>
                  <div className="flex justify-between items-center mt-2">
                    <Badge variant="outline" className="text-xs capitalize font-normal">{job.status}</Badge>
                    <Link href={job.status === "completed" ? `/results?id=${job.job_id}` : `/job?id=${job.job_id}`}>
                      <Button variant="ghost" size="sm" className="h-8 gap-1 text-xs">
                        View <ArrowRight className="h-3 w-3" />
                      </Button>
                    </Link>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )}
      </div>
      
      <div className="container max-w-5xl grid gap-6 md:grid-cols-3 pt-8">
         <Card className="border-none shadow-none bg-muted/30">
            <CardContent className="flex flex-col items-center text-center space-y-3 pt-6">
               <div className="p-3 rounded-full bg-primary text-primary-foreground">
                  <FileText className="h-6 w-6" />
               </div>
               <h3 className="text-lg font-semibold">Reports</h3>
               <p className="text-sm text-muted-foreground leading-relaxed">Generate PDF reports with citations and evidence tables</p>
            </CardContent>
         </Card>
         <Card className="border-none shadow-none bg-muted/30">
            <CardContent className="flex flex-col items-center text-center space-y-3 pt-6">
               <div className="p-3 rounded-full bg-primary text-primary-foreground">
                  <Beaker className="h-6 w-6" />
               </div>
               <h3 className="text-lg font-semibold">Multi-Agent AI</h3>
               <p className="text-sm text-muted-foreground leading-relaxed">Specialized agents for clinical, literature, patents, and market</p>
            </CardContent>
         </Card>
         <Card className="border-none shadow-none bg-muted/30">
            <CardContent className="flex flex-col items-center text-center space-y-3 pt-6">
               <div className="p-3 rounded-full bg-primary text-primary-foreground">
                  <Search className="h-6 w-6" />
               </div>
               <h3 className="text-lg font-semibold">Deep Search</h3>
               <p className="text-sm text-muted-foreground leading-relaxed">RAG-powered retrieval across millions of documents</p>
            </CardContent>
         </Card>
      </div>
    </div>
  );
}
