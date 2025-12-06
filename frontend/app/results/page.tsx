import Link from "next/link";
import { EvidenceTabs } from "../components/EvidenceTabs";
import { DEMO_JOB, MARKET_METRICS, SAMPLE_PAYLOAD } from "../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ArrowRight, BarChart3, BookOpen, FileText, Microscope, Scale, AlertTriangle, CheckCircle2, Info } from "lucide-react";

export default function ResultsPage() {
  const workers = SAMPLE_PAYLOAD.workers;
  const validation = SAMPLE_PAYLOAD.validation;
  const reportVersion = SAMPLE_PAYLOAD.report_version ?? 1;
  const structure = SAMPLE_PAYLOAD.structure;
  const quality = SAMPLE_PAYLOAD.quality;

  const qualityStatusTone: Record<string, string> = {
    pass: "bg-emerald-100 text-emerald-800",
    needs_attention: "bg-amber-100 text-amber-800",
    investigate: "bg-rose-100 text-rose-800"
  };

  return (
    <div className="space-y-8 py-8">
      <div className="grid gap-8 md:grid-cols-3">
        <div className="md:col-span-2 space-y-8">
          <Card>
            <CardHeader>
              <div className="flex items-center justify-between">
                <div className="space-y-1">
                  <CardTitle className="text-2xl">Pirfenidone Overview</CardTitle>
                  <CardDescription>Innovation Story & Recommendation</CardDescription>
                </div>
                <Badge variant={SAMPLE_PAYLOAD.recommendation === "GO" ? "default" : "destructive"} className="text-lg px-4 py-1">
                  {SAMPLE_PAYLOAD.recommendation}
                </Badge>
              </div>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="prose prose-invert max-w-none">
                <p className="text-lg leading-relaxed text-muted-foreground">
                  {SAMPLE_PAYLOAD.innovation_story}
                </p>
              </div>
              
              <div className="flex flex-wrap gap-4 pt-4">
                <Link href="/reports">
                  <Button>
                    <FileText className="mr-2 h-4 w-4" /> Download Report (V{reportVersion})
                  </Button>
                </Link>
                <Link href="/evidence/clinical">
                  <Button variant="secondary">
                    Dive into Evidence <ArrowRight className="ml-2 h-4 w-4" />
                  </Button>
                </Link>
              </div>
            </CardContent>
          </Card>

          <div className="grid gap-4 md:grid-cols-2">
             {Object.entries(workers).map(([key, worker]) => (
                <Card key={key} className="bg-card/50">
                   <CardHeader className="pb-2">
                      <CardTitle className="text-base capitalize flex items-center gap-2">
                         {key === 'clinical' && <Microscope className="h-4 w-4" />}
                         {key === 'literature' && <BookOpen className="h-4 w-4" />}
                         {key === 'patent' && <Scale className="h-4 w-4" />}
                         {key === 'market' && <BarChart3 className="h-4 w-4" />}
                         {key} Analysis
                      </CardTitle>
                   </CardHeader>
                   <CardContent>
                      <p className="text-sm text-muted-foreground line-clamp-3">
                         {worker.summary}
                      </p>
                      <Link href={`/evidence/${key}`} className="text-xs text-primary hover:underline mt-2 inline-block">
                         View details
                      </Link>
                   </CardContent>
                </Card>
             ))}
          </div>
        </div>

        <div className="space-y-8">
          <Card>
            <CardHeader>
              <CardTitle>Key Figures</CardTitle>
            </CardHeader>
            <CardContent className="space-y-6">
              {MARKET_METRICS.map((metric) => (
                <div key={metric.label} className="space-y-1">
                  <p className="text-sm font-medium text-muted-foreground">{metric.label}</p>
                  <p className="text-2xl font-bold">{metric.value}</p>
                  <p className="text-xs text-muted-foreground">{metric.description}</p>
                </div>
              ))}
            </CardContent>
          </Card>

          {structure && (
             <Card>
                <CardHeader>
                   <CardTitle>Structure</CardTitle>
                </CardHeader>
                <CardContent className="flex flex-col items-center justify-center bg-white/5 rounded-md p-4 space-y-4">
                   {structure.svg ? (
                      <div
                        className="w-full h-48 flex items-center justify-center bg-white rounded p-2"
                        dangerouslySetInnerHTML={{ __html: structure.svg }}
                      />
                   ) : (
                      <div className="text-center text-muted-foreground text-sm py-8">
                         Structure Preview Unavailable
                      </div>
                   )}
                   <div className="text-center w-full">
                      <div className="text-xs font-mono break-all opacity-50 bg-black/20 p-2 rounded">
                         {structure.smiles}
                      </div>
                      {structure.source_type && (
                         <p className="text-[10px] text-muted-foreground mt-2">
                            Source: {structure.source_type.toUpperCase()}
                         </p>
                      )}
                   </div>
                </CardContent>
             </Card>
          )}
        </div>
      </div>

      {quality && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-1">
                <CardTitle>Retrieval Health</CardTitle>
                <CardDescription>Quality guardrails and evidence coverage</CardDescription>
              </div>
              <Badge variant={quality.status === "pass" ? "outline" : "destructive"} className={quality.status === "pass" ? "text-green-500 border-green-500" : ""}>
                {quality.status === "pass" ? "Pass" : "Needs Attention"}
              </Badge>
            </div>
          </CardHeader>
          <CardContent className="space-y-6">
            {Object.keys(quality.alerts).length > 0 ? (
              <div className="space-y-3">
                {Object.entries(quality.alerts).map(([worker, alerts]) => (
                  <div key={worker} className="rounded-md border border-amber-500/20 bg-amber-500/10 p-4">
                    <div className="flex items-center gap-2 text-sm font-semibold text-amber-500 mb-2">
                      <AlertTriangle className="h-4 w-4" />
                      <span className="capitalize">{worker} Alerts</span>
                    </div>
                    <ul className="list-disc list-inside text-sm text-muted-foreground">
                      {alerts.map((alert, index) => (
                        <li key={index}>{alert}</li>
                      ))}
                    </ul>
                  </div>
                ))}
              </div>
            ) : (
              <div className="rounded-md border border-green-500/20 bg-green-500/10 p-4 flex items-center gap-2 text-green-500">
                <CheckCircle2 className="h-5 w-5" />
                <span className="text-sm font-medium">All workers met minimum evidence thresholds.</span>
              </div>
            )}

            <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
              {Object.entries(quality.metrics).map(([worker, metrics]) => (
                <div key={worker} className="rounded-lg border bg-card p-4 space-y-3">
                  <div className="flex items-center justify-between">
                    <Badge variant="secondary" className="capitalize">{worker}</Badge>
                    <span className="text-xs text-muted-foreground">{metrics.evidence_count} items</span>
                  </div>
                  <div className="space-y-1">
                    <div className="flex justify-between text-sm">
                      <span className="text-muted-foreground">Coverage</span>
                      <span className="font-medium">{(metrics.coverage_ratio * 100).toFixed(0)}%</span>
                    </div>
                    <div className="flex justify-between text-sm">
                      <span className="text-muted-foreground">Precision</span>
                      <span className="font-medium">{(metrics.precision_proxy * 100).toFixed(0)}%</span>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
      )}
      
      <div className="pt-8">
         <EvidenceTabs workers={workers} reportJobId={DEMO_JOB.id} reportVersion={reportVersion} />
      </div>

      {validation && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-1">
                <CardTitle>Claim Traceability</CardTitle>
                <CardDescription>Validation pass for innovation story claims</CardDescription>
              </div>
              <Badge variant="outline">
                {validation.claims_linked}/{validation.claims_total} claims linked
              </Badge>
            </div>
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              {validation.claim_links.map((claim) => (
                <div key={claim.claim_id} className="rounded-lg border bg-card/50 p-4 space-y-2">
                  <div className="flex items-center justify-between">
                    <Badge variant="secondary" className="capitalize">{claim.worker}</Badge>
                    <Badge variant={claim.status === "linked" ? "outline" : "destructive"} className={claim.status === "linked" ? "text-green-500 border-green-500" : ""}>
                      {claim.status === "linked" ? "Linked" : "Needs Review"}
                    </Badge>
                  </div>
                  <p className="text-sm">{claim.claim_text}</p>
                  <div className="flex items-center gap-2 text-xs text-muted-foreground">
                    <Info className="h-3 w-3" />
                    Evidence IDs: {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                  </div>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
