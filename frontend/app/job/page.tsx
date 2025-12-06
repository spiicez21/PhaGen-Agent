import Link from "next/link";
import { JobTimeline } from "../components/JobTimeline";
import { DEMO_JOB, JOB_LOGS, JOB_TIMELINE } from "../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ArrowRight, FileText, Terminal } from "lucide-react";

export default function JobStatusPage() {
  return (
    <div className="space-y-8 py-8">
      <div className="grid gap-8 md:grid-cols-2">
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div className="space-y-1">
                <CardTitle>Job {DEMO_JOB.id}</CardTitle>
                <CardDescription>Orchestration Status</CardDescription>
              </div>
              <Badge variant={DEMO_JOB.status === "RUNNING" ? "default" : "secondary"}>
                {DEMO_JOB.status}
              </Badge>
            </div>
          </CardHeader>
          <CardContent className="space-y-6">
            <dl className="grid grid-cols-2 gap-4 text-sm">
              <div>
                <dt className="text-muted-foreground">Depth</dt>
                <dd className="font-medium capitalize">{DEMO_JOB.depth}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Started (UTC)</dt>
                <dd className="font-medium">{new Date(DEMO_JOB.startedAt).toLocaleTimeString()}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">ETA</dt>
                <dd className="font-medium">{DEMO_JOB.etaMinutes} min</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Molecule</dt>
                <dd className="font-medium">Pirfenidone</dd>
              </div>
            </dl>
            
            <div className="flex flex-col gap-3 pt-4">
              <Link href="/results">
                <Button className="w-full">
                  View Provisional Results <ArrowRight className="ml-2 h-4 w-4" />
                </Button>
              </Link>
              <Link href="/reports">
                <Button variant="outline" className="w-full">
                  <FileText className="mr-2 h-4 w-4" /> Prepare Report
                </Button>
              </Link>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardContent className="pt-6">
            <JobTimeline steps={JOB_TIMELINE} />
          </CardContent>
        </Card>
      </div>

      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <Terminal className="h-5 w-5 text-muted-foreground" />
              <CardTitle>Live Logs</CardTitle>
            </div>
            <Badge variant="outline" className="animate-pulse text-green-500 border-green-500/50">
              Streaming
            </Badge>
          </div>
        </CardHeader>
        <CardContent>
          <div className="rounded-md bg-black/50 border border-border p-4 font-mono text-xs h-[300px] overflow-y-auto space-y-1">
            {JOB_LOGS.map((log, i) => (
              <div key={i} className="flex gap-3 text-muted-foreground hover:text-foreground transition-colors">
                <span className="text-muted-foreground/50 shrink-0">{new Date().toLocaleTimeString()}</span>
                <span>{typeof log === 'string' ? log : JSON.stringify(log)}</span>
              </div>
            ))}
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
