import { CLINICAL_TRIALS } from "../../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ExternalLink, Filter, Microscope, Activity, Users } from "lucide-react";

export default function ClinicalEvidencePage() {
  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2 text-sm text-muted-foreground mb-2">
            <span className="px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">Evidence Dashboard</span>
            <span>Clinical</span>
          </div>
          <CardTitle className="text-2xl">Clinical Trials</CardTitle>
          <CardDescription>
            Filtered view of trial and registry evidence powering the recommendation.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="flex flex-wrap gap-2 pb-4 border-b">
            <div className="flex items-center gap-2 mr-2 text-sm font-medium text-muted-foreground">
              <Filter className="h-4 w-4" /> Filters:
            </div>
            <Badge variant="secondary">Phase: All</Badge>
            <Badge variant="secondary">Status: Recruiting + Completed</Badge>
            <Badge variant="secondary">Condition: PF-ILD</Badge>
          </div>

          <div className="rounded-md border">
            <table className="w-full text-sm text-left">
              <thead className="bg-muted/50 text-muted-foreground font-medium">
                <tr>
                  <th className="h-10 px-4 align-middle">NCT ID</th>
                  <th className="h-10 px-4 align-middle">Phase</th>
                  <th className="h-10 px-4 align-middle">Status</th>
                  <th className="h-10 px-4 align-middle">Condition</th>
                  <th className="h-10 px-4 align-middle">Outcome</th>
                </tr>
              </thead>
              <tbody className="divide-y">
                {CLINICAL_TRIALS.map((trial, i) => (
                  <tr key={`${trial.nctId}-${i}`} className="hover:bg-muted/50">
                    <td className="p-4 font-mono font-medium text-primary">{trial.nctId}</td>
                    <td className="p-4">{trial.phase}</td>
                    <td className="p-4">
                      <Badge variant={trial.status === "Completed" ? "default" : "outline"}>
                        {trial.status}
                      </Badge>
                    </td>
                    <td className="p-4">{trial.condition}</td>
                    <td className="p-4 text-muted-foreground">{trial.outcome}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>

      <div className="grid gap-6 md:grid-cols-2">
        {CLINICAL_TRIALS.map((trial) => (
          <Card key={`${trial.nctId}-card`} className="flex flex-col">
            <CardHeader>
              <div className="flex items-center justify-between mb-2">
                <Badge variant="outline" className="font-mono">{trial.nctId}</Badge>
                <Badge variant={trial.status === "Completed" ? "default" : "secondary"}>{trial.status}</Badge>
              </div>
              <CardTitle className="text-lg">{trial.condition}</CardTitle>
              <CardDescription>Phase {trial.phase}</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4 flex-1">
              <div className="grid grid-cols-2 gap-4 text-sm">
                <div className="space-y-1">
                  <div className="flex items-center gap-2 text-muted-foreground">
                    <Users className="h-3 w-3" /> Population
                  </div>
                  <p>{trial.population}</p>
                </div>
                <div className="space-y-1">
                  <div className="flex items-center gap-2 text-muted-foreground">
                    <Activity className="h-3 w-3" /> Sample Size
                  </div>
                  <p>{trial.sampleSize}</p>
                </div>
              </div>
              
              <div className="pt-4 border-t">
                <p className="text-sm font-medium mb-1">Primary Outcome</p>
                <p className="text-sm text-muted-foreground">{trial.outcome}</p>
              </div>
              
              <div className="pt-2">
                <Button variant="link" className="h-auto p-0 text-xs" asChild>
                  <a href={trial.source} target="_blank" rel="noreferrer">
                    View Source <ExternalLink className="ml-1 h-3 w-3" />
                  </a>
                </Button>
              </div>
            </CardContent>
          </Card>
        ))}
      </div>
    </div>
  );
}
