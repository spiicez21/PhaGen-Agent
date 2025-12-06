import { HISTORY_RUNS } from "../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { FileText, ExternalLink, Clock, CheckCircle2, AlertCircle } from "lucide-react";
import Link from "next/link";

export default function HistoryPage() {
  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2 text-sm text-muted-foreground mb-2">
            <span className="px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">Workspace</span>
            <span>Archives</span>
          </div>
          <CardTitle className="text-2xl">Saved Reports</CardTitle>
          <CardDescription>
            Access your past analysis runs, download PDF reports, and review recommendations.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="rounded-md border">
            <table className="w-full text-sm text-left">
              <thead className="bg-muted/50 text-muted-foreground font-medium">
                <tr>
                  <th className="h-12 px-4 align-middle">Molecule</th>
                  <th className="h-12 px-4 align-middle">Date</th>
                  <th className="h-12 px-4 align-middle">Recommendation</th>
                  <th className="h-12 px-4 align-middle">Status</th>
                  <th className="h-12 px-4 align-middle text-right">Actions</th>
                </tr>
              </thead>
              <tbody className="divide-y">
                {HISTORY_RUNS.map((run, i) => (
                  <tr key={`${run.molecule}-${run.date}-${i}`} className="hover:bg-muted/50 transition-colors">
                    <td className="p-4 font-medium">{run.molecule}</td>
                    <td className="p-4 text-muted-foreground">
                      <div className="flex items-center gap-2">
                        <Clock className="h-3 w-3" /> {run.date}
                      </div>
                    </td>
                    <td className="p-4">
                      <Badge variant={run.recommendation === "GO" ? "default" : "destructive"}>
                        {run.recommendation}
                      </Badge>
                    </td>
                    <td className="p-4">
                      <div className="flex items-center gap-2">
                        {run.status === "Complete" ? (
                          <CheckCircle2 className="h-4 w-4 text-green-500" />
                        ) : (
                          <AlertCircle className="h-4 w-4 text-amber-500" />
                        )}
                        <span>{run.status}</span>
                      </div>
                    </td>
                    <td className="p-4 text-right">
                      <div className="flex items-center justify-end gap-2">
                        <Button variant="ghost" size="sm" asChild>
                          <Link href="/results">
                            <ExternalLink className="h-4 w-4 mr-2" /> Open
                          </Link>
                        </Button>
                        <Button variant="outline" size="sm">
                          <FileText className="h-4 w-4 mr-2" /> PDF
                        </Button>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
