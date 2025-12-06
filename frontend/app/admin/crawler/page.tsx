import { CRAWLER_QUEUE, ROBOTS_STATUS } from "../../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Activity, CheckCircle2, XCircle, Plus, Globe, ShieldAlert, ShieldCheck } from "lucide-react";

export default function CrawlerStatusPage() {
  const stats = {
    active: 4,
    completed: 320,
    failed: 2
  };

  return (
    <div className="space-y-8 py-8">
      <div className="grid gap-4 md:grid-cols-3">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Active Crawls</CardTitle>
            <Activity className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{stats.active}</div>
            <p className="text-xs text-muted-foreground">Currently processing</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Completed</CardTitle>
            <CheckCircle2 className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{stats.completed}</div>
            <p className="text-xs text-muted-foreground">Pages indexed</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
            <CardTitle className="text-sm font-medium">Failed</CardTitle>
            <XCircle className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{stats.failed}</div>
            <p className="text-xs text-muted-foreground">Errors encountered</p>
          </CardContent>
        </Card>
      </div>

      <div className="grid gap-8 md:grid-cols-2">
        <Card className="h-full">
          <CardHeader className="flex flex-row items-center justify-between">
            <div>
              <CardTitle>Crawler Queue</CardTitle>
              <CardDescription>Pending and active URLs</CardDescription>
            </div>
            <Button size="sm" variant="outline">
              <Plus className="mr-2 h-4 w-4" /> Add URL
            </Button>
          </CardHeader>
          <CardContent>
            <div className="rounded-md border">
              <table className="w-full text-sm text-left">
                <thead className="bg-muted/50 text-muted-foreground font-medium">
                  <tr>
                    <th className="h-10 px-4 align-middle">URL</th>
                    <th className="h-10 px-4 align-middle">Status</th>
                    <th className="h-10 px-4 align-middle text-right">Retries</th>
                  </tr>
                </thead>
                <tbody className="divide-y">
                  {CRAWLER_QUEUE.map((item, i) => (
                    <tr key={`${item.url}-${i}`} className="hover:bg-muted/50">
                      <td className="p-4 truncate max-w-[200px]" title={item.url}>
                        <div className="flex items-center gap-2">
                          <Globe className="h-3 w-3 text-muted-foreground" />
                          {item.url}
                        </div>
                      </td>
                      <td className="p-4">
                        <Badge variant="secondary" className="capitalize">
                          {item.status}
                        </Badge>
                      </td>
                      <td className="p-4 text-right font-mono">{item.retries}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>

        <Card className="h-full">
          <CardHeader>
            <CardTitle>Robots.txt Compliance</CardTitle>
            <CardDescription>Domain access permissions</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="rounded-md border">
              <table className="w-full text-sm text-left">
                <thead className="bg-muted/50 text-muted-foreground font-medium">
                  <tr>
                    <th className="h-10 px-4 align-middle">Domain</th>
                    <th className="h-10 px-4 align-middle">Access</th>
                    <th className="h-10 px-4 align-middle">Note</th>
                  </tr>
                </thead>
                <tbody className="divide-y">
                  {ROBOTS_STATUS.map((row, i) => (
                    <tr key={`${row.domain}-${i}`} className="hover:bg-muted/50">
                      <td className="p-4 font-medium">{row.domain}</td>
                      <td className="p-4">
                        <div className="flex items-center gap-2">
                          {row.access === "Allowed" ? (
                            <ShieldCheck className="h-4 w-4 text-green-500" />
                          ) : (
                            <ShieldAlert className="h-4 w-4 text-amber-500" />
                          )}
                          <span>{row.access}</span>
                        </div>
                      </td>
                      <td className="p-4 text-muted-foreground text-xs">{row.note}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
