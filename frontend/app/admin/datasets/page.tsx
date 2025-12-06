import { INDEX_VERSIONS } from "../../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Database, RefreshCw, Trash2, HardDrive, CheckCircle2, Clock } from "lucide-react";

export default function DatasetManagerPage() {
  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2 text-sm text-muted-foreground mb-2">
            <span className="px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">System Intelligence</span>
            <span>Vector Store</span>
          </div>
          <CardTitle className="text-2xl">Dataset & Index Manager</CardTitle>
          <CardDescription>
            Manage vector embeddings, rebuild indices, and monitor storage usage.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-4">
            <Button>
              <RefreshCw className="mr-2 h-4 w-4" /> Rebuild Active Index
            </Button>
            <Button variant="destructive">
              <Trash2 className="mr-2 h-4 w-4" /> Purge Stale Embeddings
            </Button>
          </div>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle>Index Versions</CardTitle>
          <CardDescription>History of vector store snapshots</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="rounded-md border">
            <table className="w-full text-sm text-left">
              <thead className="bg-muted/50 text-muted-foreground font-medium">
                <tr>
                  <th className="h-10 px-4 align-middle">Version</th>
                  <th className="h-10 px-4 align-middle">Date Created</th>
                  <th className="h-10 px-4 align-middle">Size</th>
                  <th className="h-10 px-4 align-middle">Status</th>
                  <th className="h-10 px-4 align-middle text-right">Actions</th>
                </tr>
              </thead>
              <tbody className="divide-y">
                {INDEX_VERSIONS.map((version, i) => (
                  <tr key={`${version.version}-${i}`} className="hover:bg-muted/50">
                    <td className="p-4 font-mono font-medium">
                      <div className="flex items-center gap-2">
                        <Database className="h-3 w-3 text-muted-foreground" />
                        {version.version}
                      </div>
                    </td>
                    <td className="p-4 text-muted-foreground">
                      <div className="flex items-center gap-2">
                        <Clock className="h-3 w-3" /> {version.date}
                      </div>
                    </td>
                    <td className="p-4 font-mono text-xs">
                      <div className="flex items-center gap-2">
                        <HardDrive className="h-3 w-3 text-muted-foreground" />
                        {version.size}
                      </div>
                    </td>
                    <td className="p-4">
                      <div className="flex items-center gap-2">
                        {version.status === "Active" ? (
                          <CheckCircle2 className="h-4 w-4 text-green-500" />
                        ) : (
                          <div className="h-2 w-2 rounded-full bg-muted-foreground/30" />
                        )}
                        <span className={version.status === "Active" ? "font-medium text-foreground" : "text-muted-foreground"}>
                          {version.status}
                        </span>
                      </div>
                    </td>
                    <td className="p-4 text-right">
                      <Button variant="ghost" size="sm" disabled={version.status === "Active"}>
                        Restore
                      </Button>
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
