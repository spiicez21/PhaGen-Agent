import { LITERATURE_EVIDENCE } from "../../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ExternalLink, Filter, BookOpen, FlaskConical } from "lucide-react";

export default function LiteratureEvidencePage() {
  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader>
          <div className="flex items-center gap-2 text-sm text-muted-foreground mb-2">
            <span className="px-2 py-0.5 rounded-full bg-primary/10 text-primary text-xs font-medium">Evidence Dashboard</span>
            <span>Literature</span>
          </div>
          <CardTitle className="text-2xl">Literature Intelligence</CardTitle>
          <CardDescription>
            Mechanism-of-action insights across preclinical and clinical papers.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <div className="flex flex-wrap gap-2">
            <div className="flex items-center gap-2 mr-2 text-sm font-medium text-muted-foreground">
              <Filter className="h-4 w-4" /> Filters:
            </div>
            <Badge variant="secondary">Type: Clinical + Mechanistic</Badge>
            <Badge variant="secondary">Strength ≥ 0.6</Badge>
          </div>
        </CardContent>
      </Card>

      <div className="space-y-4">
        {LITERATURE_EVIDENCE.map((entry, i) => (
          <Card key={`${entry.title}-${i}`}>
            <CardHeader>
              <div className="flex items-start justify-between gap-4">
                <div className="space-y-1">
                  <div className="flex items-center gap-2 text-xs text-muted-foreground font-mono">
                    <BookOpen className="h-3 w-3" /> {entry.doi}
                  </div>
                  <CardTitle className="text-lg leading-tight">{entry.title}</CardTitle>
                </div>
                <Badge variant={entry.strength >= 0.8 ? "default" : "secondary"}>
                  {Math.round(entry.strength * 100)}% Strength
                </Badge>
              </div>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="flex items-center gap-2 text-sm font-medium text-primary">
                <FlaskConical className="h-4 w-4" />
                Mechanism: {entry.mechanism}
              </div>
              
              <div className="p-4 rounded-md bg-muted/50 text-sm italic text-muted-foreground border-l-2 border-primary/50">
                "{entry.snippet}"
              </div>
              
              <div className="pt-2">
                <Button variant="link" className="h-auto p-0 text-xs" asChild>
                  <a href={`https://doi.org/${entry.doi}`} target="_blank" rel="noreferrer">
                    View DOI <ExternalLink className="ml-1 h-3 w-3" />
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
