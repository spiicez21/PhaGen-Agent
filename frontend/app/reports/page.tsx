"use client";

import { useState } from "react";
import { REPORT_SECTIONS, SAMPLE_PAYLOAD } from "../sample-data";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { FileText, Download, Code, CheckCircle2, AlertCircle, Loader2, FileJson } from "lucide-react";

const downloadBlob = (blob: Blob, filename: string) => {
  const url = window.URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  window.URL.revokeObjectURL(url);
};

export default function ReportsPage() {
  const [jobId, setJobId] = useState<string>("");
  const [isDownloading, setIsDownloading] = useState(false);
  const [isExportingJson, setIsExportingJson] = useState(false);
  const [status, setStatus] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const normalizedJobId = jobId.trim();

  const requireJobId = (): string | null => {
    if (!normalizedJobId) {
      setError("Enter a job ID before exporting.");
      return null;
    }
    setError(null);
    return normalizedJobId;
  };

  const handlePdfDownload = async () => {
    const id = requireJobId();
    if (!id) return;
    try {
      setIsDownloading(true);
      setStatus("Preparing PDF...");
      const response = await fetch(`/api/jobs/${id}/report.pdf`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const pdfBlob = await response.blob();
      downloadBlob(pdfBlob, `phagen-report-${id}.pdf`);
      setStatus("PDF downloaded successfully.");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to download PDF.");
      setStatus(null);
    } finally {
      setIsDownloading(false);
    }
  };

  const handleJsonExport = async () => {
    const id = requireJobId();
    if (!id) return;
    try {
      setIsExportingJson(true);
      setStatus("Fetching run payload...");
      const response = await fetch(`/api/jobs/${id}`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const payload = await response.json();
      const blob = new Blob([JSON.stringify(payload, null, 2)], {
        type: "application/json",
      });
      downloadBlob(blob, `phagen-run-${id}.json`);
      setStatus("JSON exported successfully.");
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to export JSON.");
      setStatus(null);
    } finally {
      setIsExportingJson(false);
    }
  };

  const validation = SAMPLE_PAYLOAD.validation;

  return (
    <div className="space-y-8 py-8">
      <Card>
        <CardHeader className="space-y-3">
          <div className="h-1 w-8 bg-primary rounded-full" />
          <CardTitle className="text-2xl">Report Workspace</CardTitle>
          <CardDescription className="leading-relaxed">
            Review, edit, and export the structured report. PDF + JSON export hooks connect directly to the backend jobs API.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-6">
          <div className="flex flex-col md:flex-row gap-4 items-end">
            <div className="grid w-full max-w-sm items-center gap-1.5">
              <Label htmlFor="jobId">Job ID</Label>
              <Input
                id="jobId"
                placeholder="JOB-XXXX"
                value={jobId}
                onChange={(event) => setJobId(event.target.value)}
              />
            </div>
            <div className="flex gap-2">
              <Button onClick={handlePdfDownload} disabled={isDownloading}>
                {isDownloading ? <Loader2 className="mr-2 h-4 w-4 animate-spin" /> : <FileText className="mr-2 h-4 w-4" />}
                Download PDF
              </Button>
              <Button variant="outline" onClick={handleJsonExport} disabled={isExportingJson}>
                {isExportingJson ? <Loader2 className="mr-2 h-4 w-4 animate-spin" /> : <FileJson className="mr-2 h-4 w-4" />}
                Export JSON
              </Button>
            </div>
          </div>

          {(status || error) && (
            <div className={`flex items-center gap-2 text-sm p-3 rounded-md ${status ? 'bg-emerald-50 text-emerald-600' : 'bg-destructive/10 text-destructive'}`}>
              {status ? <CheckCircle2 className="h-4 w-4" /> : <AlertCircle className="h-4 w-4" />}
              {status || error}
            </div>
          )}
        </CardContent>
      </Card>

      <div className="grid gap-8 md:grid-cols-[300px_1fr]">
        <div className="space-y-6">
          <Card className="sticky top-6">
            <CardHeader>
              <CardTitle className="text-lg">Table of Contents</CardTitle>
            </CardHeader>
            <CardContent>
              <nav className="flex flex-col space-y-1">
                {REPORT_SECTIONS.map((section, i) => (
                  <a 
                    key={section.title} 
                    href={`#section-${i}`}
                    className="text-sm text-muted-foreground hover:text-primary hover:underline py-1"
                  >
                    {i + 1}. {section.title}
                  </a>
                ))}
              </nav>
            </CardContent>
          </Card>
        </div>

        <div className="space-y-8">
          {validation && (
            <Card>
              <CardHeader>
                <div className="flex items-center justify-between">
                  <CardTitle className="text-lg">Validation Snapshot</CardTitle>
                  <Badge variant={validation.status === "pass" ? "default" : "secondary"}>
                    {validation.claims_linked}/{validation.claims_total} Linked
                  </Badge>
                </div>
                <CardDescription>Claims ↔ Evidence Traceability</CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                {validation.claim_links.map((claim) => (
                  <div key={claim.claim_id} className="p-4 rounded-lg bg-muted/50 space-y-2">
                    <div className="flex items-center justify-between">
                      <Badge variant="outline" className="capitalize">{claim.worker}</Badge>
                      <span className="text-xs text-muted-foreground font-mono">ID: {claim.claim_id}</span>
                    </div>
                    <p className="text-sm">{claim.claim_text}</p>
                    <div className="text-xs text-muted-foreground">
                      Evidence IDs: {claim.evidence_ids.length ? claim.evidence_ids.join(", ") : "None"}
                    </div>
                  </div>
                ))}
              </CardContent>
            </Card>
          )}

          {REPORT_SECTIONS.map((section, i) => (
            <Card key={section.title} id={`section-${i}`}>
              <CardHeader>
                <CardTitle>{section.title}</CardTitle>
              </CardHeader>
              <CardContent>
                <div className="prose prose-sm max-w-none text-muted-foreground">
                  <p>{section.body}</p>
                </div>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>
    </div>
  );
}
