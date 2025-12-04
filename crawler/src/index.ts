import fetch from "node-fetch";
import { readFile } from "node:fs/promises";

import { normalizeHtmlDocument } from "./normalize.ts";
import { checkRobots } from "./robots.ts";

const { Dataset } = await import("@crawlee/basic");

type SourceTask = {
  id: string;
  crawlUrl: string;
  source_type: string;
  apiMock?: string;
};

const SOURCES: SourceTask[] = [
  {
    id: "ctgov-nct01907900",
    crawlUrl: "https://clinicaltrials.gov/study/NCT01907900",
    source_type: "clinical",
    apiMock: "../mock-data/ctgov.json",
  },
  {
    id: "pmc-pirfenidone",
    crawlUrl: "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4418487/",
    source_type: "literature",
    apiMock: "../mock-data/pmc.json",
  },
  {
    id: "fda-esbriet",
    crawlUrl:
      "https://www.fda.gov/drugs/postmarket-drug-safety-information-patients-and-providers/pirfenidone-esbriet-information",
    source_type: "regulatory",
    apiMock: "../mock-data/fda.json",
  },
  {
    id: "patent-wo2013185843",
    crawlUrl: "https://patents.google.com/patent/WO2013185843A1/en",
    source_type: "patent",
  },
];

const USER_AGENT = "PhaGenBot/1.0 (+mailto:admin@phagen.ai)";
const SNIPPET_LIMIT = 5000;
const DEFAULT_DOMAIN_BUDGET = 2; // max HTML pages per host when no override is provided

const DOMAIN_BUDGET_OVERRIDES: Record<string, number> = {
  "clinicaltrials.gov": 5,
  "www.fda.gov": 3,
  "www.ncbi.nlm.nih.gov": 4,
  "patents.google.com": 2,
};

type BudgetStatus = {
  allowed: boolean;
  used: number;
  budget: number;
  remaining: number;
};

const domainUsage = new Map<string, number>();

function resolveBudget(host: string): number {
  const override = DOMAIN_BUDGET_OVERRIDES[host];
  if (typeof override === "number") {
    return override;
  }
  return DEFAULT_DOMAIN_BUDGET;
}

function checkDomainBudget(host: string): BudgetStatus {
  const budget = resolveBudget(host);
  if (budget <= 0) {
    return { allowed: true, used: 0, budget, remaining: Number.POSITIVE_INFINITY };
  }
  const used = domainUsage.get(host) ?? 0;
  const remaining = Math.max(budget - used, 0);
  return { allowed: used < budget, used, budget, remaining };
}

function incrementDomainUsage(host: string): void {
  const used = domainUsage.get(host) ?? 0;
  domainUsage.set(host, used + 1);
}

async function fetchApiData(task: SourceTask): Promise<Record<string, unknown> | null> {
  if (!task.apiMock) {
    return null;
  }
  try {
    const filePath = new URL(task.apiMock, import.meta.url);
    const raw = await readFile(filePath, "utf-8");
    return JSON.parse(raw);
  } catch (error) {
    console.warn(`API mock missing for ${task.id}:`, (error as Error).message);
    return null;
  }
}

async function sleep(ms: number) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

for (const task of SOURCES) {
  const apiData = await fetchApiData(task);
  if (apiData) {
    await Dataset.pushData({
      id: task.id,
      via: "api",
      payload: apiData,
      source_type: task.source_type,
    });
    console.log(`Stored API result for ${task.id}`);
    continue;
  }

  const host = new URL(task.crawlUrl).host;
  const budgetStatus = checkDomainBudget(host);
  if (!budgetStatus.allowed) {
    console.warn(
      `Skipped ${task.crawlUrl} — domain budget exhausted (${budgetStatus.used}/${budgetStatus.budget} pages)`,
    );
    await Dataset.pushData({
      id: `${task.id}#budget-cap`,
      via: "html",
      url: task.crawlUrl,
      source_type: task.source_type,
      normalization: {
        chunked: false,
        reason: "domain-budget-exhausted",
      },
      crawl_budget: {
        host,
        used: budgetStatus.used,
        budget: budgetStatus.budget,
        remaining: budgetStatus.remaining,
      },
    });
    continue;
  }

  const robots = await checkRobots(task.crawlUrl);
  if (!robots.allowed) {
    console.warn(`Skipped ${task.crawlUrl} — blocked by robots.txt`);
    continue;
  }
  if (robots.delay > 0) {
    await sleep(robots.delay * 1000);
  }

  incrementDomainUsage(host);
  const response = await fetch(task.crawlUrl, { headers: { "User-Agent": USER_AGENT } });
  if (!response.ok) {
    console.warn(`Failed to fetch ${task.crawlUrl}: ${response.status}`);
    continue;
  }
  const html = await response.text();
  const normalized = normalizeHtmlDocument(html);

  if (!normalized.chunks.length) {
    await Dataset.pushData({
      id: task.id,
      via: "html",
      url: task.crawlUrl,
      source_type: task.source_type,
      snippet: html.slice(0, SNIPPET_LIMIT),
      normalization: {
        chunked: false,
        reason: "empty-content",
      },
    });
    console.log(`Crawled ${task.crawlUrl} (normalization fallback)`);
    continue;
  }

  for (const chunk of normalized.chunks) {
    await Dataset.pushData({
      id: `${task.id}#${String(chunk.chunkIndex + 1).padStart(3, "0")}`,
      via: "html",
      url: task.crawlUrl,
      source_type: task.source_type,
      chunk_index: chunk.chunkIndex,
      chunk_count: normalized.chunks.length,
      char_count: chunk.charCount,
      word_count: chunk.wordCount,
      text: chunk.text,
      redactions: normalized.redactions,
      crawl_budget: {
        host,
        used: domainUsage.get(host) ?? 0,
        budget: resolveBudget(host),
        remaining: Math.max(resolveBudget(host) - (domainUsage.get(host) ?? 0), 0),
      },
    });
  }

  console.log(
    `Crawled ${task.crawlUrl} → ${normalized.chunks.length} normalized chunks (emails redacted: ${normalized.redactions.emails}, phones redacted: ${normalized.redactions.phones})`,
  );
}

console.log("Ingestion finished with API-first + robots fallback");
if (domainUsage.size) {
  const summary = Array.from(domainUsage.entries())
    .map(([host, used]) => `${host}=${used}/${resolveBudget(host)}`)
    .join(", ");
  console.log(`Domain usage summary: ${summary}`);
}
