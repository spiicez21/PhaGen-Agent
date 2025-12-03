import fetch from "node-fetch";
import { readFile } from "node:fs/promises";

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

  const robots = await checkRobots(task.crawlUrl);
  if (!robots.allowed) {
    console.warn(`Skipped ${task.crawlUrl} — blocked by robots.txt`);
    continue;
  }
  if (robots.delay > 0) {
    await sleep(robots.delay * 1000);
  }

  const response = await fetch(task.crawlUrl, { headers: { "User-Agent": USER_AGENT } });
  if (!response.ok) {
    console.warn(`Failed to fetch ${task.crawlUrl}: ${response.status}`);
    continue;
  }
  const html = await response.text();
  await Dataset.pushData({
    id: task.id,
    via: "html",
    url: task.crawlUrl,
    source_type: task.source_type,
    snippet: html.slice(0, SNIPPET_LIMIT),
  });
  console.log(`Crawled ${task.crawlUrl}`);
}

console.log("Ingestion finished with API-first + robots fallback");
