import { BasicCrawler, Configuration, Dataset } from "@crawlee/basic";

const seedUrls = [
  "https://clinicaltrials.gov/study/NCT00000000",
  "https://www.ncbi.nlm.nih.gov/pmc/articles/mock",
  "https://patents.google.com/patent/US0000000A/en"
];

const config = new Configuration({ defaultDatasetId: "phagen-seed" });

const crawler = new BasicCrawler({
  requestHandler: async ({ request, log, sendRequest }) => {
    const response = await sendRequest({ useExtendedUniqueKey: true });
    const text = response.body?.toString() ?? "";
    log.info(`Fetched ${request.loadedUrl}`);
    await Dataset.pushData({
      url: request.loadedUrl,
      text: text.slice(0, 5000),
      source_type: "mock",
    });
  },
  maxRequestsPerCrawl: seedUrls.length,
  requestList: await Configuration.merge({
    sources: seedUrls.map((url) => ({ url }))
  })
});

await crawler.run();
console.log("Seed crawl completed");
