import fetch from "node-fetch";
import robotsParser from "robots-parser";
import { URL } from "node:url";

const USER_AGENT = "PhaGenBot/1.0 (+mailto:admin@phagen.ai)";
const cache = new Map<string, ReturnType<typeof robotsParser>>();

export type RobotsCheckResult = {
  allowed: boolean;
  delay: number;
};

export async function checkRobots(url: string): Promise<RobotsCheckResult> {
  const parsed = new URL(url);
  const robotsUrl = `${parsed.protocol}//${parsed.host}/robots.txt`;

  let parser = cache.get(parsed.host);
  if (!parser) {
    try {
      const res = await fetch(robotsUrl, { headers: { "User-Agent": USER_AGENT } });
      if (!res.ok) {
        parser = robotsParser(robotsUrl, "");
      } else {
        const body = await res.text();
        parser = robotsParser(robotsUrl, body);
      }
      cache.set(parsed.host, parser);
    } catch (error) {
      console.warn(`robots.txt fetch error for ${parsed.host}:`, (error as Error).message);
      return { allowed: false, delay: 0 };
    }
  }

  const allowed = parser.isAllowed(url, USER_AGENT) ?? false;
  const delay = parser.getCrawlDelay(USER_AGENT) ?? 0;
  return { allowed, delay };
}
