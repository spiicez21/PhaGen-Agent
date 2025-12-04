import { load } from "cheerio";

export type NormalizedChunk = {
  text: string;
  charCount: number;
  wordCount: number;
  chunkIndex: number;
};

export type NormalizedDocument = {
  chunks: NormalizedChunk[];
  redactions: {
    emails: number;
    phones: number;
  };
  originalLength: number;
};

const DEFAULT_CHUNK_SIZE = 1200;
const MIN_CHUNK_SIZE = 400;

const EMAIL_REGEX = /[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}/gi;
const PHONE_REGEX = /(?:(?:\+\d{1,3}[ \-.()]*)?(?:\d[ \-.()]*){7,}\d)/g;

function stripBoilerplate(html: string): string {
  const $ = load(html);
  $("script, style, noscript, iframe, svg, form").remove();
  $("header, footer, nav, aside").remove();

  const blocks: string[] = [];
  const selectors = "article, main, section, p, li, h1, h2, h3, h4";
  $(selectors).each((_, el) => {
    const text = $(el)
      .text()
      .replace(/\s+/g, " ")
      .trim();
    if (text.length >= 40) {
      blocks.push(text);
    }
  });

  const joined = blocks.join("\n\n");
  return joined || $("body").text().replace(/\s+/g, " ").trim();
}

function redactPII(text: string) {
  let emailRedactions = 0;
  let phoneRedactions = 0;

  const redactedEmails = text.replace(EMAIL_REGEX, () => {
    emailRedactions += 1;
    return "[redacted-email]";
  });

  const redactedPhones = redactedEmails.replace(PHONE_REGEX, (match) => {
    const normalized = match.replace(/\D+/g, "");
    if (normalized.length < 8 || normalized.length > 15) {
      return match;
    }
    phoneRedactions += 1;
    return "[redacted-phone]";
  });

  return {
    text: redactedPhones,
    emailRedactions,
    phoneRedactions,
  };
}

function chunkText(text: string, chunkSize = DEFAULT_CHUNK_SIZE): string[] {
  const paragraphs = text.split(/\n{2,}/);
  const chunks: string[] = [];
  let buffer = "";

  const flushBuffer = () => {
    const trimmed = buffer.trim();
    if (!trimmed) {
      buffer = "";
      return;
    }
    if (trimmed.length < MIN_CHUNK_SIZE && chunks.length) {
      const prev = chunks.pop() as string;
      chunks.push(`${prev}\n\n${trimmed}`.trim());
    } else {
      chunks.push(trimmed);
    }
    buffer = "";
  };

  const appendParagraph = (paragraph: string) => {
    if (!buffer) {
      buffer = paragraph;
    } else {
      buffer = `${buffer}\n\n${paragraph}`;
    }
  };

  for (const paragraph of paragraphs) {
    if (!paragraph.trim()) {
      continue;
    }
    const candidate = buffer ? `${buffer}\n\n${paragraph}` : paragraph;
    if (candidate.length > chunkSize) {
      if (buffer) {
        flushBuffer();
      }
      if (paragraph.length > chunkSize) {
        let start = 0;
        while (start < paragraph.length) {
          const slice = paragraph.slice(start, start + chunkSize).trim();
          if (slice) {
            chunks.push(slice);
          }
          start += chunkSize;
        }
        buffer = "";
        continue;
      }
      buffer = paragraph;
    } else {
      appendParagraph(paragraph);
    }
  }

  if (buffer) {
    flushBuffer();
  }

  return chunks.filter((chunk) => chunk.length >= 60);
}

export function normalizeHtmlDocument(html: string): NormalizedDocument {
  const stripped = stripBoilerplate(html);
  const { text, emailRedactions, phoneRedactions } = redactPII(stripped);
  const chunkTexts = chunkText(text);

  const chunks: NormalizedChunk[] = chunkTexts.map((chunk, index) => ({
    text: chunk,
    charCount: chunk.length,
    wordCount: chunk.split(/\s+/).length,
    chunkIndex: index,
  }));

  return {
    chunks,
    originalLength: text.length,
    redactions: {
      emails: emailRedactions,
      phones: phoneRedactions,
    },
  };
}
