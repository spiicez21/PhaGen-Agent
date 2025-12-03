import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "PhaGen Agentic",
  description: "Agentic molecule repurposing assistant"
};

export default function RootLayout({
  children
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="bg-slate-950 text-slate-100">
        <main className="min-h-screen max-w-4xl mx-auto p-6">{children}</main>
      </body>
    </html>
  );
}
