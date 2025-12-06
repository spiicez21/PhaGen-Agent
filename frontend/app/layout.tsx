import type { Metadata } from "next";
import { Poppins } from "next/font/google";
import "./globals.css";
import { Header } from "@/components/Header";
import { cn } from "@/lib/utils";

const poppins = Poppins({
  subsets: ["latin"],
  weight: ["400", "500", "600", "700"],
  variable: "--font-poppins",
});

export const metadata: Metadata = {
  title: "PhaGen Agentic",
  description: "Agentic molecule repurposing copilot UI",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body className={cn("min-h-screen bg-background font-sans antialiased", poppins.variable)}>
        <div className="relative flex min-h-screen flex-col">
          <Header />
          <main className="flex-1 container max-w-screen-2xl mx-auto px-4 py-6">{children}</main>
        </div>
      </body>
    </html>
  );
}

