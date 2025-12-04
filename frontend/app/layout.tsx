import type { Metadata } from "next";
import Link from "next/link";
import { Poppins } from "next/font/google";
import "./globals.css";

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
    <html lang="en">
      <body className={`${poppins.variable} antialiased page-shell`}>
        <div className="app-shell">
          <header className="global-header">
            <div>
              <p className="logo-mark">PhaGen</p>
              <p className="logo-subhead">Agentic Intelligence</p>
            </div>
            <nav className="nav-links">
              <Link href="/">Home</Link>
              <Link href="/molecule">Molecule Search</Link>
              <Link href="/job">Job Status</Link>
              <Link href="/results">Results</Link>
              <Link href="/comparison">Comparison</Link>
              <Link href="/evidence/clinical">Evidence</Link>
              <Link href="/reports">Reports</Link>
              <Link href="/history">History</Link>
              <Link href="/admin/crawler">Admin</Link>
            </nav>
          </header>
          <main className="app-main">{children}</main>
        </div>
      </body>
    </html>
  );
}
