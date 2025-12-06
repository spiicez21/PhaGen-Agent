import Link from "next/link"
import { Button } from "@/components/ui/button"

export function Header() {
  return (
    <header className="sticky top-0 z-50 w-full border-b border-border/40 bg-background/95 backdrop-blur supports-backdrop-filter:bg-background/60">
      <div className="container flex h-14 max-w-screen-2xl items-center px-4">
        <div className="mr-4 hidden md:flex">
          <Link href="/" className="mr-6 flex items-center space-x-2">
            <span className="hidden font-bold sm:inline-block text-lg tracking-tight">
              PhaGen Agentic
            </span>
          </Link>
          <nav className="flex items-center gap-6 text-sm font-medium">
            <Link href="/molecule" className="transition-colors hover:text-foreground/80 text-foreground/60">Molecule Search</Link>
            <Link href="/job" className="transition-colors hover:text-foreground/80 text-foreground/60">Job Status</Link>
            <Link href="/results" className="transition-colors hover:text-foreground/80 text-foreground/60">Results</Link>
            <Link href="/comparison" className="transition-colors hover:text-foreground/80 text-foreground/60">Comparison</Link>
            <Link href="/history" className="transition-colors hover:text-foreground/80 text-foreground/60">History</Link>
          </nav>
        </div>
        <div className="flex flex-1 items-center justify-between space-x-2 md:justify-end">
          <div className="w-full flex-1 md:w-auto md:flex-none">
          </div>
          <nav className="flex items-center">
             <Link href="/admin/crawler">
                <Button variant="ghost" size="sm">Admin</Button>
             </Link>
          </nav>
        </div>
      </div>
    </header>
  )
}
