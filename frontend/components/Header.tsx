"use client";

import Link from "next/link"
import { usePathname } from "next/navigation"
import { Button } from "@/components/ui/button"
import { useAuth } from "@/contexts/AuthContext"
import { LogOut, UserCircle2, FlaskConical, Activity, FileText, Clock, Search } from "lucide-react"
import { ThemeToggle } from "./ThemeToggle"
import { cn } from "@/lib/utils"

const NAV_ITEMS = [
  { href: "/molecule", label: "Analyze", icon: Search },
  { href: "/job",      label: "Pipeline", icon: Activity },
  { href: "/results",  label: "Results",  icon: FileText },
  { href: "/history",  label: "History",  icon: Clock },
];

export function Header() {
  const { user, isAuthenticated, logout } = useAuth();
  const pathname = usePathname();

  return (
    <header className="sticky top-0 z-50 w-full border-b border-border bg-card">
      <div className="mx-auto flex h-10 max-w-screen-xl items-center px-4 sm:px-6 lg:px-8 gap-0">

        {/* Logo */}
        <Link href="/" className="flex items-center gap-2 pr-4 border-r border-border mr-4 flex-shrink-0">
          <div className="flex h-6 w-6 items-center justify-center rounded bg-primary">
            <FlaskConical className="h-3.5 w-3.5 text-primary-foreground" />
          </div>
          <span className="text-sm font-bold tracking-tight text-foreground">PhaGen</span>
          <span className="text-[10px] font-mono text-muted-foreground hidden sm:block">R&amp;D</span>
        </Link>

        {/* Nav */}
        <nav className="hidden md:flex items-stretch h-full flex-1">
          {NAV_ITEMS.map(({ href, label, icon: Icon }) => {
            const active = pathname === href || pathname?.startsWith(href + "/");
            return (
              <Link
                key={href}
                href={href}
                className={cn(
                  "relative flex items-center gap-1.5 px-3 text-xs font-medium transition-colors h-full border-r border-border",
                  active
                    ? "text-primary bg-primary/6"
                    : "text-muted-foreground hover:text-foreground hover:bg-muted/40"
                )}
              >
                {active && (
                  <span className="absolute bottom-0 left-0 right-0 h-[2px] bg-primary" />
                )}
                <Icon className="h-3 w-3" />
                {label}
              </Link>
            );
          })}
        </nav>

        {/* Spacer */}
        <div className="flex-1 md:flex-none" />

        {/* Right side */}
        <div className="flex items-center border-l border-border pl-3 gap-2">
          <ThemeToggle />

          {isAuthenticated ? (
            <>
              <div className="hidden sm:flex items-center gap-1.5 text-[11px] text-muted-foreground">
                <UserCircle2 className="h-3.5 w-3.5" />
                <span className="font-medium max-w-[100px] truncate">{user?.name || user?.email}</span>
              </div>
              <div className="w-px h-4 bg-border" />
              <button
                onClick={logout}
                className="flex items-center gap-1 text-[11px] text-muted-foreground hover:text-destructive transition-colors"
              >
                <LogOut className="h-3 w-3" />
                <span className="hidden sm:block">Logout</span>
              </button>
            </>
          ) : (
            <>
              <Link href="/login">
                <Button variant="ghost" size="sm" className="text-xs h-7 px-2">Log in</Button>
              </Link>
              <Link href="/signup">
                <Button size="sm" className="text-xs h-7 px-3 bg-primary text-primary-foreground">
                  Sign Up
                </Button>
              </Link>
            </>
          )}
        </div>

      </div>
    </header>
  )
}
