"use client";

import Link from "next/link"
import { Button } from "@/components/ui/button"
import { useAuth } from "@/contexts/AuthContext"
import { LogOut, User as UserIcon, Beaker } from "lucide-react"
import { ThemeToggle } from "./ThemeToggle"

export function Header() {
  const { user, isAuthenticated, logout } = useAuth();

  return (
    <header className="sticky top-0 z-50 w-full border-b bg-background/80 backdrop-blur-md supports-backdrop-filter:bg-background/60">
      <div className="container flex h-16 max-w-screen-2xl items-center px-6">
        <div className="mr-8 flex">
          <Link href="/" className="flex items-center space-x-2 group">
            <div className="rounded-lg bg-primary p-1.5">
              <Beaker className="h-5 w-5 text-primary-foreground" />
            </div>
            <span className="font-semibold text-lg tracking-tight">
              PhaGen
            </span>
          </Link>
        </div>
        
        <nav className="hidden md:flex items-center gap-1 text-sm font-medium flex-1">
          <Link href="/molecule" className="px-3 py-2 rounded-md transition-colors hover:bg-muted text-muted-foreground hover:text-foreground">
            Molecule
          </Link>
          <Link href="/job" className="px-3 py-2 rounded-md transition-colors hover:bg-muted text-muted-foreground hover:text-foreground">
            Job Status
          </Link>
          <Link href="/results" className="px-3 py-2 rounded-md transition-colors hover:bg-muted text-muted-foreground hover:text-foreground">
            Results
          </Link>
          <Link href="/history" className="px-3 py-2 rounded-md transition-colors hover:bg-muted text-muted-foreground hover:text-foreground">
            History
          </Link>
        </nav>
        
        <div className="flex items-center gap-3">
          <ThemeToggle />
          
          {isAuthenticated ? (
            <>
              <div className="hidden sm:flex items-center gap-2 px-3 py-1.5 rounded-md bg-muted text-sm">
                <UserIcon className="h-4 w-4 text-muted-foreground" />
                <span className="font-medium">{user?.name || user?.email}</span>
              </div>
              <Link href="/admin/crawler">
                <Button variant="ghost" size="sm" className="hidden lg:flex">
                  Admin
                </Button>
              </Link>
              <Button variant="ghost" size="sm" onClick={logout} className="gap-2">
                <LogOut className="h-4 w-4" />
                <span className="hidden sm:inline">Logout</span>
              </Button>
            </>
          ) : (
            <>
              <Link href="/login">
                <Button variant="ghost" size="sm">
                  Login
                </Button>
              </Link>
              <Link href="/signup">
                <Button size="sm" className="font-medium">
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
