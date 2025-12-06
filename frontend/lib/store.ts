/**
 * Client-side state management for job history and authentication
 * Uses localStorage for persistence
 */

export interface StoredJob {
  job_id: string;
  smiles: string;
  molecule?: string;
  status: string;
  created_at: string;
  updated_at?: string;
  recommendation?: string;
}

const JOBS_KEY = 'phagen_job_history';
const AUTH_KEY = 'phagen_auth_token';
const USER_KEY = 'phagen_user';

/**
 * Job History Management
 */
export const jobStore = {
  // Get all stored jobs
  getAll(): StoredJob[] {
    if (typeof window === 'undefined') return [];
    try {
      const stored = localStorage.getItem(JOBS_KEY);
      return stored ? JSON.parse(stored) : [];
    } catch {
      return [];
    }
  },

  // Add a new job to history
  add(job: StoredJob): void {
    if (typeof window === 'undefined') return;
    try {
      const jobs = this.getAll();
      // Avoid duplicates
      if (!jobs.find(j => j.job_id === job.job_id)) {
        jobs.unshift(job);
        // Keep only last 50 jobs
        if (jobs.length > 50) jobs.pop();
        localStorage.setItem(JOBS_KEY, JSON.stringify(jobs));
      }
    } catch (error) {
      console.error('Failed to save job:', error);
    }
  },

  // Update existing job
  update(job_id: string, updates: Partial<StoredJob>): void {
    if (typeof window === 'undefined') return;
    try {
      const jobs = this.getAll();
      const index = jobs.findIndex(j => j.job_id === job_id);
      if (index !== -1) {
        jobs[index] = { ...jobs[index], ...updates, updated_at: new Date().toISOString() };
        localStorage.setItem(JOBS_KEY, JSON.stringify(jobs));
      }
    } catch (error) {
      console.error('Failed to update job:', error);
    }
  },

  // Get specific job
  get(job_id: string): StoredJob | null {
    return this.getAll().find(j => j.job_id === job_id) || null;
  },

  // Clear all jobs
  clear(): void {
    if (typeof window === 'undefined') return;
    localStorage.removeItem(JOBS_KEY);
  }
};

/**
 * Authentication Management
 */
export interface User {
  id: string;
  email: string;
  name: string;
  role?: string;
}

export const authStore = {
  // Get auth token
  getToken(): string | null {
    if (typeof window === 'undefined') return null;
    return localStorage.getItem(AUTH_KEY);
  },

  // Set auth token
  setToken(token: string): void {
    if (typeof window === 'undefined') return;
    localStorage.setItem(AUTH_KEY, token);
  },

  // Remove auth token
  clearToken(): void {
    if (typeof window === 'undefined') return;
    localStorage.removeItem(AUTH_KEY);
    localStorage.removeItem(USER_KEY);
  },

  // Get user info
  getUser(): User | null {
    if (typeof window === 'undefined') return null;
    try {
      const stored = localStorage.getItem(USER_KEY);
      return stored ? JSON.parse(stored) : null;
    } catch {
      return null;
    }
  },

  // Set user info
  setUser(user: User): void {
    if (typeof window === 'undefined') return;
    localStorage.setItem(USER_KEY, JSON.stringify(user));
  },

  // Check if authenticated
  isAuthenticated(): boolean {
    return !!this.getToken();
  }
};
