"use client";

import React, { createContext, useContext, useEffect, useState } from 'react';
import { authStore, User } from '@/lib/store';
import { api } from '@/lib/api';

interface AuthContextType {
  user: User | null;
  isAuthenticated: boolean;
  login: (email: string, password: string) => Promise<void>;
  signup: (email: string, name: string, password: string) => Promise<void>;
  logout: () => void;
  loading: boolean;
}

const AuthContext = createContext<AuthContextType | undefined>(undefined);

export function AuthProvider({ children }: { children: React.ReactNode }) {
  const [user, setUser] = useState<User | null>(null);
  const [loading, setLoading] = useState(true);

  // Initialize auth state from localStorage and validate token
  useEffect(() => {
    const initAuth = async () => {
      const storedUser = authStore.getUser();
      const token = authStore.getToken();

      if (storedUser && token) {
        try {
          // Validate token by fetching current user
          const currentUser = await api.getCurrentUser();
          const user: User = {
            id: currentUser.id,
            email: currentUser.email,
            name: currentUser.name,
            role: 'user'
          };
          authStore.setUser(user);
          setUser(user);
        } catch {
          // Token invalid or expired, clear auth
          authStore.clearToken();
          setUser(null);
        }
      }
      setLoading(false);
    };

    initAuth();
  }, []);

  const login = async (email: string, password: string) => {
    setLoading(true);
    try {
      // Call backend login API
      const tokenResponse = await api.login({ email, password });
      
      // Store token
      authStore.setToken(tokenResponse.access_token);
      
      // Fetch user details
      const currentUser = await api.getCurrentUser();
      const user: User = {
        id: currentUser.id,
        email: currentUser.email,
        name: currentUser.name,
        role: 'user'
      };
      
      authStore.setUser(user);
      setUser(user);
    } finally {
      setLoading(false);
    }
  };

  const signup = async (email: string, name: string, password: string) => {
    setLoading(true);
    try {
      // Call backend signup API
      await api.signup({ email, name, password });
      
      // Auto-login after signup
      await login(email, password);
    } finally {
      setLoading(false);
    }
  };

  const logout = () => {
    authStore.clearToken();
    setUser(null);
  };

  return (
    <AuthContext.Provider value={{ user, isAuthenticated: !!user, login, signup, logout, loading }}>
      {children}
    </AuthContext.Provider>
  );
}

export function useAuth() {
  const context = useContext(AuthContext);
  if (context === undefined) {
    throw new Error('useAuth must be used within an AuthProvider');
  }
  return context;
}
