"""
Redis caching layer for retrieval results, LLM responses, and worker outputs.
"""
from __future__ import annotations

import json
import hashlib
from typing import Any, Optional
from functools import wraps
import logging

try:
    import redis
    REDIS_AVAILABLE = True
except ImportError:
    REDIS_AVAILABLE = False

logger = logging.getLogger(__name__)


class CacheClient:
    """Redis-backed cache with graceful fallback when Redis unavailable."""
    
    def __init__(self, redis_url: str = "redis://localhost:6379/0", enabled: bool = True):
        self.enabled = enabled and REDIS_AVAILABLE
        self.client: Optional[redis.Redis] = None
        
        if self.enabled:
            try:
                self.client = redis.from_url(redis_url, decode_responses=True)
                self.client.ping()
                logger.info(f"Redis cache connected: {redis_url}")
            except (redis.ConnectionError, redis.TimeoutError) as exc:
                logger.warning(f"Redis unavailable, caching disabled: {exc}")
                self.enabled = False
        else:
            if not REDIS_AVAILABLE:
                logger.info("Redis library not installed, caching disabled")
            else:
                logger.info("Caching explicitly disabled via config")
    
    def _make_key(self, prefix: str, *args: Any) -> str:
        """Generate deterministic cache key from prefix and arguments."""
        key_parts = [str(arg) for arg in args]
        key_str = ":".join([prefix] + key_parts)
        if len(key_str) > 200:
            # Hash long keys
            hash_suffix = hashlib.sha256(key_str.encode()).hexdigest()[:16]
            key_str = f"{prefix}:{hash_suffix}"
        return key_str
    
    def get(self, key: str) -> Optional[str]:
        """Get cached value by key."""
        if not self.enabled or not self.client:
            return None
        try:
            return self.client.get(key)
        except (redis.ConnectionError, redis.TimeoutError) as exc:
            logger.warning(f"Cache GET failed for {key}: {exc}")
            return None
    
    def set(self, key: str, value: str, ttl: int = 3600) -> bool:
        """Set cache value with TTL in seconds (default 1 hour)."""
        if not self.enabled or not self.client:
            return False
        try:
            self.client.setex(key, ttl, value)
            return True
        except (redis.ConnectionError, redis.TimeoutError) as exc:
            logger.warning(f"Cache SET failed for {key}: {exc}")
            return False
    
    def delete(self, key: str) -> bool:
        """Delete cached value."""
        if not self.enabled or not self.client:
            return False
        try:
            self.client.delete(key)
            return True
        except (redis.ConnectionError, redis.TimeoutError) as exc:
            logger.warning(f"Cache DELETE failed for {key}: {exc}")
            return False
    
    def get_json(self, key: str) -> Optional[dict]:
        """Get JSON-encoded cached value."""
        raw = self.get(key)
        if raw:
            try:
                return json.loads(raw)
            except json.JSONDecodeError as exc:
                logger.warning(f"Cache JSON decode failed for {key}: {exc}")
        return None
    
    def set_json(self, key: str, value: dict, ttl: int = 3600) -> bool:
        """Set JSON-encoded cache value."""
        try:
            serialized = json.dumps(value)
            return self.set(key, serialized, ttl)
        except (TypeError, ValueError) as exc:
            logger.warning(f"Cache JSON encode failed for {key}: {exc}")
            return False
    
    def invalidate_pattern(self, pattern: str) -> int:
        """Delete all keys matching pattern (e.g., 'retrieval:*')."""
        if not self.enabled or not self.client:
            return 0
        try:
            keys = self.client.keys(pattern)
            if keys:
                return self.client.delete(*keys)
            return 0
        except (redis.ConnectionError, redis.TimeoutError) as exc:
            logger.warning(f"Cache pattern invalidation failed for {pattern}: {exc}")
            return 0


def cached_retrieval(cache_client: CacheClient, ttl: int = 7200):
    """
    Decorator for caching retrieval results (2 hour default TTL).
    
    Usage:
        @cached_retrieval(cache_client, ttl=3600)
        def retrieve_passages(molecule: str, query: str) -> list:
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Build cache key from function name and arguments
            key_parts = [func.__name__] + [str(arg) for arg in args] + [f"{k}={v}" for k, v in sorted(kwargs.items())]
            cache_key = cache_client._make_key("retrieval", *key_parts)
            
            # Try cache first
            cached = cache_client.get_json(cache_key)
            if cached is not None:
                logger.debug(f"Cache HIT: {cache_key}")
                return cached
            
            # Cache miss - call function
            logger.debug(f"Cache MISS: {cache_key}")
            result = func(*args, **kwargs)
            
            # Store in cache
            if result is not None:
                cache_client.set_json(cache_key, result, ttl=ttl)
            
            return result
        return wrapper
    return decorator


def cached_llm_response(cache_client: CacheClient, ttl: int = 86400):
    """
    Decorator for caching LLM responses (24 hour default TTL).
    
    Usage:
        @cached_llm_response(cache_client, ttl=86400)
        def generate_summary(prompt: str, model: str) -> str:
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Build cache key from function name and arguments
            key_parts = [func.__name__] + [str(arg)[:500] for arg in args] + [f"{k}={str(v)[:500]}" for k, v in sorted(kwargs.items())]
            cache_key = cache_client._make_key("llm", *key_parts)
            
            # Try cache first
            cached = cache_client.get(cache_key)
            if cached is not None:
                logger.debug(f"Cache HIT: {cache_key}")
                return cached
            
            # Cache miss - call function
            logger.debug(f"Cache MISS: {cache_key}")
            result = func(*args, **kwargs)
            
            # Store in cache
            if result:
                cache_client.set(cache_key, result, ttl=ttl)
            
            return result
        return wrapper
    return decorator


# Global cache instance (initialized in main.py)
_cache_client: Optional[CacheClient] = None


def get_cache() -> CacheClient:
    """Get global cache client instance."""
    global _cache_client
    if _cache_client is None:
        from .config import get_settings
        settings = get_settings()
        redis_url = getattr(settings, "redis_url", "redis://localhost:6379/0")
        cache_enabled = getattr(settings, "cache_enabled", True)
        _cache_client = CacheClient(redis_url=redis_url, enabled=cache_enabled)
    return _cache_client
