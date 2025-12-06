"""
Health check endpoints for load balancer integration.
"""
from fastapi import APIRouter, Response, status
from ..database import SessionLocal
from ..cache import get_cache
import logging

logger = logging.getLogger(__name__)

router = APIRouter(tags=["health"])


@router.get("/health")
async def health_check():
    """
    Basic health check - returns 200 if service is up.
    Used by load balancers for routing decisions.
    """
    return {"status": "healthy", "service": "phagen-api"}


@router.get("/ready")
async def readiness_check(response: Response):
    """
    Readiness check - validates critical dependencies.
    Returns 200 only when service can handle requests.
    Used by Kubernetes/load balancers to determine if pod should receive traffic.
    """
    checks = {
        "database": False,
        "cache": False,
    }
    
    # Check database connectivity
    try:
        db = SessionLocal()
        db.execute("SELECT 1")
        db.close()
        checks["database"] = True
    except Exception as exc:
        logger.warning(f"Database readiness check failed: {exc}")
    
    # Check cache connectivity (non-critical)
    try:
        cache = get_cache()
        if cache.enabled:
            cache.set("readiness:test", "ok", ttl=10)
            checks["cache"] = cache.get("readiness:test") == "ok"
        else:
            checks["cache"] = True  # Cache is optional
    except Exception as exc:
        logger.warning(f"Cache readiness check failed: {exc}")
        checks["cache"] = True  # Non-critical
    
    # Service is ready only if database is accessible
    all_ready = checks["database"]
    
    if not all_ready:
        response.status_code = status.HTTP_503_SERVICE_UNAVAILABLE
        return {"status": "not_ready", "checks": checks}
    
    return {"status": "ready", "checks": checks}


@router.get("/metrics")
async def metrics_endpoint():
    """
    Basic metrics endpoint for Prometheus scraping.
    Returns simple metrics about the service.
    """
    from ..jobs import job_store
    
    try:
        # Get job statistics
        total_jobs = len(job_store._store) if hasattr(job_store, '_store') else 0
        
        metrics_data = {
            "total_jobs": total_jobs,
            "cache_enabled": get_cache().enabled,
        }
        
        # Format as Prometheus-style metrics
        metrics_lines = [
            "# HELP phagen_total_jobs Total number of jobs in store",
            "# TYPE phagen_total_jobs gauge",
            f"phagen_total_jobs {total_jobs}",
            "",
            "# HELP phagen_cache_enabled Whether Redis cache is enabled",
            "# TYPE phagen_cache_enabled gauge",
            f"phagen_cache_enabled {1 if get_cache().enabled else 0}",
        ]
        
        return Response(content="\n".join(metrics_lines), media_type="text/plain")
    except Exception as exc:
        logger.error(f"Metrics endpoint error: {exc}")
        return Response(content="# Error generating metrics\n", media_type="text/plain")
