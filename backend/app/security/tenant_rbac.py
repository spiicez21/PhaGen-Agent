"""
Tenant Isolation & Role-Based Access Control (RBAC)

Provides multi-tenant data isolation and fine-grained permission management
for pharmaceutical organizations using PhaGen.
"""
import logging
from typing import Optional, List, Set, Dict, Any
from enum import Enum
from functools import wraps
from datetime import datetime, timezone

from fastapi import HTTPException, Depends, Request
from sqlalchemy import Column, String, DateTime, Boolean, Integer, ForeignKey, Index
from sqlalchemy.orm import Session, relationship, declarative_base

logger = logging.getLogger(__name__)

Base = declarative_base()


class Role(str, Enum):
    """User roles with hierarchical permissions"""
    SUPER_ADMIN = "super_admin"  # Cross-tenant admin
    TENANT_ADMIN = "tenant_admin"  # Tenant administrator
    ANALYST = "analyst"  # Can create/view jobs, reports
    VIEWER = "viewer"  # Read-only access
    GUEST = "guest"  # Limited preview access


class Permission(str, Enum):
    """Granular permissions for operations"""
    # Job permissions
    JOB_CREATE = "job:create"
    JOB_VIEW = "job:view"
    JOB_VIEW_ALL = "job:view:all"  # View other users' jobs
    JOB_DELETE = "job:delete"
    
    # Report permissions
    REPORT_GENERATE = "report:generate"
    REPORT_VIEW = "report:view"
    REPORT_DOWNLOAD = "report:download"
    REPORT_SHARE = "report:share"
    
    # Data permissions
    DATA_UPLOAD = "data:upload"
    DATA_VIEW = "data:view"
    DATA_DELETE = "data:delete"
    
    # Admin permissions
    USER_MANAGE = "user:manage"
    TENANT_MANAGE = "tenant:manage"
    SETTINGS_MANAGE = "settings:manage"
    
    # Crawler permissions
    CRAWLER_RUN = "crawler:run"
    CRAWLER_VIEW = "crawler:view"


# Role-Permission Mapping
ROLE_PERMISSIONS: Dict[Role, Set[Permission]] = {
    Role.SUPER_ADMIN: set(Permission),  # All permissions
    
    Role.TENANT_ADMIN: {
        Permission.JOB_CREATE,
        Permission.JOB_VIEW,
        Permission.JOB_VIEW_ALL,
        Permission.JOB_DELETE,
        Permission.REPORT_GENERATE,
        Permission.REPORT_VIEW,
        Permission.REPORT_DOWNLOAD,
        Permission.REPORT_SHARE,
        Permission.DATA_UPLOAD,
        Permission.DATA_VIEW,
        Permission.DATA_DELETE,
        Permission.USER_MANAGE,
        Permission.SETTINGS_MANAGE,
        Permission.CRAWLER_VIEW,
    },
    
    Role.ANALYST: {
        Permission.JOB_CREATE,
        Permission.JOB_VIEW,
        Permission.REPORT_GENERATE,
        Permission.REPORT_VIEW,
        Permission.REPORT_DOWNLOAD,
        Permission.REPORT_SHARE,
        Permission.DATA_UPLOAD,
        Permission.DATA_VIEW,
    },
    
    Role.VIEWER: {
        Permission.JOB_VIEW,
        Permission.REPORT_VIEW,
        Permission.REPORT_DOWNLOAD,
        Permission.DATA_VIEW,
    },
    
    Role.GUEST: {
        Permission.JOB_VIEW,
        Permission.REPORT_VIEW,
    },
}


class Tenant(Base):
    """Tenant/Organization model for multi-tenancy"""
    __tablename__ = "tenants"
    
    tenant_id = Column(String(36), primary_key=True)
    name = Column(String(255), nullable=False)
    domain = Column(String(255), unique=True)  # e.g., "acme-pharma"
    
    # Settings
    is_active = Column(Boolean, default=True)
    max_users = Column(Integer, default=50)
    max_jobs_per_month = Column(Integer, default=1000)
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime(timezone=True), onupdate=lambda: datetime.now(timezone.utc))
    
    # Relationships
    users = relationship("TenantUser", back_populates="tenant")
    
    __table_args__ = (
        Index("idx_tenant_domain", "domain"),
        Index("idx_tenant_active", "is_active"),
    )


class TenantUser(Base):
    """User model with tenant association"""
    __tablename__ = "tenant_users"
    
    user_id = Column(String(36), primary_key=True)
    tenant_id = Column(String(36), ForeignKey("tenants.tenant_id"), nullable=False)
    email = Column(String(255), nullable=False, unique=True)
    name = Column(String(255))
    
    # Role & permissions
    role = Column(String(50), nullable=False, default=Role.VIEWER.value)
    
    # Status
    is_active = Column(Boolean, default=True)
    last_login_at = Column(DateTime(timezone=True))
    
    # Timestamps
    created_at = Column(DateTime(timezone=True), default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime(timezone=True), onupdate=lambda: datetime.now(timezone.utc))
    
    # Relationships
    tenant = relationship("Tenant", back_populates="users")
    
    __table_args__ = (
        Index("idx_user_tenant", "tenant_id"),
        Index("idx_user_email", "email"),
        Index("idx_user_active", "is_active"),
    )


class TenantContext:
    """
    Thread-local tenant context for request-scoped operations
    
    Usage:
        with TenantContext(tenant_id="tenant-123", user_id="user-456"):
            # All DB operations are tenant-scoped
            jobs = get_jobs()
    """
    
    def __init__(self, tenant_id: str, user_id: str, role: Role):
        self.tenant_id = tenant_id
        self.user_id = user_id
        self.role = role
        self.permissions = ROLE_PERMISSIONS.get(role, set())
    
    def has_permission(self, permission: Permission) -> bool:
        """Check if current user has permission"""
        return permission in self.permissions
    
    def __enter__(self):
        # Store in request context (implementation depends on framework)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # Clean up context
        pass


# Global context (should use contextvars in production)
_current_context: Optional[TenantContext] = None


def get_current_context() -> TenantContext:
    """Get current tenant context"""
    if _current_context is None:
        raise HTTPException(status_code=401, detail="No authentication context")
    return _current_context


def set_current_context(context: TenantContext):
    """Set current tenant context (for testing/mocking)"""
    global _current_context
    _current_context = context


def require_permission(permission: Permission):
    """
    Decorator to enforce permission requirements
    
    Usage:
        @require_permission(Permission.JOB_CREATE)
        def create_job(request: JobCreateRequest):
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            context = get_current_context()
            
            if not context.has_permission(permission):
                logger.warning(
                    f"[RBAC] Permission denied: {permission.value} "
                    f"for user {context.user_id} (role: {context.role.value})"
                )
                raise HTTPException(
                    status_code=403,
                    detail=f"Permission denied: {permission.value}"
                )
            
            return func(*args, **kwargs)
        return wrapper
    return decorator


def require_role(min_role: Role):
    """
    Decorator to enforce minimum role requirement
    
    Usage:
        @require_role(Role.ANALYST)
        def advanced_analysis():
            ...
    """
    role_hierarchy = {
        Role.GUEST: 0,
        Role.VIEWER: 1,
        Role.ANALYST: 2,
        Role.TENANT_ADMIN: 3,
        Role.SUPER_ADMIN: 4,
    }
    
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            context = get_current_context()
            
            if role_hierarchy.get(context.role, 0) < role_hierarchy.get(min_role, 99):
                logger.warning(
                    f"[RBAC] Role requirement not met: requires {min_role.value}, "
                    f"user has {context.role.value}"
                )
                raise HTTPException(
                    status_code=403,
                    detail=f"Requires role: {min_role.value}"
                )
            
            return func(*args, **kwargs)
        return wrapper
    return decorator


class TenantIsolationMixin:
    """
    Mixin for SQLAlchemy models to add tenant isolation
    
    Usage:
        class Job(Base, TenantIsolationMixin):
            __tablename__ = "jobs"
            job_id = Column(String(36), primary_key=True)
            # tenant_id column added automatically
    """
    
    tenant_id = Column(String(36), nullable=False, index=True)
    created_by = Column(String(36), nullable=False)  # user_id
    
    @classmethod
    def query_for_tenant(cls, session: Session, tenant_id: str):
        """Query with automatic tenant filtering"""
        return session.query(cls).filter(cls.tenant_id == tenant_id)
    
    @classmethod
    def get_by_id(cls, session: Session, record_id: str, tenant_id: str):
        """Get single record with tenant check"""
        return session.query(cls).filter(
            cls.id == record_id,
            cls.tenant_id == tenant_id
        ).first()


def apply_tenant_filter(query, model, context: Optional[TenantContext] = None):
    """
    Apply tenant isolation filter to SQLAlchemy query
    
    Usage:
        query = session.query(Job)
        query = apply_tenant_filter(query, Job)
        jobs = query.all()
    """
    if context is None:
        context = get_current_context()
    
    # Super admins can see all tenants
    if context.role == Role.SUPER_ADMIN:
        return query
    
    # Apply tenant filter
    if hasattr(model, 'tenant_id'):
        query = query.filter(model.tenant_id == context.tenant_id)
    
    # Apply user filter if not admin
    if context.role not in {Role.SUPER_ADMIN, Role.TENANT_ADMIN}:
        if hasattr(model, 'created_by'):
            query = query.filter(model.created_by == context.user_id)
    
    return query


class TenantManager:
    """Manager for tenant operations"""
    
    def __init__(self, session: Session):
        self.session = session
    
    def create_tenant(
        self,
        tenant_id: str,
        name: str,
        domain: str,
        max_users: int = 50,
        max_jobs_per_month: int = 1000,
    ) -> Tenant:
        """Create new tenant"""
        tenant = Tenant(
            tenant_id=tenant_id,
            name=name,
            domain=domain,
            max_users=max_users,
            max_jobs_per_month=max_jobs_per_month,
            is_active=True,
        )
        self.session.add(tenant)
        self.session.commit()
        
        logger.info(f"[TENANT] Created tenant: {tenant_id} ({name})")
        return tenant
    
    def create_user(
        self,
        user_id: str,
        tenant_id: str,
        email: str,
        name: str,
        role: Role = Role.VIEWER,
    ) -> TenantUser:
        """Create new user in tenant"""
        user = TenantUser(
            user_id=user_id,
            tenant_id=tenant_id,
            email=email,
            name=name,
            role=role.value,
            is_active=True,
        )
        self.session.add(user)
        self.session.commit()
        
        logger.info(f"[TENANT] Created user: {email} in tenant {tenant_id} with role {role.value}")
        return user
    
    def update_user_role(self, user_id: str, new_role: Role) -> TenantUser:
        """Update user role"""
        user = self.session.query(TenantUser).filter(
            TenantUser.user_id == user_id
        ).first()
        
        if not user:
            raise ValueError(f"User not found: {user_id}")
        
        old_role = user.role
        user.role = new_role.value
        user.updated_at = datetime.now(timezone.utc)
        self.session.commit()
        
        logger.info(f"[TENANT] Updated user {user_id} role: {old_role} -> {new_role.value}")
        return user
    
    def get_tenant_users(self, tenant_id: str) -> List[TenantUser]:
        """Get all users in tenant"""
        return self.session.query(TenantUser).filter(
            TenantUser.tenant_id == tenant_id
        ).all()
    
    def deactivate_user(self, user_id: str):
        """Deactivate user"""
        user = self.session.query(TenantUser).filter(
            TenantUser.user_id == user_id
        ).first()
        
        if user:
            user.is_active = False
            user.updated_at = datetime.now(timezone.utc)
            self.session.commit()
            logger.info(f"[TENANT] Deactivated user: {user_id}")


# FastAPI dependency for tenant context
async def get_tenant_context(request: Request) -> TenantContext:
    """
    FastAPI dependency to extract tenant context from request
    
    Usage:
        @router.get("/jobs")
        def get_jobs(context: TenantContext = Depends(get_tenant_context)):
            ...
    """
    # Extract from JWT token or session
    # This is a placeholder - implement based on your auth system
    
    tenant_id = request.headers.get("X-Tenant-ID")
    user_id = request.headers.get("X-User-ID")
    role_str = request.headers.get("X-User-Role", "viewer")
    
    if not tenant_id or not user_id:
        raise HTTPException(
            status_code=401,
            detail="Missing tenant or user context"
        )
    
    try:
        role = Role(role_str)
    except ValueError:
        role = Role.VIEWER
    
    context = TenantContext(tenant_id=tenant_id, user_id=user_id, role=role)
    set_current_context(context)
    
    return context


# PostgreSQL Row-Level Security (RLS) policy helper
def create_rls_policies(engine):
    """
    Create PostgreSQL Row-Level Security policies for tenant isolation
    
    Example SQL:
        ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
        
        CREATE POLICY tenant_isolation ON jobs
        USING (tenant_id = current_setting('app.tenant_id')::TEXT);
    """
    rls_sql = """
    -- Enable RLS on all tenant-scoped tables
    ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
    ALTER TABLE reports ENABLE ROW LEVEL SECURITY;
    ALTER TABLE documents ENABLE ROW LEVEL SECURITY;
    
    -- Create tenant isolation policies
    CREATE POLICY tenant_isolation_jobs ON jobs
    USING (tenant_id = current_setting('app.tenant_id', TRUE)::TEXT);
    
    CREATE POLICY tenant_isolation_reports ON reports
    USING (tenant_id = current_setting('app.tenant_id', TRUE)::TEXT);
    
    CREATE POLICY tenant_isolation_documents ON documents
    USING (tenant_id = current_setting('app.tenant_id', TRUE)::TEXT);
    
    -- Super admin bypass (optional)
    CREATE POLICY super_admin_bypass ON jobs
    USING (current_setting('app.role', TRUE) = 'super_admin');
    """
    
    with engine.connect() as conn:
        conn.execute(rls_sql)
        conn.commit()
    
    logger.info("[RLS] Row-Level Security policies created")


def set_tenant_context_in_session(session: Session, tenant_id: str, role: str):
    """Set tenant context variables in PostgreSQL session for RLS"""
    session.execute(f"SET app.tenant_id = '{tenant_id}'")
    session.execute(f"SET app.role = '{role}'")
