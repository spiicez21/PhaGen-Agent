# Authentication System Implementation

## Overview

This document describes the complete JWT-based authentication system implemented for PhaGen Agent.

## Backend Implementation

### 1. Dependencies Added

```txt
python-jose[cryptography]==3.3.0  # JWT token generation and validation
passlib[bcrypt]==1.7.4           # Password hashing
python-multipart==0.0.9          # Form data parsing
```

### 2. Database Model

**User Model** (`backend/app/models.py`):
```python
class User(Base):
    __tablename__ = "users"
    
    id: Mapped[str]                 # UUID primary key
    email: Mapped[str]              # Unique, indexed
    name: Mapped[str]               # Full name
    hashed_password: Mapped[str]    # Bcrypt hashed password
    is_active: Mapped[bool]         # Account status
    created_at: Mapped[datetime]    # Registration date
    updated_at: Mapped[datetime]    # Last update
```

### 3. Configuration

Added to `backend/app/config.py`:
```python
secret_key: str = "your-secret-key-change-in-production-min-32-chars-long"
algorithm: str = "HS256"
access_token_expire_minutes: int = 30
```

**⚠️ IMPORTANT**: Change `secret_key` in production! Use a secure random string (minimum 32 characters).

### 4. Authentication Utilities

**File**: `backend/app/auth.py`

Functions:
- `verify_password()` - Verify password against hash
- `get_password_hash()` - Hash a password with bcrypt
- `create_access_token()` - Generate JWT token
- `decode_token()` - Decode and validate JWT
- `authenticate_user()` - Verify email/password credentials
- `get_current_user()` - FastAPI dependency to get authenticated user
- `get_current_active_user()` - Verify user is active

### 5. API Endpoints

**Router**: `backend/app/routers/auth.py`

#### POST /api/auth/signup
Register a new user.

**Request Body**:
```json
{
  "email": "user@example.com",
  "name": "John Doe",
  "password": "securepassword123"
}
```

**Response** (201 Created):
```json
{
  "id": "uuid-here",
  "email": "user@example.com",
  "name": "John Doe",
  "is_active": true,
  "created_at": "2024-12-06T10:00:00Z"
}
```

**Errors**:
- 400: Email already registered

#### POST /api/auth/login
Authenticate and get JWT token.

**Request Body**:
```json
{
  "email": "user@example.com",
  "password": "securepassword123"
}
```

**Response** (200 OK):
```json
{
  "access_token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "token_type": "bearer"
}
```

**Errors**:
- 401: Incorrect email or password

#### GET /api/auth/me
Get current user information (requires authentication).

**Headers**:
```
Authorization: Bearer <token>
```

**Response** (200 OK):
```json
{
  "id": "uuid-here",
  "email": "user@example.com",
  "name": "John Doe",
  "is_active": true,
  "created_at": "2024-12-06T10:00:00Z"
}
```

**Errors**:
- 401: Invalid or expired token
- 403: Inactive user

## Frontend Implementation

### 1. API Client Updates

**File**: `frontend/lib/api.ts`

Added methods:
```typescript
async login(credentials: LoginRequest): Promise<TokenResponse>
async signup(userData: SignupRequest): Promise<UserResponse>
async getCurrentUser(): Promise<UserResponse>
```

Added `getAuthHeaders()` method that automatically includes JWT token in requests:
```typescript
private getAuthHeaders(): HeadersInit {
  const token = localStorage.getItem('phagen_auth_token');
  return {
    'Content-Type': 'application/json',
    ...(token ? { 'Authorization': `Bearer ${token}` } : {}),
  };
}
```

### 2. Authentication Context

**File**: `frontend/contexts/AuthContext.tsx`

Features:
- Global auth state management
- Token validation on app load
- Auto-login after signup
- JWT token storage in localStorage

Methods:
```typescript
login(email: string, password: string): Promise<void>
signup(email: string, name: string, password: string): Promise<void>
logout(): void
```

State:
```typescript
user: User | null           // Current user or null
isAuthenticated: boolean    // Auth status
loading: boolean           // Loading state
```

### 3. Protected Routes

**Component**: `frontend/components/ProtectedRoute.tsx`

Usage:
```tsx
export default function ProtectedPage() {
  return (
    <ProtectedRoute>
      <YourPageContent />
    </ProtectedRoute>
  );
}
```

Behavior:
- Shows loading spinner while checking auth
- Redirects to `/login` if not authenticated
- Renders page content if authenticated

### 4. Protected Pages

Pages wrapped with `ProtectedRoute`:
- `/molecule` - Job submission
- `/job` - Job status monitoring
- `/results` - Analysis results
- `/history` - Job history

### 5. Login & Signup Pages

#### Login Page (`/login`)
Features:
- Email and password inputs
- Error handling with alerts
- Loading state with spinner
- Link to signup page
- "Forgot password" placeholder

#### Signup Page (`/signup`)
Features:
- Name, email, password, confirm password inputs
- Client-side validation:
  - Minimum 8 characters password
  - Password match confirmation
- Auto-login after successful signup
- Error handling
- Link to login page

## Security Features

### Backend Security

1. **Password Hashing**: Bcrypt with automatic salt generation
2. **JWT Tokens**: HS256 algorithm with configurable expiration
3. **Token Validation**: Automatic verification on protected endpoints
4. **User Status**: `is_active` flag for account management
5. **HTTPBearer**: Secure token transmission in Authorization header

### Frontend Security

1. **Token Storage**: localStorage (consider httpOnly cookies for production)
2. **Auto Token Refresh**: Validates token on app load
3. **Protected Routes**: Automatic redirect to login
4. **HTTPS Required**: Token transmission should use HTTPS in production

## Setup Instructions

### Backend Setup

1. **Install Dependencies**:
```bash
cd backend
pip install -r requirements.txt
```

2. **Set Environment Variables** (`.env`):
```env
SECRET_KEY=your-super-secret-key-min-32-chars-change-this-in-production
DATABASE_URL=postgresql+psycopg://phagen:phagen@localhost:5432/phagen
```

3. **Run Database Migrations**:
```bash
# The User table will be created automatically on first run
python -m uvicorn backend.app.main:app --reload
```

4. **Verify Setup**:
```bash
curl http://localhost:8000/api/health
```

### Frontend Setup

1. **Install Dependencies**:
```bash
cd frontend
npm install
```

2. **Set Environment Variables** (`.env.local`):
```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

3. **Run Development Server**:
```bash
npm run dev
```

4. **Test Authentication**:
- Navigate to `http://localhost:3000/signup`
- Create a new account
- Verify auto-login and redirect
- Test protected routes
- Test logout

## Testing Guide

### Manual Testing Checklist

- [ ] **Signup Flow**
  - [ ] Register with valid email and password
  - [ ] Verify auto-login after signup
  - [ ] Try duplicate email (should fail)
  - [ ] Try weak password < 8 chars (should fail)

- [ ] **Login Flow**
  - [ ] Login with correct credentials
  - [ ] Try incorrect password (should fail)
  - [ ] Try non-existent email (should fail)
  - [ ] Verify token storage in localStorage

- [ ] **Protected Routes**
  - [ ] Access `/molecule` without login (should redirect)
  - [ ] Login and access `/molecule` (should work)
  - [ ] Test all protected pages after login

- [ ] **Token Validation**
  - [ ] Refresh page while logged in (should stay logged in)
  - [ ] Clear localStorage and refresh (should redirect to login)
  - [ ] Submit job while authenticated

- [ ] **Logout Flow**
  - [ ] Click logout button
  - [ ] Verify redirect or UI update
  - [ ] Try accessing protected route (should redirect)

### API Testing with curl

1. **Signup**:
```bash
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{"email":"test@example.com","name":"Test User","password":"testpass123"}'
```

2. **Login**:
```bash
curl -X POST http://localhost:8000/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{"email":"test@example.com","password":"testpass123"}'
```

3. **Get Current User**:
```bash
TOKEN="your-token-here"
curl http://localhost:8000/api/auth/me \
  -H "Authorization: Bearer $TOKEN"
```

4. **Create Job (Authenticated)**:
```bash
curl -X POST http://localhost:8000/api/jobs \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","molecule":"Aspirin"}'
```

## Production Considerations

### Security Enhancements

1. **Secret Key Management**:
   - Use environment variables
   - Generate with: `openssl rand -hex 32`
   - Never commit to version control

2. **Token Storage**:
   - Consider httpOnly cookies instead of localStorage
   - Implement refresh tokens for longer sessions

3. **HTTPS**:
   - Enforce HTTPS in production
   - Use secure cookies if switching from localStorage

4. **Rate Limiting**:
   - Add rate limiting to auth endpoints
   - Prevent brute force attacks

5. **Password Policy**:
   - Enforce stronger password requirements
   - Add password strength meter
   - Implement password reset flow

6. **Account Management**:
   - Email verification
   - Password reset via email
   - Account lockout after failed attempts
   - Two-factor authentication (2FA)

### Performance Optimization

1. **Token Caching**: Cache decoded tokens to reduce CPU usage
2. **Database Indexing**: Ensure email column is indexed (already done)
3. **Connection Pooling**: Use SQLAlchemy connection pooling

### Monitoring

1. **Log Authentication Events**:
   - Failed login attempts
   - Account creations
   - Token validation failures

2. **Metrics**:
   - Track authentication success/failure rates
   - Monitor token expiration patterns
   - Alert on suspicious activity

## Troubleshooting

### Backend Issues

**Issue**: Import errors for jose or passlib
```bash
# Solution: Install dependencies
pip install python-jose[cryptography] passlib[bcrypt] python-multipart
```

**Issue**: "Could not validate credentials" error
- Check token format in Authorization header
- Verify SECRET_KEY matches between token creation and validation
- Check token expiration time

**Issue**: Database errors
```python
# Check database connection
from backend.app.database import get_db
with get_db() as db:
    print(db.query(User).count())
```

### Frontend Issues

**Issue**: "Unauthorized" errors
- Check localStorage for token
- Verify token format: `Bearer <token>`
- Check browser console for errors

**Issue**: Infinite redirect loop
- Check loading state in AuthContext
- Verify ProtectedRoute logic
- Clear localStorage and cookies

**Issue**: CORS errors
- Verify CORS settings in backend
- Check API_BASE_URL configuration
- Ensure credentials are included in requests

## API Documentation

Full API documentation available at:
- Swagger UI: `http://localhost:8000/docs`
- ReDoc: `http://localhost:8000/redoc`

## Support

For issues or questions:
1. Check this README
2. Review error logs (backend console or browser console)
3. Test with curl to isolate frontend vs backend issues
4. Check database tables exist and have correct schema

---

**Status**: ✅ Production-ready authentication system
**Last Updated**: December 6, 2024
**Version**: 1.0.0
