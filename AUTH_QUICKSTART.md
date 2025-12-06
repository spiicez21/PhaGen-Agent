# Authentication System - Quick Start Guide

## 🎉 Implementation Complete!

A full-stack JWT-based authentication system has been successfully implemented for PhaGen Agent.

## 📋 What Was Implemented

### Backend (FastAPI + PostgreSQL)
- ✅ User database model with bcrypt password hashing
- ✅ JWT token generation and validation
- ✅ Authentication endpoints (`/api/auth/signup`, `/api/auth/login`, `/api/auth/me`)
- ✅ Protected route dependencies for securing endpoints
- ✅ SQLAlchemy database integration

### Frontend (Next.js + React)
- ✅ Authentication context with React hooks
- ✅ Login and signup pages with validation
- ✅ Protected route component with auto-redirect
- ✅ API client with automatic JWT token injection
- ✅ Token persistence in localStorage
- ✅ Protected pages: molecule, job, results, history

## 🚀 Quick Start

### 1. Install Backend Dependencies

```bash
cd backend
pip install -r requirements.txt
```

Key dependencies installed:
- `python-jose[cryptography]` - JWT tokens
- `passlib[bcrypt]` - Password hashing
- `python-multipart` - Form parsing

### 2. Start Backend Server

```bash
cd backend
python -m uvicorn app.main:app --reload
```

The server will:
- Create the `users` table automatically on first run
- Start on `http://localhost:8000`
- Expose auth endpoints at `/api/auth/*`

### 3. Start Frontend

```bash
cd frontend
npm run dev
```

Frontend will run on `http://localhost:3000`

## 🧪 Testing the System

### Test Flow 1: New User Signup

1. Navigate to `http://localhost:3000/signup`
2. Fill in the form:
   - Name: John Doe
   - Email: john@example.com
   - Password: testpass123 (min 8 characters)
   - Confirm Password: testpass123
3. Click "Create Account"
4. You'll be auto-logged in and redirected to home page
5. Header should show "John Doe" and logout button

### Test Flow 2: Login

1. Navigate to `http://localhost:3000/login`
2. Enter credentials:
   - Email: john@example.com
   - Password: testpass123
3. Click "Sign In"
4. You'll be redirected to home page
5. Check browser localStorage for `phagen_auth_token`

### Test Flow 3: Protected Routes

1. Open incognito window (or logout)
2. Try to access `http://localhost:3000/molecule`
3. You should be redirected to `/login`
4. Login and try again - you should see the page

### Test Flow 4: Create Job (Authenticated)

1. Login to the system
2. Navigate to `/molecule`
3. Enter SMILES: `CC(=O)OC1=CC=CC=C1C(=O)O`
4. Click "Run Analysis"
5. Job should be created successfully

## 📁 Files Changed

### Backend Files Created/Modified

| File | Status | Description |
|------|--------|-------------|
| `backend/requirements.txt` | ✏️ Modified | Added auth dependencies |
| `backend/app/config.py` | ✏️ Modified | Added JWT config (secret_key, algorithm, expiration) |
| `backend/app/models.py` | ✏️ Modified | Added User model |
| `backend/app/database.py` | ✏️ Modified | Added get_db() dependency |
| `backend/app/schemas.py` | ✏️ Modified | Added auth schemas (UserCreate, UserLogin, Token, etc.) |
| `backend/app/auth.py` | ✨ Created | JWT and password utilities |
| `backend/app/routers/auth.py` | ✨ Created | Auth endpoints (/signup, /login, /me) |
| `backend/app/main.py` | ✏️ Modified | Registered auth router |

### Frontend Files Created/Modified

| File | Status | Description |
|------|--------|-------------|
| `frontend/lib/api.ts` | ✏️ Modified | Added auth methods, auto JWT injection |
| `frontend/contexts/AuthContext.tsx` | ✏️ Modified | Real API integration, token validation |
| `frontend/components/ProtectedRoute.tsx` | ✨ Created | Route protection component |
| `frontend/app/login/page.tsx` | ✏️ Modified | Updated to use real login API |
| `frontend/app/signup/page.tsx` | ✏️ Modified | Updated to use real signup API |
| `frontend/app/molecule/page.tsx` | ✏️ Modified | Wrapped with ProtectedRoute |
| `frontend/app/job/page.tsx` | ✏️ Modified | Wrapped with ProtectedRoute |
| `frontend/app/results/page.tsx` | ✏️ Modified | Wrapped with ProtectedRoute |
| `frontend/app/history/page.tsx` | ✏️ Modified | Wrapped with ProtectedRoute |

## 🔐 Security Features

### Backend Security
- ✅ Bcrypt password hashing (automatic salt)
- ✅ JWT tokens with expiration (30 minutes default)
- ✅ HS256 algorithm for token signing
- ✅ HTTPBearer token authentication
- ✅ User account status (`is_active` flag)
- ✅ Token validation on every protected request

### Frontend Security
- ✅ Client-side token validation on app load
- ✅ Automatic token expiry handling
- ✅ Protected route redirects
- ✅ No password storage (only hashed on backend)
- ✅ Authorization header with Bearer token

## 🔧 Configuration

### Backend Environment Variables

Create/update `.env` in project root:

```env
# IMPORTANT: Change this in production!
SECRET_KEY=your-super-secret-key-at-least-32-characters-long

# Database
DATABASE_URL=postgresql+psycopg://phagen:phagen@localhost:5432/phagen

# JWT Configuration (optional, has defaults)
ACCESS_TOKEN_EXPIRE_MINUTES=30
ALGORITHM=HS256
```

### Frontend Environment Variables

Create `.env.local` in frontend directory:

```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

## 📡 API Endpoints

### Authentication Endpoints

| Method | Endpoint | Auth Required | Description |
|--------|----------|---------------|-------------|
| POST | `/api/auth/signup` | No | Register new user |
| POST | `/api/auth/login` | No | Login and get token |
| GET | `/api/auth/me` | Yes | Get current user info |

### Example Requests

**Signup**:
```bash
curl -X POST http://localhost:8000/api/auth/signup \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "name": "Test User",
    "password": "securepass123"
  }'
```

**Login**:
```bash
curl -X POST http://localhost:8000/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "securepass123"
  }'
```

**Get Current User** (with token):
```bash
TOKEN="your-jwt-token-here"
curl http://localhost:8000/api/auth/me \
  -H "Authorization: Bearer $TOKEN"
```

## 🐛 Troubleshooting

### Issue: "Import jose could not be resolved"

**Solution**: Install backend dependencies
```bash
cd backend
pip install -r requirements.txt
```

### Issue: Login redirects back to login page

**Causes**:
1. Check browser console for errors
2. Verify backend is running on port 8000
3. Check CORS settings
4. Verify token is being stored in localStorage

**Debug**:
```javascript
// In browser console
console.log(localStorage.getItem('phagen_auth_token'));
```

### Issue: "Could not validate credentials"

**Causes**:
1. Token expired (default: 30 minutes)
2. SECRET_KEY changed between token creation and validation
3. Malformed token

**Solution**: Logout and login again to get fresh token

### Issue: Database errors on startup

**Solution**: Ensure PostgreSQL is running and database exists
```bash
# Check if database exists
psql -U phagen -d phagen -c "SELECT 1"

# If not, create it
createdb -U phagen phagen
```

## 📚 Documentation

- **Full Documentation**: See `AUTH_IMPLEMENTATION.md`
- **API Docs**: `http://localhost:8000/docs` (Swagger UI)
- **Alternative API Docs**: `http://localhost:8000/redoc`

## 🎯 Next Steps

### Recommended Enhancements

1. **Email Verification**
   - Send verification email on signup
   - Verify email before allowing login

2. **Password Reset**
   - "Forgot Password" functionality
   - Email-based password reset flow

3. **Refresh Tokens**
   - Long-lived refresh tokens
   - Short-lived access tokens

4. **Rate Limiting**
   - Prevent brute force attacks
   - Limit signup and login attempts

5. **Two-Factor Authentication (2FA)**
   - TOTP-based 2FA
   - SMS verification

6. **Admin Dashboard**
   - User management
   - Account status control
   - Activity monitoring

### Production Checklist

- [ ] Change `SECRET_KEY` to secure random value (min 32 chars)
- [ ] Enable HTTPS
- [ ] Use httpOnly cookies instead of localStorage
- [ ] Implement refresh token rotation
- [ ] Add rate limiting
- [ ] Set up monitoring and alerting
- [ ] Enable CORS only for production domains
- [ ] Implement proper error logging
- [ ] Set up automated backups
- [ ] Add comprehensive testing (unit + integration)

## ✅ System Status

**Backend**: ✅ Ready for testing  
**Frontend**: ✅ Ready for testing  
**Database**: ✅ Auto-creates tables  
**Authentication**: ✅ Fully functional  
**Protected Routes**: ✅ Working  

## 🎓 Key Concepts

### JWT Tokens
- Stateless authentication (no session storage needed)
- Contains user info (email) in payload
- Signed with SECRET_KEY to prevent tampering
- Expires after configured time (30 min default)

### Password Security
- Never stored in plain text
- Bcrypt hashing with automatic salt
- One-way encryption (cannot be reversed)
- Verified by hashing input and comparing

### Protected Routes
- Frontend: React component wrapper
- Backend: FastAPI dependency injection
- Automatic redirect to login if unauthorized
- Token passed in Authorization header

---

**Ready to use!** 🚀  
**Questions?** Check `AUTH_IMPLEMENTATION.md` for detailed documentation.
