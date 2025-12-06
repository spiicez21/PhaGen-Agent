@echo off
REM PhaGen-Agent Startup Script
REM Opens each service in a separate PowerShell window

echo Starting PhaGen-Agent services...
echo.

REM Check if Docker Desktop is running
echo Checking Docker Desktop...
docker info >nul 2>&1
if errorlevel 1 (
    echo [WARNING] Docker Desktop is not running!
    echo Please start Docker Desktop first, then run this script again.
    echo.
    echo Skipping MinIO and continuing with other services...
    echo.
    timeout /t 3 /nobreak >nul
    goto SKIP_MINIO
)

REM Start MinIO (Docker)
echo [1/4] Starting MinIO S3 Storage...
start "MinIO S3" powershell -NoExit -Command "cd D:\PhaGen-Agent\infra; docker compose up minio"
timeout /t 3 /nobreak >nul

:SKIP_MINIO

REM Start Ollama (if not already running)
echo [2/4] Starting Ollama LLM Server...
start "Ollama LLM" powershell -NoExit -Command "ollama serve"

REM Wait a bit for Ollama to start
timeout /t 3 /nobreak >nul

REM Start Backend API
echo [3/4] Starting FastAPI Backend...
start "FastAPI Backend" powershell -NoExit -Command "cd D:\PhaGen-Agent\backend; & D:/PhaGen-Agent/.venv/Scripts/Activate.ps1; python -m uvicorn app.main:app --reload"

REM Wait a bit for backend to start
timeout /t 5 /nobreak >nul

REM Start Frontend
echo [4/4] Starting Next.js Frontend...
start "Next.js Frontend" powershell -NoExit -Command "cd D:\PhaGen-Agent\frontend; npm run dev"

echo.
echo ========================================
echo All services started in separate windows!
echo ========================================
echo.
echo Services:
echo   - MinIO S3:     http://localhost:9000 (console: http://localhost:9001)
echo   - Ollama LLM:   http://localhost:11434
echo   - Backend API:  http://localhost:8000
echo   - Frontend UI:  http://localhost:3000
echo.
echo Press any key to exit this window...
pause >nul
