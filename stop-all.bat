@echo off
REM PhaGen-Agent Shutdown Script
REM Stops all services gracefully

echo Stopping PhaGen-Agent services...
echo.

REM Stop Docker MinIO
echo [1/3] Stopping MinIO S3 Storage...
cd D:\PhaGen-Agent\infra
docker compose down
echo.

REM Kill Python/Uvicorn processes (backend)
echo [2/3] Stopping Backend API...
taskkill /F /FI "WINDOWTITLE eq FastAPI Backend*" 2>nul
taskkill /F /IM uvicorn.exe 2>nul
echo.

REM Kill Node/Next.js processes (frontend)
echo [3/3] Stopping Frontend...
taskkill /F /FI "WINDOWTITLE eq Next.js Frontend*" 2>nul
for /f "tokens=5" %%a in ('netstat -aon ^| findstr :3000') do taskkill /F /PID %%a 2>nul
echo.

echo ========================================
echo All services stopped!
echo ========================================
echo.
echo Note: Ollama server is still running (started independently)
echo To stop Ollama: taskkill /F /IM ollama.exe
echo.
pause
