# Local Startup Cheat Sheet

Fast way to spin up every moving part (database, backend API, frontend UI) on a Windows workstation. Run commands from a PowerShell prompt unless noted otherwise.

## 1. Prerequisites (once per machine)

1. **Python virtual environment** shared by backend + index tools:
   ```powershell
   Set-Location D:\PhaGen-Agent
   python -m venv .venv
   .\.venv\Scripts\activate
   pip install -r requirements.txt
   ```
2. **Node/NPM deps for the frontend** (run inside `frontend/` after cloning):
   ```powershell
   Set-Location D:\PhaGen-Agent\frontend
   npm install
   ```
3. **Docker Desktop** running locally so `docker compose` can start Postgres (and optional MinIO/Ollama/etc.).

Re-activate the virtualenv in every new shell with `& D:/PhaGen-Agent/.venv/Scripts/Activate.ps1`.

## 2. Start the database (Postgres via Docker)

```powershell
Set-Location D:\PhaGen-Agent\infra
$env:COMPOSE_PROJECT_NAME = "phagen"
docker compose up postgres -d
```

- Wait for the `postgres` container to report `healthy`. Check logs with `docker compose logs -f postgres`.
- When finished for the day: `docker compose down postgres` (still from `infra/`).
- If Docker is unavailable, point the backend at SQLite (`DATABASE_URL=sqlite:///./phagen.db`) and skip this section.

## 3. Start the backend API (FastAPI + Uvicorn)

```powershell
Set-Location D:\PhaGen-Agent
& .venv\Scripts\Activate.ps1
Set-Location .\backend
$env:STRUCTURE_CATALOG_PATH = "D:\PhaGen-Agent\indexes\data\structures\structures.manifest.json"
$env:DATABASE_URL = "postgresql+psycopg://phagen:phagen@localhost:5432/phagen"  # swap for sqlite:///./phagen.db when Postgres is offline
$env:STORAGE_PROVIDER = "none"  # optional: keep uploads local during dev
python -m uvicorn app.main:app --host 127.0.0.1 --port 8000 --reload
```

Notes:
- The backend reads `.env` values automatically; export env vars inline when experimenting.
- `STRUCTURE_CATALOG_PATH` should point at the latest manifest created by `indexes/build_index.py`.
- If you need MinIO/S3, set `S3_*` env vars before running uvicorn.

## 4. Start the frontend (Next.js)

In a second terminal:

```powershell
Set-Location D:\PhaGen-Agent\frontend
npm run dev
```

- App serves at `http://localhost:3000` by default.
- The UI expects the backend at `http://127.0.0.1:8000`; adjust `.env.local` if you proxy through a different port.

## 5. Optional helpers

- **Indexes / RDKit refresh**: `python indexes\build_index.py --structure-records indexes\data\normalized_smiles.jsonl`
- **Crawler seed run**: `cd crawler; npm run crawl`
- **Stop everything**: Ctrl+C in frontend/backend terminals, then `docker compose down` in `infra/`.

Keep this file close when onboarding teammates or rebooting the stack after a long break. Update it whenever startup commands change.
