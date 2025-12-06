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
3. **Supabase project + connection string**. Grab the `psycopg` URI from **Supabase → Project Settings → Database → Connection string** and paste it into `.env` as `SUPABASE_URL=postgresql://...`.

Docker Desktop is optional—install it only if you plan to run the ancillary services in `infra/docker-compose.yml` (MinIO, rdkit-service, etc.).

Re-activate the virtualenv in every new shell with `& D:/PhaGen-Agent/.venv/Scripts/Activate.ps1`.

## 2. Configure the database (Supabase first)

1. Edit `.env` in the repo root and make sure the `SUPABASE_URL` line contains your hosted Postgres URI:
   ```dotenv
   SUPABASE_URL=postgresql://<user>:<password>@<host>:5432/postgres
   ```
2. (Optional) Export it inline if you need to override during experiments:
   ```powershell
   $env:SUPABASE_URL = "postgresql://..."
   ```
3. Test the connection quickly (optional):
   ```powershell
   $env:SUPABASE_URL = (Get-Content .\.env | Select-String -Pattern '^SUPABASE_URL=').Line.Split('=', 2)[1].Trim()
   python -c "import os, psycopg; psycopg.connect(os.environ['SUPABASE_URL']).close(); print('Supabase reachable ✅')"
   ```

If Supabase is unavailable, point the backend at SQLite temporarily with `DATABASE_URL=sqlite:///./phagen.db`.

## 3. Start the backend API (FastAPI + Uvicorn)

```powershell
Set-Location D:\PhaGen-Agent
& .venv\Scripts\Activate.ps1
Set-Location .\backend
$env:STRUCTURE_CATALOG_PATH = "D:\PhaGen-Agent\indexes\data\structures\structures.manifest.json"
$env:STORAGE_PROVIDER = "none"  # optional: keep uploads local during dev
python -m uvicorn app.main:app --host 127.0.0.1 --port 8000 --reload
```

Notes:
- The backend reads `.env` values automatically; export env vars inline when experimenting.
- To point the API at hosted Supabase, set `SUPABASE_URL` (or `DATABASE_URL`) to your Supabase connection string—the backend will pick up whichever is present.
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
