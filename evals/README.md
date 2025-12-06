# Evaluation Suite

This folder carries a lightweight regression harness so we can spot check the
agent stack against a known set of sample molecules. The suite currently ships
with **five** molecules (Pirfenidone, Metformin, Nintedanib, Lenabasum,
Pamrevlumab) and can be extended as we add more benchmarks.

## Layout

```
/ evals
  ├─ cases.json           # master list of molecules + expectations
  ├─ fixtures/*.json      # captured payloads from trusted runs
  └─ run_eval.py          # CLI harness (fixtures or live API)
```

Each entry inside `cases.json` describes:

- the molecule input (name, synonyms, optional SMILES)
- a pointer to its fixture payload for offline checks
- guardrail expectations (recommendation label, acceptable market score band,
  minimum evidence counts per worker, acceptable quality status, minimum story
  length)

## Running against fixtures

Fixtures allow everyone to validate the harness without standing up the full
stack. From the repo root:

```powershell
python evals/run_eval.py
```

Use `--case pirfenidone` (repeatable) to run a subset, or `--cases-file` if you
clone the suite elsewhere.

## Running against the live API

Point the harness at a running FastAPI instance (e.g., `uvicorn app.main:app` or
the Docker Compose stack). The tool will POST each molecule, poll
`/api/jobs/{id}`, and compare the finished payload to the expectations in
`cases.json`.

```powershell
python evals/run_eval.py --mode live --api-base http://localhost:8000/api
```

Tune execution with:

- `--poll-interval` (seconds between GET requests, default 3s)
- `--timeout` (overall job deadline, default 180s)

## Extending the suite

1. Capture a new payload (download from the UI or directly from `/api/jobs`).
2. Drop it under `evals/fixtures/<molecule>.json` (either the raw payload or a
   wrapper `{ "payload": { ... } }`).
3. Add an entry to `cases.json` with reasonable expectations (score bands,
   required workers, etc.).
4. Commit both files. CI/contributors can now run `python evals/run_eval.py`
   before shipping guardrail-sensitive changes.
