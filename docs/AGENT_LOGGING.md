# Agent Logging Configuration

## Overview
Agent logging has been configured to display detailed execution information in the backend console. This helps track agent progress, debug issues, and monitor performance in real-time.

## What Was Added

### 1. Backend Main Configuration (`backend/app/main.py`)
- Configured `logging.basicConfig()` with INFO level
- Added formatted output: `timestamp - module - level - message`
- Enabled agent-specific loggers:
  - `agents.master` - Master agent orchestration
  - `agents.llm` - LLM API calls and responses
  - `agents.workers` - Worker execution details
  - `agents.retrieval` - Knowledge base retrieval

### 2. Master Agent Logging (`agents/master.py`)
Added comprehensive logging at key stages:
- **Job Start**: Displays molecule name, synonyms, and SMILES
- **Synonym Expansion**: Shows expanded synonym list
- **Worker Execution**: Logs each worker start and completion with confidence scores
- **Worker Failures**: Error logging with failure reasons
- **Synthesis Phase**: Market score, fallback generation, final synthesis
- **Completion**: Success summary with story length and worker statistics

### 3. Worker Base Logging (`agents/workers/base.py`)
Added detailed execution tracking:
- **Query Building**: Number of search queries generated
- **Retrieval**: Passages retrieved from knowledge base
- **Fallback**: When fallback queries are used
- **LLM Processing**: Summary generation status
- **Completion**: Final confidence scores

### 4. Job Router Logging (`backend/app/routers/jobs.py`)
Enhanced job execution visibility:
- **Job Start**: Clear job ID and molecule identification
- **Job Success**: Recommendation and report version
- **Job Failure**: Error details with full traceback

## Log Format

```
2025-12-06 14:11:15,261 - agents.master - INFO - [START] Agent analysis for molecule: Nintedanib
2025-12-06 14:11:15,262 - agents.master - INFO -    Provided synonyms: ['BIBF 1120']
2025-12-06 14:11:15,263 - agents.master - INFO -    SMILES: COC(=O)c1ccc(c(c1)NC(=O)c2ccc(cc2)CN3CCN(CC3)C)Nc4nc5ccccc5c(n4)C
2025-12-06 14:11:15,264 - agents.master - INFO -    Expanded synonyms: ['Nintedanib', 'BIBF 1120', 'Ofev']

2025-12-06 14:11:15,265 - agents.master - INFO - 
[WORKERS] Executing 4 workers: ['clinical', 'patent', 'literature', 'market']

2025-12-06 14:11:15,266 - agents.master - INFO -    [RUN] Running clinical worker...
2025-12-06 14:11:15,267 - agents.workers.base - INFO -       [QUERY] [clinical] Building search queries for: Nintedanib
2025-12-06 14:11:15,268 - agents.workers.base - INFO -       [QUERY] [clinical] Generated 3 search queries
2025-12-06 14:11:16,789 - agents.workers.base - INFO -       [RETRIEVAL] [clinical] Retrieved 10 passages from knowledge base
2025-12-06 14:11:16,790 - agents.workers.base - INFO -       [LLM] [clinical] Generating summary with LLM...
2025-12-06 14:11:18,234 - agents.workers.base - INFO -       [OK] [clinical] Summary complete (confidence: 0.85)
2025-12-06 14:11:18,235 - agents.master - INFO -    [OK] clinical worker completed (confidence: 0.85)

[... similar logs for other workers ...]

2025-12-06 14:11:45,123 - agents.master - INFO - 
[SCORE] Market score: 75

2025-12-06 14:11:45,124 - agents.master - INFO - 
[FALLBACK] Generating fallback story and recommendation...

2025-12-06 14:11:45,125 - agents.master - INFO - 
[SYNTHESIS] Synthesizing final innovation story and recommendation...
2025-12-06 14:11:47,456 - agents.master - INFO -    Final recommendation: Approved

2025-12-06 14:11:47,457 - agents.master - INFO - 
[COMPLETE] Agent analysis completed successfully for Nintedanib
2025-12-06 14:11:47,458 - agents.master - INFO -    Story length: 2847 chars
2025-12-06 14:11:47,459 - agents.master - INFO -    Workers completed: 4/4

================================================================================
[SUCCESS] Job abc123 completed successfully
   Recommendation: Approved
   Report version: v1
================================================================================
```

## Log Tags
- `[START]` - Job/agent start
- `[WORKERS]` - List/collection operations
- `[RUN]` - Worker execution start
- `[QUERY]` - Search/query building
- `[RETRIEVAL]` - Knowledge base retrieval
- `[FALLBACK]` - Retry/fallback operations
- `[LLM]` - LLM/synthesis operations
- `[SCORE]` - Scoring/evaluation
- `[SYNTHESIS]` - Final story synthesis
- `[OK]` - Success/completion
- `[COMPLETE]` - Overall completion
- `[SUCCESS]` - Job success
- `[ERROR]` - Error/failure
- `[FAIL]` - Worker failure

## Benefits

1. **Real-time Progress Tracking**: See exactly what the agent is doing at each step
2. **Performance Monitoring**: Identify slow workers or retrieval operations
3. **Debugging**: Detailed error information with context
4. **Confidence Visibility**: Track confidence scores for each worker
5. **Production Monitoring**: Easy to spot issues in live deployments

## Testing

Run the test script to verify logging configuration:
```bash
python test_agent_logging.py
```

Start the backend to see live logs:
```bash
cd backend
python -m uvicorn app.main:app --reload
```

Then submit a job through the API and watch the console output!

## Configuration

To adjust log levels, modify `backend/app/main.py`:

```python
# For DEBUG level (very verbose)
logging.getLogger('agents.master').setLevel(logging.DEBUG)

# For WARNING level only (minimal output)
logging.getLogger('agents.master').setLevel(logging.WARNING)

# For ERROR level only
logging.getLogger('agents.master').setLevel(logging.ERROR)
```

## Notes

- Logs appear in the backend console (not browser)
- When using `--reload`, logs will restart on file changes
- Production deployments should use proper log aggregation (e.g., ELK stack, CloudWatch)
- Sensitive data (SMILES, molecule names) may appear in logs - secure accordingly
