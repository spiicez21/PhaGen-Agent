# Frontend-Backend Integration

## Overview

The PhaGen frontend is now fully connected to the backend API. All mock data has been removed and replaced with real API calls.

## API Configuration

### Environment Variables

Create a `.env.local` file in the `frontend` directory:

```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api/v1
```

### API Client

The API client is located at `lib/api.ts` and provides methods for:

- **createJob(request)**: Submit a new molecule analysis
- **getJob(jobId)**: Get job status and results
- **getJobsForComparison(jobIds)**: Get multiple jobs for comparison
- **downloadReport(jobId)**: Download PDF report
- **getHealth()**: Check API health
- **submitFeedback(data)**: Submit user feedback
- **getJobFeedback(jobId)**: Get feedback for a job

## Updated Pages

### 1. `/molecule` - Job Submission
- **Changed**: Removed mock form fields, simplified to SMILES + optional molecule name and query
- **API Call**: `POST /api/v1/jobs` with `{ smiles, molecule?, query? }`
- **Response**: Returns `job_id` and redirects to `/job?id={job_id}`

### 2. `/job` - Job Status
- **Changed**: Client-side polling every 3 seconds to check job status
- **API Call**: `GET /api/v1/jobs/{job_id}`
- **Behavior**: 
  - Shows loading spinner while job is pending/running
  - Auto-redirects to results page when completed
  - Shows error message if job fails

### 3. `/results` - Analysis Results
- **Changed**: Fetches real job data from backend
- **API Call**: `GET /api/v1/jobs/{job_id}`
- **Displays**: 
  - Molecule name, recommendation, innovation story
  - Worker results (clinical, literature, patent, market)
  - Quality metrics and validation data
  - Structure visualization

## Backend Requirements

### Start Backend Server

```bash
cd backend
uvicorn app.main:app --reload --port 8000
```

### Backend Endpoints Used

- `POST /api/v1/jobs` - Create new analysis job
- `GET /api/v1/jobs/{job_id}` - Get job status/results
- `GET /api/v1/jobs/compare?job_ids=...` - Compare multiple jobs
- `GET /api/v1/jobs/{job_id}/report.pdf` - Download PDF report
- `GET /api/v1/health` - Health check
- `POST /api/v1/feedback` - Submit feedback
- `GET /api/v1/feedback/job/{job_id}` - Get job feedback

### CORS Configuration

Backend is configured to allow requests from:
- `http://localhost:3000`
- `http://127.0.0.1:3000`

## Development Workflow

1. Start backend server:
   ```bash
   cd backend
   uvicorn app.main:app --reload --port 8000
   ```

2. Start frontend dev server:
   ```bash
   cd frontend
   npm run dev
   ```

3. Open browser to `http://localhost:3000`

## Job Flow

1. User submits SMILES on `/molecule` page
2. Frontend calls `POST /api/v1/jobs` with request data
3. Backend creates job and returns `job_id`
4. Frontend redirects to `/job?id={job_id}`
5. Job page polls `GET /api/v1/jobs/{job_id}` every 3 seconds
6. When status becomes "completed", auto-redirects to `/results?id={job_id}`
7. Results page fetches final job data and displays analysis

## Removed Mock Data

The following mock data files are no longer used:
- `app/sample-data.ts` - Still exists but not imported in connected pages
- Mock data constants like `SAMPLE_PAYLOAD`, `DEMO_JOB`, `HISTORY_RUNS`

## Error Handling

All API calls include proper error handling:
- Network errors show user-friendly error messages
- Failed jobs display error details
- Loading states with spinners during API calls
- Graceful fallbacks for missing data

## Production Deployment

For production, update the API URL in `.env.local`:

```env
NEXT_PUBLIC_API_URL=https://your-production-api.com/api/v1
```

Make sure backend CORS settings include production frontend URL.
