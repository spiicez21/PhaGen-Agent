# Testing Guide: Automatic Navigation

## Quick Test Steps

### Test 1: Normal Flow (Most Important)
1. Start the backend: `cd backend && .venv\Scripts\activate && uvicorn app.main:app --reload`
2. Start the frontend: `cd frontend && npm run dev`
3. Login and go to `/molecule`
4. Enter SMILES: `CC(=O)OC1=CC=CC=C1C(=O)O` (Aspirin)
5. Submit
6. **Expected:** Redirects to `/job?id=XXX` automatically
7. Watch the job status page update every 3 seconds
8. **Expected:** When complete, auto-redirects to `/results?id=XXX` after 1.5 seconds

### Test 2: Header Navigation
1. With a job running, click "Job Status" in header
2. **Expected:** Automatically loads the running job
3. Open DevTools Console (F12)
4. **Expected:** See `[Job Status] Auto-detection: Found X jobs in history`

### Test 3: Direct URL Access
1. Navigate to `/job` (no ID)
2. **Expected:** 
   - If there's a running/pending job → shows that job
   - If only completed jobs → redirects to `/results`
   - If no jobs → shows error with helpful buttons

### Test 4: Results Auto-Detection
1. Navigate to `/results` (no ID)
2. **Expected:** Automatically loads most recent completed job
3. Check console for: `[Results] Auto-detected completed job: XXX`

### Test 5: History Page Links
1. Go to `/history`
2. Click "View Results" on any completed job
3. **Expected:** Shows results with correct job ID

## Console Debugging

Open DevTools Console and look for these messages:

### Job Status Page
- ✅ `[Job Status] Auto-detection: Found X jobs in history`
- ✅ `[Job Status] Auto-detected ongoing job: XXX`
- ✅ `[Job Status] Poll result: pending/running/completed for job XXX`
- ✅ `[Job Status] Job completed, redirecting to results...`

### Results Page
- ✅ `[Results] Auto-detection: Found X jobs in history`
- ✅ `[Results] Auto-detected completed job: XXX`

### Errors to Watch For
- ❌ `Failed to fetch job`
- ❌ `No active jobs found`
- ❌ `No completed jobs found`

## localStorage Inspection

Check job history in browser console:
```javascript
// View all stored jobs
JSON.parse(localStorage.getItem('phagen_job_history'))

// Clear all jobs (for testing)
localStorage.removeItem('phagen_job_history')

// Check specific job
JSON.parse(localStorage.getItem('phagen_job_history'))
  .find(j => j.job_id === 'YOUR_JOB_ID')
```

## Common Issues

### Issue: "No active jobs found" immediately after submission
**Solution:** Check if jobStore.add() is being called in molecule/page.tsx

### Issue: Job status page doesn't redirect when complete
**Solution:** 
1. Check console for polling messages
2. Verify job.status === "COMPLETED" in API response (uppercase)
3. Frontend should convert to lowercase before comparison
4. Check if router.push() is being called

### Issue: Results page shows "Job is still running"
**Solution:** Wait for job to complete, or check job status first

### Issue: Results not displaying even though job is COMPLETED
**Solution:** Case-sensitivity bug - fixed by converting status to lowercase before comparison. Verify the fix is applied:
```typescript
const status = job.status.toLowerCase();
if (status === "completed") { ... }
```

## Expected Polling Behavior

1. Job submitted → status: "PENDING" (backend returns uppercase)
2. Agents start → status: "RUNNING"
3. Poll every 3 seconds
4. Workers complete one by one
5. All workers done → status: "COMPLETED"
6. Frontend converts to lowercase for comparison
7. Wait 1.5 seconds → redirect to results

## Browser DevTools Tips

1. **Network Tab:** Watch for `GET /api/jobs/{id}` requests every 3 seconds
2. **Console:** Look for `[Job Status]` and `[Results]` logs
3. **Application Tab → Local Storage:** Check `phagen_job_history`
4. **React DevTools:** Inspect jobId state in components
5. **Status Values:** Backend returns UPPERCASE (`COMPLETED`, `RUNNING`, `PENDING`, `FAILED`)
6. **Frontend Handling:** Status is converted to lowercase before comparison
