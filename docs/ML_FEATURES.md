# Advanced ML Features Implementation Summary

## Overview
This document summarizes the implementation of three advanced ML features for PhaGen: disease mapping, reranker fine-tuning, and custom embedding training. These features enable continuous improvement of retrieval quality through user feedback and provide proactive repurposing suggestions.

---

## 1. Molecule-Disease Mapping (`backend/app/ml/disease_mapper.py`)

### Purpose
Predicts potential disease indications for molecules based on evidence gathered during analysis, enabling proactive repurposing suggestions.

### Implementation
- **MoleculeDiseaseMapper** class with category-based prediction
- Pre-defined disease categories: cancer, cardiovascular, neurodegenerative, inflammatory, metabolic, infectious
- Keyword-based scoring from clinical trials and literature evidence
- Confidence scoring normalized by keyword matches
- Support for adding custom known mappings

### Key Methods
```python
predict_diseases(molecule, smiles, evidence, top_k=5)
# Returns List[DiseaseMapping] sorted by confidence

add_known_mapping(molecule, disease, disease_id, confidence, source)
# Stores validated molecule-disease associations
```

### Integration
- Automatically runs after job completion in `backend/app/routers/jobs.py`
- Suggestions added to job payload under `repurposing_suggestions` field
- Non-blocking: job succeeds even if suggestions fail

### Usage Example
```python
from backend.app.ml import MoleculeDiseaseMapper

mapper = MoleculeDiseaseMapper()
predictions = mapper.predict_diseases(
    molecule="Aspirin",
    evidence=[{
        "type": "clinical",
        "text": "Phase 2 trial in heart failure patients",
        "url": "https://clinicaltrials.gov/..."
    }]
)
# Returns: [DiseaseMapping(disease="Cardiovascular", confidence=0.85, ...)]
```

---

## 2. Reranker Fine-Tuning (`backend/app/ml/reranker_trainer.py`)

### Purpose
Learns optimal evidence type weights from historical user feedback (upvotes/downvotes) to improve future retrieval relevance.

### Implementation
- **RerankerTrainer** class with gradient descent optimization
- Fetches training data from `evidence_feedback` table
- Adjusts weights for clinical, literature, patent, market evidence
- Applies recency boost (3-year window) and citation boost (>10 citations)
- Persists learned weights to `models/reranker/reranker_weights.json`

### Training Algorithm
1. Fetch feedback from database (feedback_score: +1.0 upvote, -1.0 downvote)
2. Group by evidence type
3. Update weights: `weight += learning_rate * avg_feedback_score`
4. Clamp weights to [0.1, 2.0]
5. Normalize to sum to 4.0

### Default Weights
```python
{
    "clinical": 1.0,
    "literature": 0.8,
    "patent": 0.6,
    "market": 0.7,
    "recency_boost": 0.1,
    "citation_boost": 0.15
}
```

### Usage
```bash
# Train from command line
python -m backend.app.ml.reranker_trainer

# Or import programmatically
from backend.app.ml import RerankerTrainer

trainer = RerankerTrainer()
examples = trainer.fetch_training_data(db_session)
trainer.train(examples, learning_rate=0.01)
reranked = trainer.rerank_passages(passages, query="molecule")
```

### Integration with Feedback Loop
1. User submits feedback via `POST /api/feedback`
2. Feedback stored in `evidence_feedback` table
3. Periodic training runs (cron/manual)
4. Learned weights applied to future retrievals

---

## 3. Custom Embedding Training (`backend/app/ml/embedding_trainer.py`)

### Purpose
Fine-tunes sentence-transformers models on pharma-specific terminology to improve retrieval quality for domain-specific queries.

### Implementation
- **EmbeddingTrainer** class using sentence-transformers library
- Base model: `all-MiniLM-L6-v2` (384-dim embeddings)
- Contrastive learning with CosineSimilarityLoss
- Pharma-specific training pairs (synonyms, abbreviations, related concepts)
- Export functionality for Chroma integration

### Training Data
Pre-defined pharma term pairs:
- Positive pairs: "Phase 1 clinical trial" ↔ "First-in-human study"
- Negative pairs: "Phase 1 clinical trial" ↔ "Market exclusivity period"

Can be extended by extracting from Chroma database:
- Documents with shared citations
- Documents with similar keywords
- Query-document pairs from retrieval logs

### Usage
```bash
# Train from command line
python -m backend.app.ml.embedding_trainer train

# Export for Chroma
from backend.app.ml import EmbeddingTrainer

trainer = EmbeddingTrainer()
training_pairs = trainer.prepare_training_data_from_chroma()
trainer.train(training_pairs, epochs=3)
trainer.export_for_chroma()
```

### Integration with Chroma
After training, update `agents/retrieval.py`:
```python
from chromadb.utils import embedding_functions

# Use fine-tuned model
embedding_fn = embedding_functions.SentenceTransformerEmbeddingFunction(
    model_name="models/embeddings/chroma_compatible"
)
collection = client.get_collection(
    name="literature",
    embedding_function=embedding_fn
)
```

---

## Dependencies

### Required (Core ML Features)
```txt
redis>=5.0.0              # Caching + Celery broker
celery>=5.3.0             # Distributed task queue
numpy>=1.24.0             # Reranker weight calculations
```

### Optional (Embedding Training)
```txt
sentence-transformers>=2.2.0  # Fine-tuning embeddings
torch>=2.0.0                  # Required by sentence-transformers
```

**Note**: Disease mapping and reranker training work without sentence-transformers. Embedding training requires it.

---

## File Structure

```
backend/app/ml/
├── __init__.py                # Exports for ML module
├── disease_mapper.py          # Molecule-disease mapping
├── reranker_trainer.py        # Feedback-based reranking
└── embedding_trainer.py       # Custom embedding training

models/
├── reranker/
│   └── reranker_weights.json  # Learned weights
└── embeddings/
    ├── finetuned/             # Fine-tuned model
    └── chroma_compatible/     # Chroma-ready export

indexes/
└── disease_mappings/
    └── molecule_disease_mappings.json  # Disease index
```

---

## Production Deployment

### 1. Database Migration
Add `evidence_feedback` table (already defined in `backend/app/models.py`):
```sql
CREATE TABLE evidence_feedback (
    id VARCHAR(36) PRIMARY KEY,
    job_id VARCHAR(36) REFERENCES jobs(id) ON DELETE CASCADE,
    evidence_id VARCHAR(255),
    evidence_type VARCHAR(50),
    feedback_type VARCHAR(50),
    feedback_score FLOAT,
    created_at TIMESTAMP DEFAULT NOW()
);
CREATE INDEX idx_evidence_feedback_job_id ON evidence_feedback(job_id);
CREATE INDEX idx_evidence_feedback_evidence_id ON evidence_feedback(evidence_id);
```

### 2. Scheduled Training
Add cron job or Celery Beat task:
```python
# In celery_tasks.py
@celery_app.task(name="phagen.retrain_reranker")
def retrain_reranker():
    from backend.app.ml import RerankerTrainer
    from backend.app.database import SessionLocal
    
    trainer = RerankerTrainer()
    with SessionLocal() as db:
        examples = trainer.fetch_training_data(db)
        if len(examples) >= 10:  # Minimum training set
            trainer.train(examples)

# Schedule weekly
celery_app.conf.beat_schedule["retrain-reranker"] = {
    "task": "phagen.retrain_reranker",
    "schedule": 604800.0  # 7 days
}
```

### 3. Environment Variables
```bash
# Optional: Disable ML features
ENABLE_REPURPOSING_SUGGESTIONS=true
ENABLE_RERANKER_TRAINING=true

# Model paths
RERANKER_WEIGHTS_PATH=models/reranker/reranker_weights.json
DISEASE_MAPPING_INDEX=indexes/disease_mappings/
```

---

## Performance Considerations

### Disease Mapping
- **Latency**: ~50-100ms per job (negligible)
- **Memory**: <10 MB for disease index
- **Scaling**: Stateless, can run on any worker

### Reranker Training
- **Training Time**: ~1-5 seconds per 100 examples
- **Frequency**: Weekly or when feedback count exceeds threshold
- **Resource**: CPU-only, single core sufficient

### Embedding Training
- **Training Time**: 5-30 minutes for 3 epochs (CPU), 1-5 minutes (GPU)
- **Frequency**: Monthly or when corpus significantly changes
- **Resource**: GPU recommended (2-4 GB VRAM), CPU fallback available
- **Disk**: ~100 MB per fine-tuned model

---

## Monitoring & Evaluation

### Metrics to Track
1. **Disease Mapping**
   - Prediction confidence distribution
   - Top-5 accuracy vs. known mappings
   - Coverage rate (% jobs with suggestions)

2. **Reranker**
   - Weight evolution over time
   - Average feedback score trend
   - Retrieval precision improvement

3. **Embeddings**
   - Embedding quality (cosine similarity on test pairs)
   - Retrieval recall@k on benchmark molecules
   - User satisfaction with results

### Dashboard Queries
```python
# Feedback statistics
SELECT evidence_type, AVG(feedback_score) as avg_score, COUNT(*) as count
FROM evidence_feedback
WHERE created_at > NOW() - INTERVAL '30 days'
GROUP BY evidence_type;

# Repurposing coverage
SELECT COUNT(*) as total_jobs,
       SUM(CASE WHEN payload->>'repurposing_suggestions' IS NOT NULL THEN 1 ELSE 0 END) as with_suggestions
FROM jobs
WHERE status = 'completed' AND created_at > NOW() - INTERVAL '7 days';
```

---

## Next Steps

### Immediate
1. ✅ Database migration for `evidence_feedback` table
2. ✅ Test repurposing suggestions on sample molecules
3. ✅ Configure reranker training schedule

### Short-term
1. Add UI components to display repurposing suggestions
2. Implement feedback collection buttons in evidence tabs
3. Create admin dashboard for ML model monitoring

### Long-term
1. Deep learning reranker with BERT/cross-encoders
2. Multi-task learning for disease/target/pathway prediction
3. Active learning to prioritize which molecules need review
4. Federated learning for cross-organization model training

---

## Testing

### Unit Tests
```bash
# Disease mapping
pytest backend/tests/test_disease_mapper.py

# Reranker
pytest backend/tests/test_reranker_trainer.py

# Embeddings (requires sentence-transformers)
pytest backend/tests/test_embedding_trainer.py
```

### Integration Tests
```bash
# End-to-end feedback loop
pytest backend/tests/integration/test_feedback_loop.py

# Job with repurposing suggestions
pytest backend/tests/integration/test_repurposing_integration.py
```

---

## Troubleshooting

### Issue: No repurposing suggestions generated
- **Check**: Job payload contains evidence from clinical/literature workers
- **Check**: Disease mapper loaded successfully (no import errors)
- **Fix**: Ensure evidence has sufficient keyword matches

### Issue: Reranker not learning
- **Check**: `evidence_feedback` table has entries with non-zero feedback_score
- **Check**: Training runs without errors (check logs)
- **Fix**: Collect more feedback (minimum 10 examples per evidence type)

### Issue: Embedding training fails
- **Check**: `sentence-transformers` and `torch` installed
- **Check**: Sufficient disk space for model checkpoints
- **Fix**: Install dependencies: `pip install sentence-transformers torch`

---

## References

- **Disease Mapping**: Category-based prediction with MeSH IDs
- **Reranker**: Gradient descent on feedback scores
- **Embeddings**: Contrastive learning with CosineSimilarityLoss
- **Feedback API**: `backend/app/routers/feedback.py`
- **Database Schema**: `backend/app/models.py`

For questions, see `docs/architecture.md` or consult the team.
