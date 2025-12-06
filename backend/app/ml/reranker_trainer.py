"""
Reranker fine-tuning pipeline using historical feedback for retrieval optimization.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class TrainingExample:
    """Training example for reranker fine-tuning."""
    query: str
    passage: str
    label: float  # 1.0 for upvote, -1.0 for downvote, 0.0 for flag
    evidence_type: str
    metadata: Dict


class RerankerTrainer:
    """
    Fine-tune reranker weights based on user feedback.
    Uses historical upvote/downvote data to learn better relevance scoring.
    """
    
    def __init__(self, model_path: Path = Path("models/reranker")):
        self.model_path = model_path
        self.model_path.mkdir(parents=True, exist_ok=True)
        self.weights_file = self.model_path / "reranker_weights.json"
        self.weights = self._load_weights()
    
    def _load_weights(self) -> Dict[str, float]:
        """Load existing weights or initialize defaults."""
        if self.weights_file.exists():
            with open(self.weights_file) as f:
                return json.load(f)
        
        # Default weights favor clinical evidence
        return {
            "clinical": 1.0,
            "literature": 0.8,
            "patent": 0.6,
            "market": 0.7,
            "recency_boost": 0.1,
            "citation_boost": 0.15,
        }
    
    def _save_weights(self):
        """Persist weights to disk."""
        with open(self.weights_file, 'w') as f:
            json.dump(self.weights, f, indent=2)
        logger.info(f"Saved reranker weights to {self.weights_file}")
    
    def fetch_training_data(self, db_session) -> List[TrainingExample]:
        """
        Fetch feedback data from database.
        
        Args:
            db_session: SQLAlchemy session
        
        Returns:
            List of training examples
        """
        from backend.app.models import EvidenceFeedback, Job
        
        examples = []
        
        # Query feedback with job context
        feedback_list = db_session.query(EvidenceFeedback, Job).join(
            Job, EvidenceFeedback.job_id == Job.id
        ).filter(
            EvidenceFeedback.feedback_score != 0.0  # Exclude flags
        ).all()
        
        for feedback, job in feedback_list:
            # Extract query from job payload
            query = job.payload.get("molecule", "") if job.payload else ""
            
            # Use evidence_id as passage identifier
            passage = feedback.evidence_id
            
            examples.append(TrainingExample(
                query=query,
                passage=passage,
                label=feedback.feedback_score,
                evidence_type=feedback.evidence_type,
                metadata={"job_id": job.id, "feedback_id": feedback.id}
            ))
        
        logger.info(f"Fetched {len(examples)} training examples")
        return examples
    
    def train(self, examples: List[TrainingExample], learning_rate: float = 0.01):
        """
        Train reranker using gradient descent on feedback scores.
        
        Args:
            examples: Training examples from user feedback
            learning_rate: Learning rate for weight updates
        """
        if not examples:
            logger.warning("No training examples, skipping training")
            return
        
        # Group by evidence type
        type_scores = {t: [] for t in ["clinical", "literature", "patent", "market"]}
        
        for ex in examples:
            if ex.evidence_type in type_scores:
                type_scores[ex.evidence_type].append(ex.label)
        
        # Update weights based on average feedback
        for ev_type, scores in type_scores.items():
            if scores:
                avg_score = np.mean(scores)
                # Adjust weight: positive feedback increases weight, negative decreases
                delta = learning_rate * avg_score
                self.weights[ev_type] = max(0.1, min(2.0, self.weights[ev_type] + delta))
        
        # Normalize weights
        total = sum(self.weights[k] for k in ["clinical", "literature", "patent", "market"])
        for k in ["clinical", "literature", "patent", "market"]:
            self.weights[k] = (self.weights[k] / total) * 4.0
        
        self._save_weights()
        logger.info(f"Updated reranker weights: {self.weights}")
    
    def rerank_passages(self, passages: List[Dict], query: str) -> List[Dict]:
        """
        Rerank passages using learned weights.
        
        Args:
            passages: List of passage dicts with 'type' and 'score' fields
            query: Query string (unused in this simple version)
        
        Returns:
            Reranked list of passages
        """
        scored = []
        for p in passages:
            ev_type = p.get("type", "literature")
            base_score = p.get("score", 0.5)
            
            # Apply type weight
            type_weight = self.weights.get(ev_type, 0.5)
            
            # Apply recency boost if available
            recency_boost = 0.0
            if "year" in p:
                years_old = 2025 - p["year"]
                if years_old <= 3:
                    recency_boost = self.weights["recency_boost"]
            
            # Apply citation boost if available
            citation_boost = 0.0
            if "citations" in p and p["citations"] > 10:
                citation_boost = self.weights["citation_boost"]
            
            final_score = base_score * type_weight + recency_boost + citation_boost
            
            scored.append((final_score, p))
        
        # Sort by score descending
        scored.sort(key=lambda x: x[0], reverse=True)
        
        return [p for _, p in scored]
    
    def get_stats(self) -> Dict:
        """Get current reranker statistics."""
        return {
            "weights": self.weights,
            "model_path": str(self.model_path),
            "weights_file": str(self.weights_file),
        }


def train_reranker_from_feedback():
    """
    CLI tool to train reranker from database feedback.
    Run: python -m backend.app.ml.reranker_trainer
    """
    from backend.app.database import SessionLocal
    
    trainer = RerankerTrainer()
    db = SessionLocal()
    
    try:
        examples = trainer.fetch_training_data(db)
        
        if not examples:
            logger.info("No feedback data available for training")
            return
        
        logger.info(f"Training on {len(examples)} examples")
        trainer.train(examples, learning_rate=0.05)
        
        stats = trainer.get_stats()
        logger.info(f"Training complete. Stats: {stats}")
        
    finally:
        db.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    train_reranker_from_feedback()
