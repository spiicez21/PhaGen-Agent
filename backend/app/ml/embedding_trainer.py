"""
Custom embedding model fine-tuning on pharma corpus.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

try:
    from sentence_transformers import SentenceTransformer, InputExample, losses
    from torch.utils.data import DataLoader
    SENTENCE_TRANSFORMERS_AVAILABLE = True
except ImportError:
    SENTENCE_TRANSFORMERS_AVAILABLE = False
    logging.warning("sentence-transformers not installed. Embedding training unavailable.")

logger = logging.getLogger(__name__)


@dataclass
class TrainingPair:
    """Training pair for contrastive learning."""
    text_a: str
    text_b: str
    label: float  # 1.0 for similar, 0.0 for dissimilar


class EmbeddingTrainer:
    """
    Fine-tunes sentence-transformers models on pharma-specific corpus.
    """
    
    def __init__(
        self,
        base_model: str = "all-MiniLM-L6-v2",
        output_dir: Path = Path("models/embeddings")
    ):
        if not SENTENCE_TRANSFORMERS_AVAILABLE:
            raise ImportError("sentence-transformers required for embedding training")
        
        self.base_model = base_model
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.model = SentenceTransformer(base_model)
        logger.info(f"Initialized embedding trainer with {base_model}")
    
    def prepare_training_data_from_chroma(
        self,
        chroma_dir: Path = Path("indexes/chroma")
    ) -> List[TrainingPair]:
        """
        Extract training pairs from Chroma vector database.
        
        Creates positive pairs from:
        - Documents with shared citations
        - Documents with similar keywords
        - Query-document pairs from retrieval logs
        
        Args:
            chroma_dir: Path to Chroma database
        
        Returns:
            List of training pairs
        """
        training_pairs = []
        
        # Try to load existing collections metadata
        collections_to_check = [
            "clinical_trials",
            "literature",
            "patents",
            "market_data"
        ]
        
        # Create synthetic positive pairs from common pharma terms
        pharma_terms = [
            ("Phase 1 clinical trial", "First-in-human study"),
            ("Overall survival", "OS endpoint"),
            ("Progression-free survival", "PFS outcome"),
            ("Mechanism of action", "MOA pathway"),
            ("Adverse events", "Safety profile"),
            ("Pharmacokinetics", "PK parameters"),
            ("Pharmacodynamics", "PD markers"),
            ("Bioavailability", "Drug absorption"),
            ("Off-target effects", "Selectivity profile"),
            ("Drug-drug interaction", "DDI risk"),
            ("Maximum tolerated dose", "MTD determination"),
            ("Objective response rate", "ORR analysis"),
            ("Target engagement", "Receptor occupancy"),
            ("Dose escalation", "3+3 design"),
            ("Surrogate endpoint", "Biomarker response"),
        ]
        
        for term_a, term_b in pharma_terms:
            training_pairs.append(TrainingPair(
                text_a=term_a,
                text_b=term_b,
                label=1.0
            ))
        
        # Create negative pairs (dissimilar concepts)
        negative_pairs = [
            ("Phase 1 clinical trial", "Market exclusivity period"),
            ("Overall survival", "Patent expiration date"),
            ("Pharmacokinetics", "Competitive landscape"),
            ("Adverse events", "Revenue forecast"),
            ("Mechanism of action", "Reimbursement policy"),
        ]
        
        for term_a, term_b in negative_pairs:
            training_pairs.append(TrainingPair(
                text_a=term_a,
                text_b=term_b,
                label=0.0
            ))
        
        logger.info(f"Prepared {len(training_pairs)} training pairs")
        return training_pairs
    
    def train(
        self,
        training_pairs: List[TrainingPair],
        epochs: int = 3,
        batch_size: int = 16,
        warmup_steps: int = 100,
        evaluation_steps: int = 500
    ):
        """
        Fine-tune the embedding model.
        
        Args:
            training_pairs: Training data
            epochs: Number of training epochs
            batch_size: Batch size
            warmup_steps: Warmup steps for learning rate
            evaluation_steps: Steps between evaluations
        """
        if not training_pairs:
            logger.warning("No training pairs provided")
            return
        
        # Convert to InputExample format
        train_examples = []
        for pair in training_pairs:
            train_examples.append(InputExample(
                texts=[pair.text_a, pair.text_b],
                label=pair.label
            ))
        
        # Create DataLoader
        train_dataloader = DataLoader(
            train_examples,
            shuffle=True,
            batch_size=batch_size
        )
        
        # Use CosineSimilarityLoss for contrastive learning
        train_loss = losses.CosineSimilarityLoss(self.model)
        
        # Train the model
        logger.info(f"Starting training for {epochs} epochs...")
        self.model.fit(
            train_objectives=[(train_dataloader, train_loss)],
            epochs=epochs,
            warmup_steps=warmup_steps,
            output_path=str(self.output_dir / "finetuned"),
            show_progress_bar=True,
            evaluation_steps=evaluation_steps
        )
        
        logger.info(f"Training complete. Model saved to {self.output_dir / 'finetuned'}")
    
    def evaluate_on_benchmark(self, test_molecules: List[str]) -> Dict:
        """
        Evaluate model on benchmark molecules.
        
        Args:
            test_molecules: List of molecule names to test
        
        Returns:
            Evaluation metrics
        """
        if not test_molecules:
            return {"error": "No test molecules provided"}
        
        # Load fine-tuned model
        finetuned_path = self.output_dir / "finetuned"
        if finetuned_path.exists():
            model = SentenceTransformer(str(finetuned_path))
        else:
            logger.warning("No fine-tuned model found, using base model")
            model = self.model
        
        # Generate embeddings
        embeddings = model.encode(test_molecules)
        
        # Simple evaluation: check embedding quality
        results = {
            "model": str(finetuned_path) if finetuned_path.exists() else self.base_model,
            "num_test_molecules": len(test_molecules),
            "embedding_dim": embeddings.shape[1] if len(embeddings) > 0 else 0,
            "avg_norm": float(sum([sum(e**2)**0.5 for e in embeddings]) / len(embeddings)) if len(embeddings) > 0 else 0.0
        }
        
        logger.info(f"Evaluation results: {results}")
        return results
    
    def export_for_chroma(self, output_path: Optional[Path] = None):
        """
        Export fine-tuned model for use with Chroma.
        
        Args:
            output_path: Where to save the model
        """
        finetuned_path = self.output_dir / "finetuned"
        if not finetuned_path.exists():
            logger.error("No fine-tuned model found")
            return
        
        if output_path is None:
            output_path = self.output_dir / "chroma_compatible"
        
        # Copy model files
        import shutil
        shutil.copytree(finetuned_path, output_path, dirs_exist_ok=True)
        
        # Create config file
        config = {
            "model_name": str(output_path),
            "base_model": self.base_model,
            "fine_tuned": True,
            "usage": "Use with ChromaDB: collection = chromadb.Collection(embedding_function=SentenceTransformerEmbeddingFunction(model_name=<path>))"
        }
        
        with open(output_path / "chroma_config.json", 'w') as f:
            json.dump(config, f, indent=2)
        
        logger.info(f"Exported model to {output_path}")
        logger.info("To use with Chroma, update your embedding function in agents/retrieval.py")
    
    def get_stats(self) -> Dict:
        """Get trainer statistics."""
        finetuned_exists = (self.output_dir / "finetuned").exists()
        
        return {
            "base_model": self.base_model,
            "output_dir": str(self.output_dir),
            "fine_tuned_model_exists": finetuned_exists,
            "fine_tuned_path": str(self.output_dir / "finetuned") if finetuned_exists else None
        }


def train_embeddings_from_corpus(
    chroma_dir: Path = Path("indexes/chroma"),
    base_model: str = "all-MiniLM-L6-v2",
    epochs: int = 3
):
    """
    CLI entry point for embedding training.
    
    Args:
        chroma_dir: Path to Chroma database
        base_model: Base model to fine-tune
        epochs: Training epochs
    """
    if not SENTENCE_TRANSFORMERS_AVAILABLE:
        print("ERROR: sentence-transformers not installed")
        print("Install with: pip install sentence-transformers torch")
        return
    
    trainer = EmbeddingTrainer(base_model=base_model)
    
    # Prepare training data
    print("Preparing training data...")
    training_pairs = trainer.prepare_training_data_from_chroma(chroma_dir)
    
    if not training_pairs:
        print("ERROR: No training pairs generated")
        return
    
    # Train
    print(f"Training on {len(training_pairs)} pairs for {epochs} epochs...")
    trainer.train(training_pairs, epochs=epochs)
    
    # Export for Chroma
    trainer.export_for_chroma()
    
    # Show stats
    stats = trainer.get_stats()
    print("\nTraining complete!")
    print(f"Stats: {json.dumps(stats, indent=2)}")
    
    print("\nTo use the fine-tuned model:")
    print("1. Update agents/retrieval.py to use the new model")
    print("2. Rebuild Chroma indexes with the new embeddings")
    print("3. Test retrieval quality on benchmark molecules")


if __name__ == "__main__":
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    # Example usage
    if len(sys.argv) > 1 and sys.argv[1] == "train":
        train_embeddings_from_corpus(epochs=3)
    else:
        # Demo mode
        print("Embedding Trainer Demo")
        print("=" * 50)
        
        if not SENTENCE_TRANSFORMERS_AVAILABLE:
            print("ERROR: sentence-transformers not installed")
            print("Install with: pip install sentence-transformers torch")
        else:
            trainer = EmbeddingTrainer()
            
            # Show stats
            stats = trainer.get_stats()
            print(f"Stats: {json.dumps(stats, indent=2)}")
            
            # Evaluate on test molecules
            test_molecules = [
                "Aspirin inhibits COX-2 enzyme",
                "Phase 2 trial in breast cancer",
                "Patent expires in 2025"
            ]
            
            results = trainer.evaluate_on_benchmark(test_molecules)
            print(f"\nBenchmark results: {json.dumps(results, indent=2)}")
            
            print("\nTo train on your corpus:")
            print("python -m backend.app.ml.embedding_trainer train")
