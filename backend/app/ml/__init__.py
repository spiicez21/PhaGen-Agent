"""
Machine learning models for PhaGen.
"""
from .reranker_trainer import RerankerTrainer, train_reranker_from_feedback
from .disease_mapper import MoleculeDiseaseMapper, generate_repurposing_suggestions

try:
    from .embedding_trainer import EmbeddingTrainer, train_embeddings_from_corpus
except ImportError:
    # sentence-transformers optional
    EmbeddingTrainer = None
    train_embeddings_from_corpus = None

__all__ = [
    "RerankerTrainer",
    "train_reranker_from_feedback",
    "MoleculeDiseaseMapper",
    "generate_repurposing_suggestions",
    "EmbeddingTrainer",
    "train_embeddings_from_corpus",
]
