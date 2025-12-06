#!/usr/bin/env python
"""
CLI helper for training ML models.
"""
import sys
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def train_reranker(args):
    """Train reranker from feedback."""
    print("=" * 60)
    print("Training Reranker from Feedback")
    print("=" * 60)
    
    from backend.app.ml.reranker_trainer import train_reranker_from_feedback
    train_reranker_from_feedback()

def train_embeddings(args):
    """Train custom embeddings."""
    print("=" * 60)
    print("Training Custom Embeddings")
    print("=" * 60)
    
    try:
        from backend.app.ml.embedding_trainer import train_embeddings_from_corpus
    except ImportError:
        print("ERROR: sentence-transformers not installed")
        print("Install with: pip install sentence-transformers torch")
        sys.exit(1)
    
    from pathlib import Path
    
    train_embeddings_from_corpus(
        chroma_dir=Path(args.chroma_dir),
        base_model=args.base_model,
        epochs=args.epochs
    )

def show_stats(args):
    """Show ML model statistics."""
    print("=" * 60)
    print("PhaGen ML Model Statistics")
    print("=" * 60)
    
    # Reranker stats
    try:
        from backend.app.ml.reranker_trainer import RerankerTrainer
        reranker = RerankerTrainer()
        stats = reranker.get_stats()
        print("\nReranker:")
        print(f"  Weights file: {stats['weights_file']}")
        print(f"  Weights loaded: {stats['weights_loaded']}")
        if stats['weights_loaded']:
            print(f"  Clinical: {stats['weights'].get('clinical', 0):.2f}")
            print(f"  Literature: {stats['weights'].get('literature', 0):.2f}")
            print(f"  Patent: {stats['weights'].get('patent', 0):.2f}")
            print(f"  Market: {stats['weights'].get('market', 0):.2f}")
    except Exception as e:
        print(f"\nReranker: Error - {e}")
    
    # Disease mapper stats
    try:
        from backend.app.ml.disease_mapper import MoleculeDiseaseMapper
        mapper = MoleculeDiseaseMapper()
        stats = mapper.get_stats()
        print("\nDisease Mapper:")
        print(f"  Total categories: {stats['total_categories']}")
        print(f"  Custom mappings: {stats['custom_mappings']}")
        print(f"  Data dir: {stats['data_dir']}")
    except Exception as e:
        print(f"\nDisease Mapper: Error - {e}")
    
    # Embedding trainer stats
    try:
        from backend.app.ml.embedding_trainer import EmbeddingTrainer, SENTENCE_TRANSFORMERS_AVAILABLE
        
        if not SENTENCE_TRANSFORMERS_AVAILABLE:
            print("\nEmbedding Trainer: Not available (sentence-transformers not installed)")
        else:
            trainer = EmbeddingTrainer()
            stats = trainer.get_stats()
            print("\nEmbedding Trainer:")
            print(f"  Base model: {stats['base_model']}")
            print(f"  Fine-tuned exists: {stats['fine_tuned_model_exists']}")
            if stats['fine_tuned_model_exists']:
                print(f"  Fine-tuned path: {stats['fine_tuned_path']}")
    except Exception as e:
        print(f"\nEmbedding Trainer: Error - {e}")
    
    print("\n" + "=" * 60)

def main():
    parser = argparse.ArgumentParser(
        description="Train and manage PhaGen ML models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show current model stats
  python ml_cli.py stats
  
  # Train reranker from feedback
  python ml_cli.py reranker
  
  # Train custom embeddings
  python ml_cli.py embeddings --epochs 5
  
  # Train with custom base model
  python ml_cli.py embeddings --base-model all-mpnet-base-v2
"""
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Stats command
    subparsers.add_parser('stats', help='Show ML model statistics')
    
    # Reranker command
    subparsers.add_parser('reranker', help='Train reranker from feedback')
    
    # Embeddings command
    embeddings_parser = subparsers.add_parser('embeddings', help='Train custom embeddings')
    embeddings_parser.add_argument(
        '--chroma-dir',
        default='indexes/chroma',
        help='Path to Chroma database (default: indexes/chroma)'
    )
    embeddings_parser.add_argument(
        '--base-model',
        default='all-MiniLM-L6-v2',
        help='Base model to fine-tune (default: all-MiniLM-L6-v2)'
    )
    embeddings_parser.add_argument(
        '--epochs',
        type=int,
        default=3,
        help='Training epochs (default: 3)'
    )
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(0)
    
    if args.command == 'stats':
        show_stats(args)
    elif args.command == 'reranker':
        train_reranker(args)
    elif args.command == 'embeddings':
        train_embeddings(args)

if __name__ == '__main__':
    main()
