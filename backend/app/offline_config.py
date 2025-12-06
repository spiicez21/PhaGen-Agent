"""
PhaGen Offline Mode Configuration
Enables operation in air-gapped environments with no external network access
"""

import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class OfflineConfig:
    """Configuration for offline/air-gapped operation"""
    
    def __init__(self):
        # Check if offline mode is enabled
        self.enabled = os.getenv('PHAGEN_OFFLINE_MODE', 'false').lower() == 'true'
        
        # Paths to local resources
        self.models_dir = Path(os.getenv('PHAGEN_MODELS_DIR', '/opt/phagen/models'))
        self.indexes_dir = Path(os.getenv('PHAGEN_INDEXES_DIR', '/opt/phagen/indexes'))
        
        # Feature flags for offline mode
        self.disable_external_apis = self.enabled
        self.disable_telemetry = self.enabled
        self.disable_updates = self.enabled
        
        if self.enabled:
            logger.info("🔒 PhaGen running in OFFLINE MODE")
            logger.info(f"   Models: {self.models_dir}")
            logger.info(f"   Indexes: {self.indexes_dir}")
    
    @property
    def is_online(self) -> bool:
        """Check if external network access is allowed"""
        return not self.enabled
    
    def get_model_path(self, model_name: str) -> Optional[Path]:
        """
        Get local path for a model
        
        Args:
            model_name: Model identifier (e.g., 'sentence-transformers/all-MiniLM-L6-v2')
        
        Returns:
            Path to local model or None if not found
        """
        # Convert model name to directory name
        model_dir_name = model_name.replace('/', '_')
        model_path = self.models_dir / model_dir_name
        
        if model_path.exists():
            return model_path
        
        logger.warning(f"Model not found in offline bundle: {model_name}")
        return None
    
    def get_index_path(self, index_name: str) -> Optional[Path]:
        """
        Get local path for an index
        
        Args:
            index_name: Index identifier (e.g., 'chroma')
        
        Returns:
            Path to local index or None if not found
        """
        index_path = self.indexes_dir / index_name
        
        if index_path.exists():
            return index_path
        
        logger.warning(f"Index not found in offline bundle: {index_name}")
        return None


# Global offline config instance
offline_config = OfflineConfig()


def is_offline_mode() -> bool:
    """Check if PhaGen is running in offline mode"""
    return offline_config.enabled


def require_online(feature_name: str):
    """
    Decorator to mark functions that require online access
    
    Args:
        feature_name: Name of the feature for error message
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            if is_offline_mode():
                raise RuntimeError(
                    f"{feature_name} is not available in offline mode. "
                    f"Set PHAGEN_OFFLINE_MODE=false to enable."
                )
            return func(*args, **kwargs)
        return wrapper
    return decorator


# Example usage in code:
#
# @require_online("PubMed search")
# def search_pubmed(query: str):
#     # This will raise an error in offline mode
#     ...
#
# if is_offline_mode():
#     # Use local models
#     model_path = offline_config.get_model_path("sentence-transformers/all-MiniLM-L6-v2")
#     model = SentenceTransformer(str(model_path))
# else:
#     # Download from internet
#     model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")
