"""
Continuous crawler with change detection for automatic data refresh.
"""
from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, asdict

logger = logging.getLogger(__name__)


@dataclass
class DocumentSnapshot:
    """Snapshot of a crawled document."""
    url: str
    content_hash: str
    last_crawled: str  # ISO format
    last_modified: Optional[str] = None
    etag: Optional[str] = None
    freshness_score: float = 1.0


@dataclass
class ChangeDetectionResult:
    """Result of change detection analysis."""
    url: str
    status: str  # 'unchanged', 'modified', 'new', 'deleted'
    content_hash_old: Optional[str] = None
    content_hash_new: Optional[str] = None
    last_check: Optional[str] = None
    needs_reindex: bool = False


class ContinuousCrawler:
    """
    Monitors sources for changes and triggers incremental updates.
    """
    
    def __init__(
        self,
        snapshot_file: Path = Path("crawler/storage/.document_snapshots.json"),
        freshness_threshold_days: int = 30
    ):
        self.snapshot_file = snapshot_file
        self.snapshot_file.parent.mkdir(parents=True, exist_ok=True)
        self.freshness_threshold = freshness_threshold_days
        self.snapshots: Dict[str, DocumentSnapshot] = {}
        self._load_snapshots()
    
    def _load_snapshots(self):
        """Load existing document snapshots."""
        if self.snapshot_file.exists():
            try:
                with open(self.snapshot_file) as f:
                    data = json.load(f)
                    self.snapshots = {
                        url: DocumentSnapshot(**snap)
                        for url, snap in data.items()
                    }
                logger.info(f"Loaded {len(self.snapshots)} document snapshots")
            except Exception as e:
                logger.error(f"Failed to load snapshots: {e}")
                self.snapshots = {}
        else:
            logger.info("No existing snapshots found")
    
    def _save_snapshots(self):
        """Persist document snapshots."""
        try:
            data = {
                url: asdict(snap)
                for url, snap in self.snapshots.items()
            }
            with open(self.snapshot_file, 'w') as f:
                json.dump(data, f, indent=2)
            logger.info(f"Saved {len(self.snapshots)} snapshots")
        except Exception as e:
            logger.error(f"Failed to save snapshots: {e}")
    
    def compute_content_hash(self, content: str) -> str:
        """
        Compute SHA-256 hash of content.
        
        Args:
            content: Document content
        
        Returns:
            Hex digest of content hash
        """
        return hashlib.sha256(content.encode('utf-8')).hexdigest()
    
    def detect_changes(
        self,
        current_documents: List[Dict]
    ) -> List[ChangeDetectionResult]:
        """
        Detect changes in crawled documents.
        
        Args:
            current_documents: List of currently crawled documents
        
        Returns:
            List of change detection results
        """
        results = []
        current_urls = set()
        
        for doc in current_documents:
            url = doc.get('url', '')
            if not url:
                continue
            
            current_urls.add(url)
            content = doc.get('text', '') or doc.get('content', '')
            current_hash = self.compute_content_hash(content)
            
            # Check if document exists in snapshots
            if url in self.snapshots:
                old_snap = self.snapshots[url]
                
                if old_snap.content_hash == current_hash:
                    # Unchanged
                    results.append(ChangeDetectionResult(
                        url=url,
                        status='unchanged',
                        content_hash_old=old_snap.content_hash,
                        content_hash_new=current_hash,
                        last_check=datetime.utcnow().isoformat(),
                        needs_reindex=False
                    ))
                else:
                    # Modified
                    results.append(ChangeDetectionResult(
                        url=url,
                        status='modified',
                        content_hash_old=old_snap.content_hash,
                        content_hash_new=current_hash,
                        last_check=datetime.utcnow().isoformat(),
                        needs_reindex=True
                    ))
                    
                    # Update snapshot
                    self.snapshots[url] = DocumentSnapshot(
                        url=url,
                        content_hash=current_hash,
                        last_crawled=datetime.utcnow().isoformat(),
                        last_modified=doc.get('last_modified'),
                        etag=doc.get('etag')
                    )
            else:
                # New document
                results.append(ChangeDetectionResult(
                    url=url,
                    status='new',
                    content_hash_new=current_hash,
                    last_check=datetime.utcnow().isoformat(),
                    needs_reindex=True
                ))
                
                # Add to snapshots
                self.snapshots[url] = DocumentSnapshot(
                    url=url,
                    content_hash=current_hash,
                    last_crawled=datetime.utcnow().isoformat(),
                    last_modified=doc.get('last_modified'),
                    etag=doc.get('etag')
                )
        
        # Check for deleted documents
        deleted_urls = set(self.snapshots.keys()) - current_urls
        for url in deleted_urls:
            results.append(ChangeDetectionResult(
                url=url,
                status='deleted',
                content_hash_old=self.snapshots[url].content_hash,
                last_check=datetime.utcnow().isoformat(),
                needs_reindex=True
            ))
            # Remove from snapshots
            del self.snapshots[url]
        
        # Save updated snapshots
        self._save_snapshots()
        
        return results
    
    def compute_freshness_scores(self) -> Dict[str, float]:
        """
        Compute freshness scores for all documents.
        
        Returns:
            Mapping of URL to freshness score (0.0 to 1.0)
        """
        freshness_scores = {}
        now = datetime.utcnow()
        threshold = timedelta(days=self.freshness_threshold)
        
        for url, snap in self.snapshots.items():
            try:
                last_crawled = datetime.fromisoformat(snap.last_crawled)
                age = now - last_crawled
                
                # Linear decay over threshold period
                if age > threshold:
                    freshness = max(0.0, 1.0 - (age.days - threshold.days) / threshold.days)
                else:
                    freshness = 1.0
                
                snap.freshness_score = round(freshness, 3)
                freshness_scores[url] = snap.freshness_score
            except Exception as e:
                logger.warning(f"Failed to compute freshness for {url}: {e}")
                freshness_scores[url] = 0.5
        
        return freshness_scores
    
    def get_stale_documents(
        self,
        min_freshness: float = 0.5
    ) -> List[DocumentSnapshot]:
        """
        Get documents that need re-crawling.
        
        Args:
            min_freshness: Minimum freshness threshold
        
        Returns:
            List of stale document snapshots
        """
        self.compute_freshness_scores()
        
        stale = [
            snap for snap in self.snapshots.values()
            if snap.freshness_score < min_freshness
        ]
        
        logger.info(f"Found {len(stale)} stale documents (freshness < {min_freshness})")
        return stale
    
    def generate_change_report(
        self,
        changes: List[ChangeDetectionResult]
    ) -> Dict:
        """
        Generate summary report of changes.
        
        Args:
            changes: List of change detection results
        
        Returns:
            Change summary report
        """
        new_count = sum(1 for c in changes if c.status == 'new')
        modified_count = sum(1 for c in changes if c.status == 'modified')
        deleted_count = sum(1 for c in changes if c.status == 'deleted')
        unchanged_count = sum(1 for c in changes if c.status == 'unchanged')
        needs_reindex = [c for c in changes if c.needs_reindex]
        
        return {
            "timestamp": datetime.utcnow().isoformat(),
            "total_documents": len(changes),
            "new_documents": new_count,
            "modified_documents": modified_count,
            "deleted_documents": deleted_count,
            "unchanged_documents": unchanged_count,
            "needs_reindex_count": len(needs_reindex),
            "needs_reindex_urls": [c.url for c in needs_reindex[:10]],  # Top 10
            "recommendation": self._generate_recommendation(
                new_count, modified_count, deleted_count
            )
        }
    
    def _generate_recommendation(
        self,
        new_count: int,
        modified_count: int,
        deleted_count: int
    ) -> str:
        """Generate recommendation based on change counts."""
        total_changes = new_count + modified_count + deleted_count
        
        if total_changes == 0:
            return "No changes detected. Index is up to date."
        elif total_changes < 10:
            return f"Minor changes detected ({total_changes} documents). Consider incremental reindex."
        elif total_changes < 50:
            return f"Moderate changes detected ({total_changes} documents). Schedule incremental reindex."
        else:
            return f"Major changes detected ({total_changes} documents). Full reindex recommended."
    
    def get_stats(self) -> Dict:
        """Get crawler statistics."""
        freshness_scores = self.compute_freshness_scores()
        
        return {
            "total_snapshots": len(self.snapshots),
            "snapshot_file": str(self.snapshot_file),
            "avg_freshness": round(
                sum(freshness_scores.values()) / len(freshness_scores)
                if freshness_scores else 0.0,
                3
            ),
            "freshness_threshold_days": self.freshness_threshold,
            "stale_documents": len(self.get_stale_documents(0.5))
        }


def run_continuous_crawl_check(
    crawler_output_dir: Path = Path("crawler/storage/datasets/default"),
    snapshot_file: Path = Path("crawler/storage/.document_snapshots.json")
) -> Dict:
    """
    Run continuous crawl check on existing crawler output.
    
    Args:
        crawler_output_dir: Directory with crawler output
        snapshot_file: Path to snapshot file
    
    Returns:
        Change detection report
    """
    crawler = ContinuousCrawler(snapshot_file=snapshot_file)
    
    # Load current documents
    current_documents = []
    if crawler_output_dir.exists():
        for json_file in crawler_output_dir.glob("*.json"):
            try:
                with open(json_file) as f:
                    doc = json.load(f)
                    current_documents.append(doc)
            except Exception as e:
                logger.warning(f"Failed to load {json_file}: {e}")
    
    logger.info(f"Loaded {len(current_documents)} current documents")
    
    # Detect changes
    changes = crawler.detect_changes(current_documents)
    
    # Generate report
    report = crawler.generate_change_report(changes)
    
    # Add stats
    report["crawler_stats"] = crawler.get_stats()
    
    return report


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Run change detection
    report = run_continuous_crawl_check()
    
    print("=" * 60)
    print("Continuous Crawler Change Detection Report")
    print("=" * 60)
    print(json.dumps(report, indent=2))
    
    print("\nTo integrate with your workflow:")
    print("1. Run after each crawl: python -m agents.workers.continuous_crawler")
    print("2. Schedule periodic checks (e.g., daily cron job)")
    print("3. Trigger index rebuild when needs_reindex_count > threshold")
