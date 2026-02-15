"""
Paid market data API integration framework.
Adapters for IQVIA, Evaluate Pharma, and other commercial data sources.
"""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional
from dataclasses import dataclass
from datetime import datetime, timezone

logger = logging.getLogger(__name__)


@dataclass
class MarketDataPoint:
    """Standardized market data structure."""
    indication: str
    market_size_usd: Optional[float] = None
    cagr_percent: Optional[float] = None
    patient_population: Optional[int] = None
    market_share_leaders: Optional[List[str]] = None
    data_source: str = "synthetic"
    last_updated: Optional[str] = None
    confidence: str = "low"  # low, medium, high


class MarketDataProvider(ABC):
    """Abstract base class for market data providers."""
    
    @abstractmethod
    def get_market_data(
        self,
        indication: str,
        region: str = "global"
    ) -> Optional[MarketDataPoint]:
        """
        Fetch market data for indication.
        
        Args:
            indication: Disease/therapeutic area
            region: Geographic region (global, us, eu, etc.)
        
        Returns:
            Market data or None if unavailable
        """
        pass
    
    @abstractmethod
    def is_available(self) -> bool:
        """Check if API is accessible."""
        pass


class IQVIAAdapter(MarketDataProvider):
    """
    Adapter for IQVIA market intelligence API.
    """
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.base_url = "https://api.iqvia.com/v1"  # Placeholder
        logger.info("Initialized IQVIA adapter")
    
    def is_available(self) -> bool:
        """Check if IQVIA API is accessible."""
        if not self.api_key:
            logger.warning("IQVIA API key not configured")
            return False
        
        # In production, make health check request
        # For now, just check key presence
        return True
    
    def get_market_data(
        self,
        indication: str,
        region: str = "global"
    ) -> Optional[MarketDataPoint]:
        """
        Fetch market data from IQVIA.
        
        Args:
            indication: Disease/therapeutic area
            region: Geographic region
        
        Returns:
            Market data or None
        """
        if not self.is_available():
            return None
        
        try:
            # Placeholder for actual API call
            # In production:
            # response = requests.get(
            #     f"{self.base_url}/market-intelligence/disease/{indication}",
            #     headers={"Authorization": f"Bearer {self.api_key}"},
            #     params={"region": region}
            # )
            # data = response.json()
            
            # For now, return None to indicate API not implemented
            logger.info(f"IQVIA API call would fetch: {indication} ({region})")
            return None
        
        except Exception as e:
            logger.error(f"IQVIA API error: {e}")
            return None


class EvaluatePharmaAdapter(MarketDataProvider):
    """
    Adapter for Evaluate Pharma market data API.
    """
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.base_url = "https://api.evaluate.com/v1"  # Placeholder
        logger.info("Initialized Evaluate Pharma adapter")
    
    def is_available(self) -> bool:
        """Check if Evaluate Pharma API is accessible."""
        if not self.api_key:
            logger.warning("Evaluate Pharma API key not configured")
            return False
        
        return True
    
    def get_market_data(
        self,
        indication: str,
        region: str = "global"
    ) -> Optional[MarketDataPoint]:
        """
        Fetch market data from Evaluate Pharma.
        
        Args:
            indication: Disease/therapeutic area
            region: Geographic region
        
        Returns:
            Market data or None
        """
        if not self.is_available():
            return None
        
        try:
            # Placeholder for actual API call
            logger.info(f"Evaluate Pharma API call would fetch: {indication} ({region})")
            return None
        
        except Exception as e:
            logger.error(f"Evaluate Pharma API error: {e}")
            return None


class SyntheticMarketDataProvider(MarketDataProvider):
    """
    Fallback synthetic market data generator.
    Used when paid APIs are unavailable.
    """
    
    def __init__(self):
        self.market_estimates = self._load_synthetic_data()
        logger.info("Initialized synthetic market data provider")
    
    def _load_synthetic_data(self) -> Dict:
        """Load synthetic market estimates."""
        return {
            "cancer": {
                "market_size_usd": 150_000_000_000,
                "cagr_percent": 7.5,
                "patient_population": 19_000_000,
                "market_share_leaders": ["Keytruda", "Opdivo", "Revlimid"]
            },
            "cardiovascular": {
                "market_size_usd": 45_000_000_000,
                "cagr_percent": 4.2,
                "patient_population": 85_000_000,
                "market_share_leaders": ["Eliquis", "Xarelto", "Plavix"]
            },
            "diabetes": {
                "market_size_usd": 60_000_000_000,
                "cagr_percent": 6.8,
                "patient_population": 537_000_000,
                "market_share_leaders": ["Ozempic", "Trulicity", "Lantus"]
            },
            "alzheimers": {
                "market_size_usd": 8_000_000_000,
                "cagr_percent": 12.5,
                "patient_population": 6_500_000,
                "market_share_leaders": ["Aricept", "Namenda", "Leqembi"]
            },
            "rheumatoid_arthritis": {
                "market_size_usd": 25_000_000_000,
                "cagr_percent": 3.5,
                "patient_population": 1_500_000,
                "market_share_leaders": ["Humira", "Enbrel", "Xeljanz"]
            }
        }
    
    def is_available(self) -> bool:
        """Synthetic data is always available."""
        return True
    
    def get_market_data(
        self,
        indication: str,
        region: str = "global"
    ) -> Optional[MarketDataPoint]:
        """
        Generate synthetic market data.
        
        Args:
            indication: Disease/therapeutic area
            region: Geographic region
        
        Returns:
            Synthetic market data
        """
        # Normalize indication name
        indication_key = indication.lower().replace(" ", "_")
        
        # Look for match in synthetic data
        data = None
        for key, value in self.market_estimates.items():
            if key in indication_key or indication_key in key:
                data = value
                break
        
        if not data:
            # Generic fallback
            data = {
                "market_size_usd": 10_000_000_000,
                "cagr_percent": 5.0,
                "patient_population": 1_000_000,
                "market_share_leaders": ["Leader 1", "Leader 2", "Leader 3"]
            }
        
        return MarketDataPoint(
            indication=indication,
            market_size_usd=data["market_size_usd"],
            cagr_percent=data["cagr_percent"],
            patient_population=data["patient_population"],
            market_share_leaders=data["market_share_leaders"],
            data_source="synthetic",
            last_updated=datetime.now(timezone.utc).isoformat(),
            confidence="low"
        )


class MarketDataAggregator:
    """
    Aggregates data from multiple providers with fallback logic.
    """
    
    def __init__(
        self,
        iqvia_api_key: Optional[str] = None,
        evaluate_api_key: Optional[str] = None
    ):
        # Initialize providers in priority order
        self.providers: List[MarketDataProvider] = []
        
        # Paid providers first
        if iqvia_api_key:
            self.providers.append(IQVIAAdapter(iqvia_api_key))
        
        if evaluate_api_key:
            self.providers.append(EvaluatePharmaAdapter(evaluate_api_key))
        
        # Synthetic fallback always available
        self.providers.append(SyntheticMarketDataProvider())
        
        logger.info(f"Initialized aggregator with {len(self.providers)} providers")
    
    def get_market_data(
        self,
        indication: str,
        region: str = "global"
    ) -> MarketDataPoint:
        """
        Get market data with automatic fallback.
        
        Args:
            indication: Disease/therapeutic area
            region: Geographic region
        
        Returns:
            Market data (guaranteed via synthetic fallback)
        """
        for provider in self.providers:
            if not provider.is_available():
                continue
            
            try:
                data = provider.get_market_data(indication, region)
                if data:
                    logger.info(f"Retrieved data from {data.data_source}")
                    return data
            except Exception as e:
                logger.warning(f"Provider failed: {e}, trying next...")
                continue
        
        # Should never reach here due to synthetic fallback
        # But add safety net
        return MarketDataPoint(
            indication=indication,
            data_source="fallback",
            confidence="low"
        )
    
    def get_stats(self) -> Dict:
        """Get aggregator statistics."""
        return {
            "total_providers": len(self.providers),
            "available_providers": sum(1 for p in self.providers if p.is_available()),
            "providers": [
                {
                    "type": type(p).__name__,
                    "available": p.is_available()
                }
                for p in self.providers
            ]
        }


def format_market_summary(data: MarketDataPoint) -> str:
    """
    Format market data for worker output.
    
    Args:
        data: Market data point
    
    Returns:
        Formatted summary text
    """
    lines = [f"**Market Analysis: {data.indication}**"]
    
    if data.market_size_usd:
        lines.append(f"- Market Size: ${data.market_size_usd / 1e9:.1f}B USD")
    
    if data.cagr_percent:
        lines.append(f"- Growth Rate (CAGR): {data.cagr_percent}%")
    
    if data.patient_population:
        lines.append(f"- Patient Population: {data.patient_population / 1e6:.1f}M")
    
    if data.market_share_leaders:
        leaders = ", ".join(data.market_share_leaders[:3])
        lines.append(f"- Market Leaders: {leaders}")
    
    lines.append(f"- Data Source: {data.data_source} ({data.confidence} confidence)")
    
    if data.last_updated:
        lines.append(f"- Last Updated: {data.last_updated[:10]}")
    
    return "\n".join(lines)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Demo usage
    print("=" * 60)
    print("Market Data API Integration Demo")
    print("=" * 60)
    
    # Create aggregator (no API keys = synthetic fallback)
    aggregator = MarketDataAggregator()
    
    # Get stats
    stats = aggregator.get_stats()
    print(f"\nAggregator stats: {stats}")
    
    # Test market data retrieval
    test_indications = ["cancer", "diabetes", "alzheimers"]
    
    for indication in test_indications:
        print(f"\n{'-' * 60}")
        data = aggregator.get_market_data(indication)
        print(format_market_summary(data))
    
    print("\n" + "=" * 60)
    print("To enable paid APIs:")
    print("1. Set environment variables:")
    print("   - IQVIA_API_KEY=your_key")
    print("   - EVALUATE_PHARMA_API_KEY=your_key")
    print("2. Update backend/app/config.py with API key fields")
    print("3. Modify agents/workers/market.py to use MarketDataAggregator")
