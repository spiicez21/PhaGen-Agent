"""
Patent claim-level semantic matching and freedom-to-operate analysis.
"""
from __future__ import annotations

import re
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class PatentClaim:
    """Represents a parsed patent claim."""
    claim_number: int
    claim_type: str  # 'independent' or 'dependent'
    claim_text: str
    depends_on: Optional[int] = None
    blocking_score: float = 0.0
    semantic_similarity: float = 0.0


@dataclass
class ClaimMap:
    """Visual representation of claim dependencies."""
    patent_id: str
    independent_claims: List[PatentClaim]
    dependent_claims: List[PatentClaim]
    blocking_claims: List[PatentClaim]
    fto_risk_score: float


class PatentClaimAnalyzer:
    """
    Analyzes patent claims at granular level for FTO assessment.
    """
    
    def __init__(self, retriever=None):
        """
        Initialize claim analyzer.
        
        Args:
            retriever: Optional RAG retriever for semantic search
        """
        self.retriever = retriever
    
    def extract_claims(self, patent_text: str) -> List[PatentClaim]:
        """
        Parse patent claims from full patent text.
        
        Args:
            patent_text: Full patent document text
        
        Returns:
            List of parsed claims
        """
        claims = []
        
        # Find claims section
        claims_match = re.search(
            r'(?:CLAIMS?|What is claimed is:)(.*?)(?:ABSTRACT|DESCRIPTION|$)',
            patent_text,
            re.IGNORECASE | re.DOTALL
        )
        
        if not claims_match:
            logger.warning("No claims section found in patent text")
            return claims
        
        claims_text = claims_match.group(1)
        
        # Extract individual claims (numbered)
        claim_pattern = r'(\d+)\.\s+(.*?)(?=\d+\.|$)'
        claim_matches = re.finditer(claim_pattern, claims_text, re.DOTALL)
        
        for match in claim_matches:
            claim_num = int(match.group(1))
            claim_text = match.group(2).strip()
            
            # Determine if dependent claim
            depends_match = re.search(r'claim (\d+)', claim_text, re.IGNORECASE)
            is_dependent = depends_match is not None
            depends_on = int(depends_match.group(1)) if depends_match else None
            
            claims.append(PatentClaim(
                claim_number=claim_num,
                claim_type='dependent' if is_dependent else 'independent',
                claim_text=claim_text,
                depends_on=depends_on
            ))
        
        logger.info(f"Extracted {len(claims)} claims")
        return claims
    
    def analyze_blocking_risk(
        self,
        claims: List[PatentClaim],
        molecule_description: str,
        smiles: Optional[str] = None
    ) -> List[PatentClaim]:
        """
        Score claims by blocking risk for given molecule.
        
        Args:
            claims: List of patent claims
            molecule_description: Description of target molecule
            smiles: Optional SMILES string
        
        Returns:
            Claims sorted by blocking risk (highest first)
        """
        if not claims:
            return []
        
        # Build query from molecule info
        query = molecule_description
        if smiles:
            query += f" SMILES: {smiles}"
        
        # Score each claim
        for claim in claims:
            # Basic keyword matching
            keywords = self._extract_keywords(molecule_description)
            matches = sum(1 for kw in keywords if kw.lower() in claim.claim_text.lower())
            
            # Higher score = higher blocking risk
            claim.blocking_score = min(1.0, matches / max(1, len(keywords)))
            
            # Semantic similarity (if retriever available)
            if self.retriever:
                try:
                    # Simple embedding cosine similarity
                    claim.semantic_similarity = self._compute_similarity(
                        query, claim.claim_text
                    )
                    # Combine keyword and semantic scores
                    claim.blocking_score = (
                        0.6 * claim.blocking_score + 
                        0.4 * claim.semantic_similarity
                    )
                except Exception as e:
                    logger.warning(f"Semantic similarity failed: {e}")
        
        # Sort by blocking score
        return sorted(claims, key=lambda c: c.blocking_score, reverse=True)
    
    def generate_claim_map(
        self,
        patent_id: str,
        claims: List[PatentClaim],
        molecule_description: str
    ) -> ClaimMap:
        """
        Generate visual claim dependency map.
        
        Args:
            patent_id: Patent identifier
            claims: List of analyzed claims
            molecule_description: Target molecule description
        
        Returns:
            ClaimMap with dependencies and FTO risk
        """
        if not claims:
            return ClaimMap(
                patent_id=patent_id,
                independent_claims=[],
                dependent_claims=[],
                blocking_claims=[],
                fto_risk_score=0.0
            )
        
        # Analyze blocking risk first
        analyzed = self.analyze_blocking_risk(claims, molecule_description)
        
        # Separate independent and dependent
        independent = [c for c in analyzed if c.claim_type == 'independent']
        dependent = [c for c in analyzed if c.claim_type == 'dependent']
        
        # Identify blocking claims (score > 0.5)
        blocking = [c for c in analyzed if c.blocking_score > 0.5]
        
        # Calculate overall FTO risk
        if blocking:
            # Independent blocking claims are highest risk
            independent_blocking = [c for c in blocking if c.claim_type == 'independent']
            if independent_blocking:
                fto_risk = max(c.blocking_score for c in independent_blocking)
            else:
                fto_risk = max(c.blocking_score for c in blocking) * 0.7
        else:
            fto_risk = 0.0
        
        return ClaimMap(
            patent_id=patent_id,
            independent_claims=independent,
            dependent_claims=dependent,
            blocking_claims=blocking,
            fto_risk_score=round(fto_risk, 3)
        )
    
    def generate_fto_report(
        self,
        claim_maps: List[ClaimMap],
        molecule: str
    ) -> Dict:
        """
        Generate freedom-to-operate report.
        
        Args:
            claim_maps: List of claim maps for relevant patents
            molecule: Target molecule name
        
        Returns:
            FTO report with risk assessment
        """
        if not claim_maps:
            return {
                "molecule": molecule,
                "overall_risk": "LOW",
                "blocking_patents": [],
                "recommendations": ["No blocking claims identified"]
            }
        
        # Sort by risk
        high_risk_maps = [m for m in claim_maps if m.fto_risk_score > 0.7]
        medium_risk_maps = [m for m in claim_maps if 0.4 <= m.fto_risk_score <= 0.7]
        
        # Determine overall risk
        if high_risk_maps:
            overall_risk = "HIGH"
        elif medium_risk_maps:
            overall_risk = "MEDIUM"
        else:
            overall_risk = "LOW"
        
        # Generate recommendations
        recommendations = []
        if high_risk_maps:
            recommendations.append(
                f"Conduct detailed FTO analysis for {len(high_risk_maps)} high-risk patents"
            )
            recommendations.append("Consider design-around strategies or licensing")
        elif medium_risk_maps:
            recommendations.append("Monitor medium-risk patents for potential conflicts")
        else:
            recommendations.append("No significant FTO risks identified")
        
        # Build blocking patents list
        blocking_patents = []
        for cm in high_risk_maps + medium_risk_maps:
            blocking_patents.append({
                "patent_id": cm.patent_id,
                "risk_score": cm.fto_risk_score,
                "blocking_claims": [
                    {
                        "claim_number": c.claim_number,
                        "claim_type": c.claim_type,
                        "blocking_score": c.blocking_score,
                        "excerpt": c.claim_text[:200] + "..."
                    }
                    for c in cm.blocking_claims[:3]  # Top 3 blocking claims
                ]
            })
        
        return {
            "molecule": molecule,
            "overall_risk": overall_risk,
            "risk_score": max(m.fto_risk_score for m in claim_maps) if claim_maps else 0.0,
            "total_patents_analyzed": len(claim_maps),
            "high_risk_patents": len(high_risk_maps),
            "medium_risk_patents": len(medium_risk_maps),
            "blocking_patents": blocking_patents,
            "recommendations": recommendations
        }
    
    def _extract_keywords(self, text: str) -> List[str]:
        """Extract key terms from molecule description."""
        # Simple keyword extraction
        stopwords = {'a', 'an', 'the', 'is', 'are', 'for', 'of', 'in', 'to', 'and', 'or'}
        words = re.findall(r'\b\w{4,}\b', text.lower())
        return [w for w in words if w not in stopwords]
    
    def _compute_similarity(self, text1: str, text2: str) -> float:
        """Compute semantic similarity between two texts."""
        if not self.retriever:
            return 0.0
        
        # Placeholder for semantic similarity
        # In production, use retriever's embedding model
        try:
            # Simple Jaccard similarity as fallback
            words1 = set(self._extract_keywords(text1))
            words2 = set(self._extract_keywords(text2))
            
            if not words1 or not words2:
                return 0.0
            
            intersection = len(words1 & words2)
            union = len(words1 | words2)
            
            return intersection / union if union > 0 else 0.0
        except Exception:
            return 0.0


def analyze_patent_claims_for_molecule(
    patent_documents: List[Dict],
    molecule: str,
    molecule_description: str,
    smiles: Optional[str] = None
) -> Dict:
    """
    Entry point for patent claim analysis.
    
    Args:
        patent_documents: List of patent documents from retrieval
        molecule: Molecule name
        molecule_description: Description/mechanism
        smiles: Optional SMILES string
    
    Returns:
        FTO analysis report
    """
    analyzer = PatentClaimAnalyzer()
    claim_maps = []
    
    for doc in patent_documents[:10]:  # Analyze top 10 patents
        patent_text = doc.get('text', '') or doc.get('content', '')
        patent_id = doc.get('patent_id') or doc.get('url', 'UNKNOWN')
        
        if not patent_text:
            continue
        
        # Extract and analyze claims
        claims = analyzer.extract_claims(patent_text)
        if claims:
            claim_map = analyzer.generate_claim_map(
                patent_id=patent_id,
                claims=claims,
                molecule_description=molecule_description
            )
            claim_maps.append(claim_map)
    
    # Generate FTO report
    return analyzer.generate_fto_report(claim_maps, molecule)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Example usage
    sample_patent = """
    CLAIMS
    1. A pharmaceutical composition comprising a compound of formula X.
    2. The composition of claim 1, wherein the compound inhibits kinase Y.
    3. The composition of claim 1, for use in treating cancer.
    """
    
    analyzer = PatentClaimAnalyzer()
    claims = analyzer.extract_claims(sample_patent)
    
    print(f"Extracted {len(claims)} claims:")
    for claim in claims:
        print(f"  Claim {claim.claim_number} ({claim.claim_type})")
    
    # Analyze blocking risk
    analyzed = analyzer.analyze_blocking_risk(
        claims,
        molecule_description="kinase Y inhibitor for cancer treatment"
    )
    
    print(f"\nBlocking risk analysis:")
    for claim in analyzed[:3]:
        print(f"  Claim {claim.claim_number}: {claim.blocking_score:.2f}")
