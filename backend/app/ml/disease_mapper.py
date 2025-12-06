"""
Molecule-disease mapping model for proactive repurposing suggestions.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class DiseaseMapping:
    """Molecule-disease mapping result."""
    disease: str
    disease_id: str  # MeSH, ICD-10, or OMIM
    confidence: float
    evidence_sources: List[str]
    mechanism: Optional[str] = None


class MoleculeDiseaseMapper:
    """
    Maps molecules to potential disease indications for repurposing.
    Uses pre-computed embeddings and similarity search.
    """
    
    def __init__(self, data_dir: Path = Path("indexes/disease_mappings")):
        self.data_dir = data_dir
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.mappings_file = self.data_dir / "molecule_disease_mappings.json"
        self.disease_index = self._load_disease_index()
    
    def _load_disease_index(self) -> Dict:
        """Load disease mapping index."""
        if self.mappings_file.exists():
            with open(self.mappings_file) as f:
                return json.load(f)
        
        # Initialize with common disease categories
        return {
            "cancer": {
                "subcategories": ["breast_cancer", "lung_cancer", "leukemia", "lymphoma"],
                "keywords": ["tumor", "carcinoma", "oncology", "metastasis"],
                "mesh_ids": ["D009369", "D001943", "D008175"]
            },
            "cardiovascular": {
                "subcategories": ["hypertension", "heart_failure", "arrhythmia"],
                "keywords": ["cardiac", "heart", "vascular", "blood_pressure"],
                "mesh_ids": ["D002318", "D006973", "D006333"]
            },
            "neurodegenerative": {
                "subcategories": ["alzheimers", "parkinsons", "als", "huntingtons"],
                "keywords": ["neurodegeneration", "dementia", "cognitive_decline"],
                "mesh_ids": ["D000544", "D010300", "D016472"]
            },
            "inflammatory": {
                "subcategories": ["rheumatoid_arthritis", "crohns", "psoriasis"],
                "keywords": ["inflammation", "autoimmune", "cytokine"],
                "mesh_ids": ["D001172", "D003424", "D011565"]
            },
            "metabolic": {
                "subcategories": ["diabetes", "obesity", "nafld"],
                "keywords": ["glucose", "insulin", "metabolic_syndrome"],
                "mesh_ids": ["D003920", "D009765", "D024821"]
            },
            "infectious": {
                "subcategories": ["viral", "bacterial", "fungal"],
                "keywords": ["infection", "pathogen", "antimicrobial"],
                "mesh_ids": ["D003141", "D001424", "D009181"]
            },
        }
    
    def _save_disease_index(self):
        """Persist disease index."""
        with open(self.mappings_file, 'w') as f:
            json.dump(self.disease_index, f, indent=2)
    
    def predict_diseases(
        self,
        molecule: str,
        smiles: Optional[str] = None,
        evidence: Optional[List[Dict]] = None,
        top_k: int = 5
    ) -> List[DiseaseMapping]:
        """
        Predict disease indications for a molecule.
        
        Args:
            molecule: Molecule name
            smiles: SMILES string (for structure-based prediction)
            evidence: Existing evidence from workers
            top_k: Number of predictions to return
        
        Returns:
            List of disease mappings sorted by confidence
        """
        predictions = []
        
        if not evidence:
            logger.warning(f"No evidence provided for {molecule}, returning empty predictions")
            return predictions
        
        # Analyze evidence for disease keywords
        for category, data in self.disease_index.items():
            score = 0.0
            matched_keywords = []
            sources = []
            
            # Check clinical trial evidence
            for ev in evidence:
                if ev.get("type") == "clinical":
                    text = (ev.get("text", "") + " " + ev.get("summary", "")).lower()
                    
                    for keyword in data["keywords"]:
                        if keyword.lower() in text:
                            score += 0.3
                            matched_keywords.append(keyword)
                            sources.append(ev.get("url", ""))
                
                # Check literature evidence
                elif ev.get("type") == "literature":
                    text = (ev.get("text", "") + " " + ev.get("summary", "")).lower()
                    
                    for keyword in data["keywords"]:
                        if keyword.lower() in text:
                            score += 0.2
                            matched_keywords.append(keyword)
                            sources.append(ev.get("url", ""))
            
            if score > 0.1:
                # Normalize confidence
                confidence = min(1.0, score / len(data["keywords"]))
                
                predictions.append(DiseaseMapping(
                    disease=category.replace("_", " ").title(),
                    disease_id=data["mesh_ids"][0] if data["mesh_ids"] else "UNKNOWN",
                    confidence=round(confidence, 3),
                    evidence_sources=list(set(sources))[:3],
                    mechanism=f"Targets identified in {len(matched_keywords)} pathways"
                ))
        
        # Sort by confidence
        predictions.sort(key=lambda x: x.confidence, reverse=True)
        
        return predictions[:top_k]
    
    def add_known_mapping(
        self,
        molecule: str,
        disease: str,
        disease_id: str,
        confidence: float,
        source: str
    ):
        """
        Add a known molecule-disease mapping to the index.
        
        Args:
            molecule: Molecule name
            disease: Disease name
            disease_id: MeSH/ICD code
            confidence: Confidence score (0-1)
            source: Evidence source
        """
        # Store in custom mappings
        if "custom_mappings" not in self.disease_index:
            self.disease_index["custom_mappings"] = {}
        
        if molecule not in self.disease_index["custom_mappings"]:
            self.disease_index["custom_mappings"][molecule] = []
        
        self.disease_index["custom_mappings"][molecule].append({
            "disease": disease,
            "disease_id": disease_id,
            "confidence": confidence,
            "source": source
        })
        
        self._save_disease_index()
        logger.info(f"Added mapping: {molecule} -> {disease} ({confidence:.2f})")
    
    def get_stats(self) -> Dict:
        """Get mapper statistics."""
        custom_count = len(self.disease_index.get("custom_mappings", {}))
        
        return {
            "total_categories": len([k for k in self.disease_index.keys() if k != "custom_mappings"]),
            "custom_mappings": custom_count,
            "data_dir": str(self.data_dir),
        }


def generate_repurposing_suggestions(job_id: str, result: Dict) -> List[Dict]:
    """
    Generate repurposing suggestions for a completed job.
    
    Args:
        job_id: Job identifier
        result: Master agent result with worker outputs
    
    Returns:
        List of repurposing suggestions
    """
    mapper = MoleculeDiseaseMapper()
    
    molecule = result.get("molecule", "Unknown")
    smiles = result.get("smiles")
    
    # Extract evidence from worker outputs
    evidence = []
    for worker_type in ["clinical", "literature", "patent", "market"]:
        worker_data = result.get(worker_type, {})
        if isinstance(worker_data, dict) and "evidence" in worker_data:
            for ev in worker_data["evidence"]:
                evidence.append({
                    "type": worker_type,
                    "text": ev.get("text", ""),
                    "summary": ev.get("summary", ""),
                    "url": ev.get("url", ""),
                })
    
    # Predict diseases
    predictions = mapper.predict_diseases(
        molecule=molecule,
        smiles=smiles,
        evidence=evidence,
        top_k=5
    )
    
    suggestions = []
    for pred in predictions:
        suggestions.append({
            "disease": pred.disease,
            "disease_id": pred.disease_id,
            "confidence": pred.confidence,
            "rationale": f"Found evidence across {len(pred.evidence_sources)} sources",
            "mechanism": pred.mechanism,
            "sources": pred.evidence_sources,
        })
    
    return suggestions


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Example usage
    mapper = MoleculeDiseaseMapper()
    
    # Add known mapping
    mapper.add_known_mapping(
        molecule="Aspirin",
        disease="Cardiovascular Disease",
        disease_id="D002318",
        confidence=0.95,
        source="FDA Label"
    )
    
    # Test prediction
    test_evidence = [
        {
            "type": "clinical",
            "text": "Phase 2 trial showing efficacy in heart failure patients",
            "summary": "Improved cardiac function",
            "url": "https://clinicaltrials.gov/ct2/show/NCT12345"
        }
    ]
    
    predictions = mapper.predict_diseases(
        molecule="Test-Molecule",
        evidence=test_evidence
    )
    
    print(f"Predictions: {predictions}")
    print(f"Stats: {mapper.get_stats()}")
