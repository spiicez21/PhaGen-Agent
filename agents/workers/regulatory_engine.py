"""
Multi-region regulatory rule engine for US/EU/India compliance.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class Region(Enum):
    """Regulatory regions."""
    US = "us"
    EU = "eu"
    INDIA = "india"
    GLOBAL = "global"


class ApprovalPathway(Enum):
    """Drug approval pathways."""
    STANDARD = "standard"
    ACCELERATED = "accelerated"
    BREAKTHROUGH = "breakthrough"
    ORPHAN = "orphan"
    PRIORITY_REVIEW = "priority_review"
    FAST_TRACK = "fast_track"


@dataclass
class RegulatoryRequirement:
    """Regulatory requirement for a region."""
    region: str
    phase: str  # preclinical, phase1, phase2, phase3, nda
    description: str
    estimated_duration_months: int
    success_rate_percent: float
    key_considerations: List[str]


@dataclass
class ApprovalPathwayInfo:
    """Information about approval pathway."""
    pathway: str
    region: str
    eligibility_criteria: List[str]
    benefits: List[str]
    timeline_reduction_months: int
    approval_probability_boost: float


@dataclass
class RegulatoryAnalysis:
    """Complete regulatory analysis result."""
    region: str
    indication: str
    recommended_pathway: str
    estimated_timeline_months: int
    approval_probability: float
    requirements: List[RegulatoryRequirement]
    special_pathways: List[ApprovalPathwayInfo]
    compliance_notes: List[str]


class RegulatoryRuleEngine:
    """
    Rule engine for multi-region regulatory compliance.
    """
    
    def __init__(self):
        self.us_rules = self._initialize_us_rules()
        self.eu_rules = self._initialize_eu_rules()
        self.india_rules = self._initialize_india_rules()
        logger.info("Initialized regulatory rule engine")
    
    def _initialize_us_rules(self) -> Dict:
        """Initialize FDA/US regulatory rules."""
        return {
            "standard_pathway": [
                RegulatoryRequirement(
                    region="US",
                    phase="IND Filing",
                    description="Investigational New Drug application",
                    estimated_duration_months=6,
                    success_rate_percent=90.0,
                    key_considerations=[
                        "Chemistry, Manufacturing, and Controls (CMC) data required",
                        "Toxicology studies in two species",
                        "Clinical protocol for Phase 1 studies"
                    ]
                ),
                RegulatoryRequirement(
                    region="US",
                    phase="Phase 1",
                    description="Safety and dosing in healthy volunteers (20-100 subjects)",
                    estimated_duration_months=12,
                    success_rate_percent=70.0,
                    key_considerations=[
                        "Safety and tolerability assessment",
                        "Pharmacokinetics/pharmacodynamics",
                        "MTD determination"
                    ]
                ),
                RegulatoryRequirement(
                    region="US",
                    phase="Phase 2",
                    description="Efficacy and safety in patient population (100-300 subjects)",
                    estimated_duration_months=24,
                    success_rate_percent=35.0,
                    key_considerations=[
                        "Proof of concept for efficacy",
                        "Dose optimization",
                        "Biomarker identification"
                    ]
                ),
                RegulatoryRequirement(
                    region="US",
                    phase="Phase 3",
                    description="Confirmatory efficacy trials (300-3000+ subjects)",
                    estimated_duration_months=36,
                    success_rate_percent=60.0,
                    key_considerations=[
                        "Randomized, controlled design",
                        "Statistical significance required",
                        "Long-term safety data"
                    ]
                ),
                RegulatoryRequirement(
                    region="US",
                    phase="NDA Review",
                    description="New Drug Application review by FDA",
                    estimated_duration_months=10,
                    success_rate_percent=85.0,
                    key_considerations=[
                        "Complete dossier submission",
                        "FDA Advisory Committee meeting possible",
                        "Risk Evaluation and Mitigation Strategy (REMS) may be required"
                    ]
                )
            ],
            "accelerated_pathways": {
                "breakthrough": ApprovalPathwayInfo(
                    pathway="Breakthrough Therapy",
                    region="US",
                    eligibility_criteria=[
                        "Treats serious or life-threatening condition",
                        "Preliminary clinical evidence shows substantial improvement",
                        "No adequate alternative exists"
                    ],
                    benefits=[
                        "More frequent FDA meetings",
                        "Rolling review of NDA",
                        "Priority Review designation",
                        "Potential for accelerated approval"
                    ],
                    timeline_reduction_months=12,
                    approval_probability_boost=0.15
                ),
                "fast_track": ApprovalPathwayInfo(
                    pathway="Fast Track",
                    region="US",
                    eligibility_criteria=[
                        "Treats serious condition",
                        "Addresses unmet medical need"
                    ],
                    benefits=[
                        "Rolling review of applications",
                        "More frequent communication with FDA",
                        "Eligibility for accelerated approval/priority review"
                    ],
                    timeline_reduction_months=6,
                    approval_probability_boost=0.10
                ),
                "orphan": ApprovalPathwayInfo(
                    pathway="Orphan Drug",
                    region="US",
                    eligibility_criteria=[
                        "Disease affects <200,000 patients in US",
                        "Or affects >200,000 but no reasonable expectation of cost recovery"
                    ],
                    benefits=[
                        "7-year market exclusivity",
                        "Tax credits for clinical trial costs",
                        "FDA fee waivers",
                        "Protocol assistance"
                    ],
                    timeline_reduction_months=0,
                    approval_probability_boost=0.05
                )
            }
        }
    
    def _initialize_eu_rules(self) -> Dict:
        """Initialize EMA/EU regulatory rules."""
        return {
            "standard_pathway": [
                RegulatoryRequirement(
                    region="EU",
                    phase="CTA Submission",
                    description="Clinical Trial Application to national authority",
                    estimated_duration_months=6,
                    success_rate_percent=88.0,
                    key_considerations=[
                        "IMPD (Investigational Medicinal Product Dossier)",
                        "Ethics committee approval required",
                        "EudraCT number registration"
                    ]
                ),
                RegulatoryRequirement(
                    region="EU",
                    phase="Phase 1",
                    description="Safety and tolerability (20-100 subjects)",
                    estimated_duration_months=12,
                    success_rate_percent=68.0,
                    key_considerations=[
                        "GMP-manufactured drug required",
                        "Pharmacovigilance reporting",
                        "Subject insurance mandatory"
                    ]
                ),
                RegulatoryRequirement(
                    region="EU",
                    phase="Phase 2",
                    description="Dose-finding and efficacy (100-300 subjects)",
                    estimated_duration_months=24,
                    success_rate_percent=33.0,
                    key_considerations=[
                        "Multi-country trials common",
                        "HTA bodies early engagement recommended",
                        "Pediatric investigation plan (PIP) may be required"
                    ]
                ),
                RegulatoryRequirement(
                    region="EU",
                    phase="Phase 3",
                    description="Confirmatory trials (300-3000+ subjects)",
                    estimated_duration_months=36,
                    success_rate_percent=58.0,
                    key_considerations=[
                        "Centralized or decentralized procedure",
                        "Quality of life endpoints valued",
                        "Real-world evidence increasingly considered"
                    ]
                ),
                RegulatoryRequirement(
                    region="EU",
                    phase="MAA Review",
                    description="Marketing Authorization Application review by EMA",
                    estimated_duration_months=12,
                    success_rate_percent=82.0,
                    key_considerations=[
                        "Centralized procedure for most new drugs",
                        "CHMP opinion required",
                        "National pricing/reimbursement separate process"
                    ]
                )
            ],
            "accelerated_pathways": {
                "prime": ApprovalPathwayInfo(
                    pathway="PRIME (Priority Medicines)",
                    region="EU",
                    eligibility_criteria=[
                        "Major public health interest",
                        "Therapeutic innovation",
                        "Targets unmet medical need"
                    ],
                    benefits=[
                        "Early dialogue with EMA",
                        "Accelerated assessment (150 days vs 210)",
                        "Appointment of rapporteur"
                    ],
                    timeline_reduction_months=8,
                    approval_probability_boost=0.12
                ),
                "conditional": ApprovalPathwayInfo(
                    pathway="Conditional Authorization",
                    region="EU",
                    eligibility_criteria=[
                        "Serious disease with unmet need",
                        "Benefit-risk balance positive",
                        "Complete data will be provided post-approval"
                    ],
                    benefits=[
                        "Earlier market access",
                        "Annual renewal until full data",
                        "Conversion to standard MA possible"
                    ],
                    timeline_reduction_months=12,
                    approval_probability_boost=0.10
                )
            }
        }
    
    def _initialize_india_rules(self) -> Dict:
        """Initialize CDSCO/India regulatory rules."""
        return {
            "standard_pathway": [
                RegulatoryRequirement(
                    region="India",
                    phase="IND Application",
                    description="Investigational New Drug permission from CDSCO",
                    estimated_duration_months=4,
                    success_rate_percent=85.0,
                    key_considerations=[
                        "Prior approval in ICH countries may expedite",
                        "Ethics Committee clearance mandatory",
                        "DCGI approval required"
                    ]
                ),
                RegulatoryRequirement(
                    region="India",
                    phase="Phase 1-3",
                    description="Clinical trials in Indian population",
                    estimated_duration_months=48,
                    success_rate_percent=55.0,
                    key_considerations=[
                        "Local clinical data increasingly required",
                        "Post-marketing surveillance stringent",
                        "Price control for essential medicines"
                    ]
                ),
                RegulatoryRequirement(
                    region="India",
                    phase="NDA Review",
                    description="New Drug Application to CDSCO",
                    estimated_duration_months=12,
                    success_rate_percent=75.0,
                    key_considerations=[
                        "Manufacturing site inspection",
                        "Price negotiation with NPPA",
                        "Patent linkage considerations"
                    ]
                )
            ],
            "accelerated_pathways": {
                "fast_track": ApprovalPathwayInfo(
                    pathway="Fast Track Approval",
                    region="India",
                    eligibility_criteria=[
                        "Already approved in US/EU/Japan/Australia/Canada",
                        "Addresses unmet medical need in India"
                    ],
                    benefits=[
                        "Waiver of local clinical trials possible",
                        "Shorter review timeline",
                        "Priority review by CDSCO"
                    ],
                    timeline_reduction_months=18,
                    approval_probability_boost=0.20
                )
            }
        }
    
    def analyze_pathway(
        self,
        indication: str,
        region: Region,
        is_orphan: bool = False,
        is_life_threatening: bool = False,
        has_foreign_approval: bool = False
    ) -> RegulatoryAnalysis:
        """
        Analyze regulatory pathway for molecule.
        
        Args:
            indication: Disease indication
            region: Target regulatory region
            is_orphan: Whether qualifies as orphan drug
            is_life_threatening: Whether treats life-threatening condition
            has_foreign_approval: Whether approved in other major markets
        
        Returns:
            Complete regulatory analysis
        """
        # Get appropriate rule set
        if region == Region.US:
            rules = self.us_rules
        elif region == Region.EU:
            rules = self.eu_rules
        elif region == Region.INDIA:
            rules = self.india_rules
        else:
            # Default to US
            rules = self.us_rules
        
        # Determine recommended pathway
        special_pathways = []
        recommended_pathway = "Standard"
        timeline_reduction = 0
        probability_boost = 0.0
        
        # Check eligibility for special pathways
        for pathway_name, pathway_info in rules.get("accelerated_pathways", {}).items():
            eligible = False
            
            if pathway_name == "orphan" and is_orphan:
                eligible = True
            elif pathway_name == "breakthrough" and is_life_threatening:
                eligible = True
            elif pathway_name == "fast_track" and (is_life_threatening or has_foreign_approval):
                eligible = True
            elif pathway_name == "prime" and is_life_threatening:
                eligible = True
            elif pathway_name == "conditional" and is_life_threatening:
                eligible = True
            
            if eligible:
                special_pathways.append(pathway_info)
                # Use most beneficial pathway
                if pathway_info.timeline_reduction_months > timeline_reduction:
                    recommended_pathway = pathway_info.pathway
                    timeline_reduction = pathway_info.timeline_reduction_months
                    probability_boost = pathway_info.approval_probability_boost
        
        # Calculate timeline and probability
        standard_requirements = rules["standard_pathway"]
        base_timeline = sum(r.estimated_duration_months for r in standard_requirements)
        total_timeline = max(24, base_timeline - timeline_reduction)  # Minimum 2 years
        
        # Calculate success probability (product of phase success rates)
        base_probability = 1.0
        for req in standard_requirements:
            base_probability *= (req.success_rate_percent / 100.0)
        
        total_probability = min(0.95, base_probability + probability_boost)
        
        # Generate compliance notes
        compliance_notes = [
            f"Estimated timeline: {total_timeline} months from IND to approval",
            f"Overall success probability: {total_probability*100:.1f}%"
        ]
        
        if special_pathways:
            compliance_notes.append(
                f"Eligible for {len(special_pathways)} accelerated pathway(s)"
            )
        
        if region == Region.EU:
            compliance_notes.append(
                "Separate HTA and reimbursement process required per country"
            )
        elif region == Region.INDIA:
            compliance_notes.append(
                "Price control under DPCO may apply for essential medicines"
            )
        
        return RegulatoryAnalysis(
            region=region.value,
            indication=indication,
            recommended_pathway=recommended_pathway,
            estimated_timeline_months=total_timeline,
            approval_probability=round(total_probability, 3),
            requirements=standard_requirements,
            special_pathways=special_pathways,
            compliance_notes=compliance_notes
        )
    
    def compare_regions(
        self,
        indication: str,
        regions: List[Region],
        **kwargs
    ) -> Dict[str, RegulatoryAnalysis]:
        """
        Compare regulatory pathways across regions.
        
        Args:
            indication: Disease indication
            regions: List of regions to compare
            **kwargs: Additional parameters for analyze_pathway
        
        Returns:
            Mapping of region to analysis
        """
        results = {}
        
        for region in regions:
            results[region.value] = self.analyze_pathway(
                indication=indication,
                region=region,
                **kwargs
            )
        
        return results


def format_regulatory_summary(analysis: RegulatoryAnalysis) -> str:
    """
    Format regulatory analysis for output.
    
    Args:
        analysis: Regulatory analysis
    
    Returns:
        Formatted summary text
    """
    lines = [
        f"**Regulatory Analysis: {analysis.indication} ({analysis.region.upper()})**",
        f"",
        f"**Recommended Pathway:** {analysis.recommended_pathway}",
        f"**Estimated Timeline:** {analysis.estimated_timeline_months} months",
        f"**Approval Probability:** {analysis.approval_probability*100:.1f}%",
        f""
    ]
    
    if analysis.special_pathways:
        lines.append("**Eligible Special Pathways:**")
        for pathway in analysis.special_pathways:
            lines.append(f"- {pathway.pathway}")
            lines.append(f"  - Timeline Reduction: {pathway.timeline_reduction_months} months")
            lines.append(f"  - Benefits: {', '.join(pathway.benefits[:2])}")
        lines.append("")
    
    lines.append("**Key Compliance Notes:**")
    for note in analysis.compliance_notes:
        lines.append(f"- {note}")
    
    return "\n".join(lines)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 60)
    print("Multi-Region Regulatory Rule Engine Demo")
    print("=" * 60)
    
    engine = RegulatoryRuleEngine()
    
    # Test case: orphan oncology drug
    print("\nTest Case: Orphan Oncology Drug")
    print("-" * 60)
    
    comparison = engine.compare_regions(
        indication="Rare Leukemia",
        regions=[Region.US, Region.EU, Region.INDIA],
        is_orphan=True,
        is_life_threatening=True,
        has_foreign_approval=False
    )
    
    for region, analysis in comparison.items():
        print(f"\n{format_regulatory_summary(analysis)}")
        print("-" * 60)
