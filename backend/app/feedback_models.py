"""
Database models for evidence feedback and active learning.
"""
from __future__ import annotations

from datetime import datetime, timezone
from sqlalchemy import Column, String, Integer, Float, DateTime, Text, ForeignKey, Index
from sqlalchemy.orm import relationship
from .database import Base


class EvidenceFeedback(Base):
    """
    User feedback on evidence quality for active learning.
    Used to fine-tune reranker weights and improve retrieval relevance.
    """
    __tablename__ = "evidence_feedback"
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    job_id = Column(String, ForeignKey("jobs.job_id"), nullable=False, index=True)
    evidence_id = Column(String, nullable=False, index=True)  # Unique evidence identifier
    evidence_type = Column(String, nullable=False)  # clinical, literature, patent, market
    feedback_type = Column(String, nullable=False)  # upvote, downvote, flag
    feedback_score = Column(Float, nullable=False)  # +1.0 for upvote, -1.0 for downvote, 0.0 for flag
    user_id = Column(String, nullable=True, index=True)  # Optional user identifier
    comment = Column(Text, nullable=True)  # Optional user comment
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc), nullable=False)
    
    # Relationship to job
    job = relationship("Job", back_populates="feedback")
    
    __table_args__ = (
        Index("idx_evidence_feedback_job_evidence", "job_id", "evidence_id"),
        Index("idx_evidence_feedback_type_score", "evidence_type", "feedback_score"),
    )
    
    def __repr__(self):
        return f"<EvidenceFeedback(id={self.id}, job_id={self.job_id}, type={self.feedback_type}, score={self.feedback_score})>"
