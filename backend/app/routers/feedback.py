"""
Evidence feedback router for active learning.
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
from typing import Optional, List
from datetime import datetime

from ..database import SessionLocal
from ..models import EvidenceFeedback, Job
from ..config import get_settings

settings = get_settings()
router = APIRouter(prefix=f"{settings.api_prefix}/feedback", tags=["feedback"])


class FeedbackCreate(BaseModel):
    """Request model for creating evidence feedback."""
    job_id: str = Field(..., description="Job ID containing the evidence")
    evidence_id: str = Field(..., description="Unique identifier for the evidence")
    evidence_type: str = Field(..., description="Type of evidence: clinical, literature, patent, market")
    feedback_type: str = Field(..., description="Feedback type: upvote, downvote, flag")
    user_id: Optional[str] = Field(None, description="Optional user identifier")
    comment: Optional[str] = Field(None, description="Optional comment explaining feedback")
    
    class Config:
        json_schema_extra = {
            "example": {
                "job_id": "123e4567-e89b-12d3-a456-426614174000",
                "evidence_id": "NCT01907900",
                "evidence_type": "clinical",
                "feedback_type": "upvote",
                "user_id": "user@example.com",
                "comment": "Highly relevant Phase 2 trial"
            }
        }


class FeedbackResponse(BaseModel):
    """Response model for feedback."""
    id: str
    job_id: str
    evidence_id: str
    evidence_type: str
    feedback_type: str
    feedback_score: float
    user_id: Optional[str]
    comment: Optional[str]
    created_at: datetime
    
    class Config:
        from_attributes = True


@router.post("", response_model=FeedbackResponse, status_code=201)
def submit_feedback(request: FeedbackCreate):
    """
    Submit feedback on evidence quality.
    
    - **upvote**: Evidence is highly relevant and useful (+1.0)
    - **downvote**: Evidence is not relevant or low quality (-1.0)
    - **flag**: Evidence has issues (misleading, outdated, etc.) (0.0)
    """
    db = SessionLocal()
    try:
        # Verify job exists
        job = db.query(Job).filter(Job.id == request.job_id).first()
        if not job:
            raise HTTPException(status_code=404, detail="Job not found")
        
        # Map feedback type to score
        score_map = {
            "upvote": 1.0,
            "downvote": -1.0,
            "flag": 0.0,
        }
        
        if request.feedback_type not in score_map:
            raise HTTPException(
                status_code=400,
                detail=f"Invalid feedback_type. Must be one of: {list(score_map.keys())}"
            )
        
        feedback = EvidenceFeedback(
            job_id=request.job_id,
            evidence_id=request.evidence_id,
            evidence_type=request.evidence_type,
            feedback_type=request.feedback_type,
            feedback_score=score_map[request.feedback_type],
            user_id=request.user_id,
            comment=request.comment,
        )
        
        db.add(feedback)
        db.commit()
        db.refresh(feedback)
        
        return FeedbackResponse.model_validate(feedback)
        
    finally:
        db.close()


@router.get("/job/{job_id}", response_model=List[FeedbackResponse])
def get_job_feedback(job_id: str):
    """Get all feedback for a specific job."""
    db = SessionLocal()
    try:
        feedback_list = db.query(EvidenceFeedback).filter(
            EvidenceFeedback.job_id == job_id
        ).order_by(EvidenceFeedback.created_at.desc()).all()
        
        return [FeedbackResponse.model_validate(f) for f in feedback_list]
        
    finally:
        db.close()


@router.get("/evidence/{evidence_id}", response_model=List[FeedbackResponse])
def get_evidence_feedback(
    evidence_id: str,
    evidence_type: Optional[str] = Query(None, description="Filter by evidence type")
):
    """Get all feedback for a specific evidence across all jobs."""
    db = SessionLocal()
    try:
        query = db.query(EvidenceFeedback).filter(
            EvidenceFeedback.evidence_id == evidence_id
        )
        
        if evidence_type:
            query = query.filter(EvidenceFeedback.evidence_type == evidence_type)
        
        feedback_list = query.order_by(EvidenceFeedback.created_at.desc()).all()
        
        return [FeedbackResponse.model_validate(f) for f in feedback_list]
        
    finally:
        db.close()


@router.get("/stats")
def get_feedback_stats():
    """
    Get aggregate feedback statistics for reranker training.
    Returns evidence types with best/worst average scores.
    """
    db = SessionLocal()
    try:
        from sqlalchemy import func
        
        # Aggregate by evidence type
        type_stats = db.query(
            EvidenceFeedback.evidence_type,
            func.count(EvidenceFeedback.id).label("total_feedback"),
            func.avg(EvidenceFeedback.feedback_score).label("avg_score"),
            func.sum(func.case((EvidenceFeedback.feedback_type == "upvote", 1), else_=0)).label("upvotes"),
            func.sum(func.case((EvidenceFeedback.feedback_type == "downvote", 1), else_=0)).label("downvotes"),
            func.sum(func.case((EvidenceFeedback.feedback_type == "flag", 1), else_=0)).label("flags"),
        ).group_by(EvidenceFeedback.evidence_type).all()
        
        # Aggregate by evidence ID (top/bottom evidence)
        evidence_stats = db.query(
            EvidenceFeedback.evidence_id,
            EvidenceFeedback.evidence_type,
            func.count(EvidenceFeedback.id).label("total_feedback"),
            func.avg(EvidenceFeedback.feedback_score).label("avg_score"),
        ).group_by(
            EvidenceFeedback.evidence_id,
            EvidenceFeedback.evidence_type
        ).having(func.count(EvidenceFeedback.id) >= 3).all()
        
        # Sort by score
        top_evidence = sorted(evidence_stats, key=lambda x: x.avg_score, reverse=True)[:10]
        bottom_evidence = sorted(evidence_stats, key=lambda x: x.avg_score)[:10]
        
        return {
            "by_type": [
                {
                    "evidence_type": row.evidence_type,
                    "total_feedback": row.total_feedback,
                    "avg_score": round(float(row.avg_score), 3) if row.avg_score else 0.0,
                    "upvotes": row.upvotes,
                    "downvotes": row.downvotes,
                    "flags": row.flags,
                }
                for row in type_stats
            ],
            "top_evidence": [
                {
                    "evidence_id": row.evidence_id,
                    "evidence_type": row.evidence_type,
                    "total_feedback": row.total_feedback,
                    "avg_score": round(float(row.avg_score), 3) if row.avg_score else 0.0,
                }
                for row in top_evidence
            ],
            "bottom_evidence": [
                {
                    "evidence_id": row.evidence_id,
                    "evidence_type": row.evidence_type,
                    "total_feedback": row.total_feedback,
                    "avg_score": round(float(row.avg_score), 3) if row.avg_score else 0.0,
                }
                for row in bottom_evidence
            ],
        }
        
    finally:
        db.close()


@router.delete("/{feedback_id}")
def delete_feedback(feedback_id: str):
    """Delete a feedback entry (admin only in production)."""
    db = SessionLocal()
    try:
        feedback = db.query(EvidenceFeedback).filter(EvidenceFeedback.id == feedback_id).first()
        if not feedback:
            raise HTTPException(status_code=404, detail="Feedback not found")
        
        db.delete(feedback)
        db.commit()
        
        return {"status": "deleted", "feedback_id": feedback_id}
        
    finally:
        db.close()
