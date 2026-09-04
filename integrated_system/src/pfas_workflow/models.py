from __future__ import annotations

import hashlib
import json
from enum import StrEnum
from typing import Any

from pydantic import BaseModel, Field, model_validator


class Analysis(StrEnum):
    DEG_SUMMARY = "deg_summary"
    CONSENSUS_MODULE_LOOKUP = "consensus_module_lookup"
    GO_LOOKUP = "go_lookup"
    EXPLORATORY_PAIR_LOOKUP = "exploratory_pair_lookup"


class PlanStatus(StrEnum):
    READY_FOR_REVIEW = "ready_for_review"
    NEEDS_CLARIFICATION = "needs_clarification"
    INCOMPATIBLE = "incompatible"
    ABSTAINED = "abstained"


class FindingSeverity(StrEnum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class ReviewDecision(StrEnum):
    APPROVED = "approved"
    REJECTED = "rejected"


class BiologicalRequest(BaseModel):
    question: str = Field(min_length=8, max_length=2_000)


class InterpretedRequest(BaseModel):
    intent: Analysis | None = None
    chemicals: list[str] = Field(default_factory=list)
    gene_ids: list[str] = Field(default_factory=list)
    requested_parameters: dict[str, Any] = Field(default_factory=dict)


class EvidenceCitation(BaseModel):
    evidence_id: str
    source: str
    source_sha256: str
    supports: str


class ValidationFinding(BaseModel):
    code: str
    severity: FindingSeverity
    message: str


class WorkflowPlan(BaseModel):
    schema_version: str = "2.0"
    status: PlanStatus
    analysis: Analysis | None = None
    chemicals: list[str] = Field(default_factory=list)
    gene_ids: list[str] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    steps: list[str] = Field(default_factory=list)
    reason: str
    citations: list[EvidenceCitation] = Field(default_factory=list)
    findings: list[ValidationFinding] = Field(default_factory=list)
    requires_human_review: bool = True
    planner_version: str = "evidence-rules-2.0"
    router_metadata: dict[str, Any] = Field(default_factory=dict)

    @property
    def accepted(self) -> bool:
        return self.status == PlanStatus.READY_FOR_REVIEW

    @model_validator(mode="after")
    def validate_lifecycle_state(self):
        ready = self.status == PlanStatus.READY_FOR_REVIEW
        if ready != (self.analysis is not None):
            raise ValueError("review-ready plans require an analysis; other states must not have one")
        if ready and not self.citations:
            raise ValueError("review-ready plans require verified evidence citations")
        if ready and any(item.severity == FindingSeverity.ERROR for item in self.findings):
            raise ValueError("plans with validation errors cannot be review-ready")
        return self

    def review_digest(self) -> str:
        payload = self.model_dump(mode="json")
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(encoded).hexdigest()


class HumanReview(BaseModel):
    decision: ReviewDecision
    reviewer: str = Field(min_length=2, max_length=200)
    plan_digest: str = Field(pattern=r"^[0-9a-f]{64}$")
    comment: str = Field(default="", max_length=2_000)


class ApprovedExecution(BaseModel):
    plan: WorkflowPlan
    review: HumanReview


class ProvenanceRecord(BaseModel):
    sha256: str
    roles: list[str]


class ExecutionResult(BaseModel):
    execution_id: str
    plan_digest: str
    analysis: Analysis
    records: list[dict[str, Any]]
    provenance: dict[str, ProvenanceRecord]
    findings: list[ValidationFinding]
    review: HumanReview
    executor_version: str = "artifact-executor-2.0"
