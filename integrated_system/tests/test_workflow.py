import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from pfas_workflow.api import app, catalog
from pfas_workflow.executor import WorkflowExecutor
from pfas_workflow.models import (
    ApprovedExecution,
    BiologicalRequest,
    HumanReview,
    PlanStatus,
    ReviewDecision,
    WorkflowPlan,
)
from pfas_workflow.planner import RuleBasedPlanner

ROOT = Path(__file__).parents[1]


def approve(plan, reviewer="scientist"):
    return ApprovedExecution(plan=plan, review=HumanReview(decision=ReviewDecision.APPROVED, reviewer=reviewer, plan_digest=plan.review_digest()))


def test_catalog_integrity():
    assert catalog.validate() == {"chemicals": 10, "genes_per_contrast": 13_852}


@pytest.mark.parametrize("case", json.loads((ROOT / "evals" / "cases.json").read_text()))
def test_gold_plans(case):
    result = RuleBasedPlanner(catalog).plan(BiologicalRequest(question=case["question"]))
    assert result.accepted is case["accepted"]
    assert (result.analysis.value if result.analysis else None) == case["analysis"]


def test_plan_integrates_interpretation_evidence_validation_and_review_state():
    plan = RuleBasedPlanner(catalog).plan(BiologicalRequest(question="Summarize PFOS differential expression at FDR 0.05."))
    assert plan.status == PlanStatus.READY_FOR_REVIEW
    assert plan.chemicals == ["PFOS"]
    assert plan.parameters == {"fdr_threshold": 0.05}
    assert len(plan.steps) == 3
    assert {item.evidence_id for item in plan.citations} == {"dge-schema", "method-decisions"}
    assert all(len(item.source_sha256) == 64 for item in plan.citations)
    assert plan.router_metadata["selected_router"] == "deterministic evidence rules"
    assert plan.router_metadata["adapted_model"]["evaluation_status"].startswith("loss-only")


def test_execution_requires_matching_approval_and_is_provenanced():
    plan = RuleBasedPlanner(catalog).plan(BiologicalRequest(question="Summarize PFOS differential expression."))
    first = WorkflowExecutor(catalog).execute(approve(plan))
    second = WorkflowExecutor(catalog).execute(approve(plan))
    assert first == second
    assert first.records == [{"chemical": "PFOS", "genes_tested": 13_852, "fdr_threshold": 0.05, "genes_detected_at_threshold": 190}]
    assert first.plan_digest == plan.review_digest()
    assert {role for item in first.provenance.values() for role in item.roles} == {
        "planning_evidence", "result_input",
    }
    assert len(first.execution_id) == 20


def test_rejected_refused_and_tampered_plans_cannot_execute():
    executor = WorkflowExecutor(catalog)
    refused = RuleBasedPlanner(catalog).plan(BiologicalRequest(question="Run raw reads through DESeq2."))
    review = HumanReview(decision=ReviewDecision.APPROVED, reviewer="scientist", plan_digest=refused.review_digest())
    with pytest.raises(ValueError, match="Only review-ready"):
        executor.execute(ApprovedExecution(plan=refused, review=review))
    plan = RuleBasedPlanner(catalog).plan(BiologicalRequest(question="Summarize PFOS differential expression."))
    rejected = HumanReview(decision=ReviewDecision.REJECTED, reviewer="scientist", plan_digest=plan.review_digest())
    with pytest.raises(ValueError, match="did not approve"):
        executor.execute(ApprovedExecution(plan=plan, review=rejected))
    approval = approve(plan).review
    plan.parameters["fdr_threshold"] = 0.01
    with pytest.raises(ValueError, match="does not match"):
        executor.execute(ApprovedExecution(plan=plan, review=approval))


@pytest.mark.parametrize(("question", "status", "code"), [
    ("Compare PFHxS with PFOS.", PlanStatus.INCOMPATIBLE, "unknown_chemical"),
    ("Run single-cell analysis for PFOS.", PlanStatus.INCOMPATIBLE, "missing_required_data"),
    ("Prove PFOS causes receptor toxicity.", PlanStatus.ABSTAINED, "unsupported_claim"),
    ("Summarize PFOS at FDR 0.5.", PlanStatus.ABSTAINED, "invalid_fdr_threshold"),
    ("Find WBGene99999999 in the consensus modules.", PlanStatus.INCOMPATIBLE, "unknown_gene"),
    ("Please analyze the biology in depth.", PlanStatus.NEEDS_CLARIFICATION, "ambiguous_request"),
])
def test_scientific_guardrails(question, status, code):
    plan = RuleBasedPlanner(catalog).plan(BiologicalRequest(question=question))
    assert plan.status == status
    assert plan.analysis is None
    assert plan.findings[0].code == code


def test_api_contract_exposes_planning_review_execution_and_modeling_contract():
    client = TestClient(app)
    assert client.get("/health").status_code == 200
    contract = client.get("/modeling-contract")
    assert contract.status_code == 200
    assert contract.json()["supported_backends"] == ["pytorch", "tensorflow"]
    response = client.post("/plan", json={"question": "Summarize PFBA differential expression."})
    assert response.status_code == 200
    plan = response.json()
    assert plan["analysis"] == "deg_summary"
    parsed = WorkflowPlan.model_validate(plan)
    payload = {"plan": plan, "review": {"decision": "approved", "reviewer": "scientist", "plan_digest": parsed.review_digest(), "comment": ""}}
    executed = client.post("/execute", json=payload)
    assert executed.status_code == 200
    assert executed.json()["analysis"] == "deg_summary"
