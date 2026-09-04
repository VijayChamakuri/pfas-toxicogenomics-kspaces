from __future__ import annotations

import hashlib
import json

import pandas as pd

from .catalog import CHEMICALS, ArtifactCatalog
from .models import (
    Analysis,
    ApprovedExecution,
    ExecutionResult,
    FindingSeverity,
    PlanStatus,
    ProvenanceRecord,
    ReviewDecision,
    ValidationFinding,
    WorkflowPlan,
)


class WorkflowExecutor:
    def __init__(self, catalog: ArtifactCatalog):
        self.catalog = catalog

    def execute(self, approved: ApprovedExecution) -> ExecutionResult:
        plan, review = approved.plan, approved.review
        if plan.status != PlanStatus.READY_FOR_REVIEW or plan.analysis is None:
            raise ValueError("Only review-ready plans can be executed")
        if review.decision != ReviewDecision.APPROVED:
            raise ValueError("Human review did not approve this workflow")
        digest = plan.review_digest()
        if review.plan_digest != digest:
            raise ValueError("Approval does not match the current workflow plan")
        evidence = self._validate_evidence(plan)
        records, result_files = self._run(plan)
        provenance = {
            str(path.relative_to(self.catalog.root)): ProvenanceRecord(
                sha256=self.catalog.sha256(path), roles=["result_input"]
            )
            for path in result_files
        }
        for source, record in evidence.items():
            if source in provenance:
                provenance[source].roles.extend(
                    role for role in record.roles if role not in provenance[source].roles
                )
            else:
                provenance[source] = record
        run_payload = json.dumps(
            {"plan_digest": digest, "provenance": {key: value.sha256 for key, value in provenance.items()}},
            sort_keys=True,
        ).encode()
        execution_id = hashlib.sha256(run_payload).hexdigest()[:20]
        findings = plan.findings + [ValidationFinding(
            code="expert_interpretation", severity=FindingSeverity.INFO,
            message="The approved computation completed; biological interpretation remains subject to domain-expert review.",
        )]
        return ExecutionResult(
            execution_id=execution_id, plan_digest=digest, analysis=plan.analysis,
            records=records, provenance=provenance, findings=findings, review=review,
        )

    def _validate_evidence(self, plan: WorkflowPlan) -> dict[str, ProvenanceRecord]:
        provenance: dict[str, ProvenanceRecord] = {}
        for citation in plan.citations:
            path = (self.catalog.root / citation.source).resolve()
            if not path.is_relative_to(self.catalog.root):
                raise ValueError("Evidence path escapes the project root")
            if not path.is_file() or self.catalog.sha256(path) != citation.source_sha256:
                raise ValueError(f"Evidence changed after planning: {citation.source}")
            provenance[citation.source] = ProvenanceRecord(
                sha256=citation.source_sha256, roles=["planning_evidence"]
            )
        return provenance

    def _run(self, plan: WorkflowPlan):
        if plan.analysis == Analysis.DEG_SUMMARY:
            records, files = [], []
            threshold = float(plan.parameters["fdr_threshold"])
            for chemical in plan.chemicals or list(CHEMICALS):
                path = self.catalog.dge_path(chemical)
                frame = self.catalog.load_dge(path)
                records.append({
                    "chemical": chemical, "genes_tested": len(frame),
                    "fdr_threshold": threshold,
                    "genes_detected_at_threshold": int((frame.FDR < threshold).sum()),
                })
                files.append(path)
            return records, files
        if plan.analysis == Analysis.CONSENSUS_MODULE_LOOKUP:
            assignment_path = self.catalog.run_dir / "primary_k5_consensus_module_assignments.csv"
            strength_path = self.catalog.run_dir / "primary_k5_consensus_strength.csv"
            merged = pd.read_csv(assignment_path).merge(pd.read_csv(strength_path), on="WB_id", how="left")
            selected = merged[merged.WB_id.str.upper().isin({gene.upper() for gene in plan.gene_ids})]
            return selected.to_dict("records"), [assignment_path, strength_path]
        if plan.analysis == Analysis.GO_LOOKUP:
            path = self.catalog.run_dir / "go_enrichment" / "consensus_k5_GO_enrichment.csv"
            frame = pd.read_csv(path).sort_values("p_fdr_bh")
            top_n = int(plan.parameters["top_terms_per_module"])
            return frame.groupby("module", as_index=False).head(top_n).to_dict("records"), [path]
        if plan.analysis == Analysis.EXPLORATORY_PAIR_LOOKUP:
            path = self.catalog.run_dir / "residual_metric_cross_k_robust_pairs.csv"
            frame = pd.read_csv(path)
            wanted = set(plan.chemicals)
            frame = frame[frame.pair.apply(lambda value: set(value.split("-")) == wanted)]
            return frame.to_dict("records"), [path]
        raise ValueError(f"Unsupported analysis: {plan.analysis}")
