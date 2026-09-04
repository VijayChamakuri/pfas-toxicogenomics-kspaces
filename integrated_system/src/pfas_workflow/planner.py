from __future__ import annotations

import re
from pathlib import Path
from typing import Protocol

import pandas as pd

from .catalog import CHEMICALS, ArtifactCatalog
from .evidence import VerifiedEvidenceStore
from .modeling import router_disclosure
from .models import (
    Analysis,
    BiologicalRequest,
    FindingSeverity,
    InterpretedRequest,
    PlanStatus,
    ValidationFinding,
    WorkflowPlan,
)

RETRACTED_WARNING = (
    "Raw percent-per-module chemical correlations were retracted because their null median was "
    "about 0.98. Pair output is exploratory and uses the corrected residual analysis only."
)
INCOMPATIBLE = {
    "fastq": "Raw sequencing reads are unavailable.",
    "raw reads": "Raw sequencing reads are unavailable.",
    "deseq2": "Raw counts and sample metadata required by DESeq2 are unavailable.",
    "batch correction": "Sample-level batch metadata is unavailable.",
    "single-cell": "The available expression data are bulk RNA-seq contrasts.",
    "proteomics": "No proteomics measurements are available.",
    "metabolomics": "No metabolomics measurements are available.",
    "survival": "No survival endpoint is available.",
}
CAUSAL_TERMS = ("prove", "causes", "causal", "guarantees", "definitive mechanism")
UNSUPPORTED_TOOLS = ("seurat", "scanpy")
UNSUPPORTED_COLUMNS = ("age", "sex", "patient_id", "batch_id", "survival_time")
STEPS = {
    Analysis.DEG_SUMMARY: ["Validate the shared measured-gene universe", "Apply the declared per-contrast FDR threshold", "Summarize detection counts while preserving effect-size context"],
    Analysis.CONSENSUS_MODULE_LOOKUP: ["Validate WormBase identifiers against the consensus artifact", "Join majority-voted assignments with consensus strength", "Report module membership with uncertainty"],
    Analysis.GO_LOOKUP: ["Load the stored consensus-module GO results", "Validate adjusted probabilities and declared test family", "Report terms as exploratory associations"],
    Analysis.EXPLORATORY_PAIR_LOOKUP: ["Load the corrected residual-similarity sensitivity results", "Filter to the requested chemical pair", "Report multiplicity and clustering-uncertainty boundaries"],
}


class Planner(Protocol):
    def plan(self, request: BiologicalRequest) -> WorkflowPlan: ...


class RuleBasedPlanner:
    """Evidence-grounded deterministic planner and execution safety gate."""

    def __init__(self, catalog: ArtifactCatalog | None = None):
        root = Path(__file__).resolve().parents[3]
        self.catalog = catalog or ArtifactCatalog(root)
        self.evidence = VerifiedEvidenceStore(self.catalog)
        metrics = root / "integrated_system" / "llm_system" / "training_output" / "metrics.json"
        self.router_metadata = router_disclosure(metrics)

    def plan(self, request: BiologicalRequest) -> WorkflowPlan:
        interpreted = self.interpret(request)
        compatibility = self.check_compatibility(request, interpreted)
        if compatibility:
            return self._stop(*compatibility)
        if interpreted.intent is None:
            return self._stop(PlanStatus.NEEDS_CLARIFICATION, "The question does not map unambiguously to a supported analysis.", "ambiguous_request")
        findings = self.validate(interpreted)
        errors = [item for item in findings if item.severity == FindingSeverity.ERROR]
        if errors:
            return self._stop(PlanStatus.ABSTAINED, errors[0].message, errors[0].code, findings)
        try:
            citations = self.evidence.retrieve(interpreted.intent)
        except FileNotFoundError as exc:
            return self._stop(PlanStatus.ABSTAINED, str(exc), "missing_evidence")
        return WorkflowPlan(
            status=PlanStatus.READY_FOR_REVIEW, analysis=interpreted.intent,
            chemicals=interpreted.chemicals, gene_ids=interpreted.gene_ids,
            parameters=interpreted.requested_parameters, steps=STEPS[interpreted.intent],
            reason="The request is compatible, evidence-grounded, and passes deterministic scientific rules.",
            citations=citations, findings=findings, router_metadata=self.router_metadata,
        )

    def interpret(self, request: BiologicalRequest) -> InterpretedRequest:
        text = request.question.lower()
        chemicals = [chemical for chemical in CHEMICALS if chemical.lower() in text]
        genes = sorted(set(re.findall(r"WBGene\d{8}", request.question, flags=re.IGNORECASE)))
        alpha_match = re.search(r"(?:fdr|alpha)\s*(?:<|=|at)?\s*(0?\.\d+)", text)
        parameters = {"fdr_threshold": float(alpha_match.group(1)) if alpha_match else 0.05}
        intent: Analysis | None = None
        if genes:
            intent, parameters = Analysis.CONSENSUS_MODULE_LOOKUP, {}
        elif "go" in text or "ontology" in text or "pathway" in text:
            intent, parameters = Analysis.GO_LOOKUP, {"top_terms_per_module": 5}
        elif (len(chemicals) >= 2 or "pfas" in text) and any(word in text for word in ("compare", "similar", "shared", "pair")):
            intent, parameters = Analysis.EXPLORATORY_PAIR_LOOKUP, {}
        elif chemicals or any(word in text for word in ("deg", "differential", "summary", "response")):
            intent = Analysis.DEG_SUMMARY
        return InterpretedRequest(intent=intent, chemicals=chemicals, gene_ids=genes, requested_parameters=parameters)

    def check_compatibility(self, request: BiologicalRequest, interpreted: InterpretedRequest) -> tuple[PlanStatus, str, str] | None:
        text = request.question.lower()
        mentioned = set(re.findall(r"\bPF[A-Za-z0-9]+\b", request.question, flags=re.IGNORECASE))
        known = {chemical.lower() for chemical in CHEMICALS} | {"pfas"}
        unknown = sorted(value for value in mentioned if value.lower() not in known)
        if unknown:
            return PlanStatus.INCOMPATIBLE, "Unsupported chemical identifier(s): " + ", ".join(unknown) + ".", "unknown_chemical"
        for phrase, reason in INCOMPATIBLE.items():
            if phrase in text:
                return PlanStatus.INCOMPATIBLE, reason, "missing_required_data"
        if any(tool in text for tool in UNSUPPORTED_TOOLS):
            return PlanStatus.INCOMPATIBLE, "The requested tool is not supported by this dataset and workflow.", "unsupported_tool"
        unavailable_columns = [column for column in UNSUPPORTED_COLUMNS if column in text]
        if unavailable_columns:
            return PlanStatus.INCOMPATIBLE, "Unavailable column(s): " + ", ".join(unavailable_columns) + ".", "unknown_column"
        if any(term in text for term in CAUSAL_TERMS):
            return PlanStatus.ABSTAINED, "These observational analyses cannot establish a causal mechanism.", "unsupported_claim"
        if "retracted" in text or "percent-per-module" in text or ("module" in text and "correlation" in text):
            return PlanStatus.ABSTAINED, RETRACTED_WARNING, "retracted_method"
        if interpreted.intent == Analysis.CONSENSUS_MODULE_LOOKUP and interpreted.gene_ids:
            path = self.catalog.run_dir / "primary_k5_consensus_module_assignments.csv"
            if not path.is_file():
                return PlanStatus.ABSTAINED, "Consensus assignment evidence is unavailable.", "missing_evidence"
            known_genes = set(pd.read_csv(path, usecols=["WB_id"]).WB_id.str.upper())
            missing = [gene for gene in interpreted.gene_ids if gene.upper() not in known_genes]
            if missing:
                return PlanStatus.INCOMPATIBLE, "Gene identifier(s) are absent from the analyzed universe: " + ", ".join(missing) + ".", "unknown_gene"
        return None

    @staticmethod
    def validate(interpreted: InterpretedRequest) -> list[ValidationFinding]:
        findings: list[ValidationFinding] = []
        if interpreted.intent == Analysis.DEG_SUMMARY:
            threshold = interpreted.requested_parameters["fdr_threshold"]
            if not 0 < threshold <= 0.1:
                findings.append(ValidationFinding(code="invalid_fdr_threshold", severity=FindingSeverity.ERROR, message="The FDR threshold must be greater than 0 and no more than 0.1."))
            findings.append(ValidationFinding(code="contrast_scope", severity=FindingSeverity.WARNING, message="FDR is controlled within each contrast, not across all ten contrasts."))
        elif interpreted.intent == Analysis.GO_LOOKUP:
            findings.append(ValidationFinding(
                code="go_scope",
                severity=FindingSeverity.WARNING,
                message=(
                    "Stored GO results use one family-wide BH correction, but remain exploratory "
                    "because annotation and module uncertainty are not fully propagated."
                ),
            ))
        elif interpreted.intent == Analysis.EXPLORATORY_PAIR_LOOKUP:
            findings.extend([
                ValidationFinding(code="retracted_predecessor", severity=FindingSeverity.WARNING, message=RETRACTED_WARNING),
                ValidationFinding(code="pair_scope", severity=FindingSeverity.WARNING, message="No pair passes BH q<0.05 at k=5 or k=10; the stored cross-k result is a valid empty table."),
            ])
        elif interpreted.intent == Analysis.CONSENSUS_MODULE_LOOKUP:
            findings.append(ValidationFinding(code="module_scope", severity=FindingSeverity.WARNING, message="Module membership is provisional and does not establish mechanism."))
        return findings

    def _stop(self, status: PlanStatus, reason: str, code: str, findings: list[ValidationFinding] | None = None) -> WorkflowPlan:
        items = findings or [ValidationFinding(code=code, severity=FindingSeverity.ERROR, message=reason)]
        return WorkflowPlan(status=status, reason=reason, findings=items, router_metadata=self.router_metadata)
