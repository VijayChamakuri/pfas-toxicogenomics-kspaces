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
    # Study-design limits, distinct from missing modalities.
    "raw count": "Replicate-level count matrices are not distributed with this study.",
    "count matrix": "Replicate-level count matrices are not distributed with this study.",
    "dose-response": "Each chemical was assayed at a single experiment-specific EC50; no dose series exists.",
    "dose response": "Each chemical was assayed at a single experiment-specific EC50; no dose series exists.",
    "across concentrations": "Each chemical was assayed at a single experiment-specific EC50; no dose series exists.",
    "time course": "The study has no time-course design.",
    "time-course": "The study has no time-course design.",
    "per-replicate": "Replicate-level data are not distributed with this study.",
    "replicate variance": "Replicate-level data are not distributed with this study.",
}

# Causal and mechanistic language. Word-boundary alternation, not substring membership:
# the guard's failure mode is a fabricated causal claim reaching a scientist, so it matches
# the ordinary ways a biologist phrases causation rather than five literal strings.
CAUSAL_PATTERN = re.compile(
    r"""\b(?:caus(?:e|es|ed|ing|al|ally|ation))\b
      | \b(?:prove|proves|proved|proven|proof)\b
      | \b(?:guarantee|guarantees|guaranteed)\b
      | \bmechanism\b
      | \bresponsible\s+for\b
      | \blead(?:s|ing)?\s+to\b
      | \bdriv(?:e|es|en|ing)\b
      | \bbecause\s+of\b
      | \bdue\s+to\b
      | \b(?:effect|impact|influence)\s+(?:of|does|do|on)\b
      | \bresult(?:s|ing)?\s+(?:in|from)\b
      | \bharmful\b
      | \btoxic\s+to\b
      | \bhuman\s+health\b
      | \bhumans?\b
    """,
    re.VERBOSE | re.IGNORECASE,
)

# Regulatory grouping is only a claim when it is applied to the chemicals or to policy.
REGULATORY_PATTERN = re.compile(
    r"\bregulat\w*\b(?=.*\b(?:pfas|chemical|chemicals|class|policy|risk|limit|limits)\b)",
    re.IGNORECASE | re.DOTALL,
)

# Withdrawn or never-reconstructed output, in user vocabulary as well as internal vocabulary.
RETRACTED_PATTERN = re.compile(
    r"""\bretracted\b
      | \bpercent-per-module\b
      | \bunadjusted\b
      | \breactome\b
      | \bglutathione\s+transferase\b
      | \b(?:earlier|previous|original|prior)\s+(?:analysis|finding|findings|ranking|rankings|result|results|screen)\b
      | \btop\s+eleven\b
    """,
    re.VERBOSE | re.IGNORECASE,
)

CAUSAL_TERMS = ("prove", "causes", "causal", "guarantees", "definitive mechanism")
UNSUPPORTED_TOOLS = ("seurat", "scanpy")
UNSUPPORTED_COLUMNS = ("age", "sex", "patient_id", "batch_id", "survival_time")
# Word-boundary forms, so "average", "percentage" and "coverage" are not read as "age".
UNSUPPORTED_COLUMN_PATTERNS = {
    column: re.compile(rf"\b{re.escape(column)}\b", re.IGNORECASE) for column in UNSUPPORTED_COLUMNS
}

# "GO" is an acronym: match it case-sensitively so that "genes go up" is not an ontology request.
GO_ACRONYM = re.compile(r"\bGO\b")
GO_VOCABULARY = (
    "gene ontology", "ontology", "pathway", "pathways",
    "enriched", "enrichment", "biological process", "biological processes",
    "molecular function", "molecular functions", "cellular component", "cellular components",
)

# Threshold parsing accepts decimal and percent forms; an unparsed but clearly intended
# threshold must never fall back to the default.
THRESHOLD_PATTERN = re.compile(
    r"(?:fdr|alpha|q[-\s]?value)"
    r"(?:\s*(?:<=|<|=|of|at|to|is|below|under|cutoff|threshold|value|level))*"
    r"\s*(\d*\.?\d+)\s*(%?)",
    re.IGNORECASE,
)
THRESHOLD_MENTION = re.compile(r"\b(?:fdr|alpha|q[-\s]?value)\b", re.IGNORECASE)

COUNT_WORDS = ("how many", "number of", "count", "counts")
SIMILARITY_WORDS = ("similar", "similarity", "shared", "pair", "pairs")
# Non-WormBase gene symbols, e.g. TP53 or BRCA1.
GENE_SYMBOL_PATTERN = re.compile(r"\b[A-Z]{2,6}\d{1,3}\b")

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
        raw = request.question
        text = raw.lower()
        chemicals = [chemical for chemical in CHEMICALS if chemical.lower() in text]
        genes = sorted(set(re.findall(r"WBGene\d{8}", raw, flags=re.IGNORECASE)))
        parameters: dict = {"fdr_threshold": self.parse_threshold(raw)}

        wants_go = bool(GO_ACRONYM.search(raw)) or any(word in text for word in GO_VOCABULARY)
        wants_counts = any(word in text for word in COUNT_WORDS)
        wants_similarity = any(word in text for word in SIMILARITY_WORDS) or "compare" in text
        mentions_module = "module" in text

        intent: Analysis | None = None
        if genes:
            intent, parameters = Analysis.CONSENSUS_MODULE_LOOKUP, {}
        elif mentions_module and wants_counts and not wants_go:
            # "How many genes are in each consensus module?" is a membership question.
            intent, parameters = Analysis.CONSENSUS_MODULE_LOOKUP, {}
        elif wants_go:
            intent, parameters = Analysis.GO_LOOKUP, {"top_terms_per_module": 5}
        elif wants_counts and (chemicals or "differential" in text or "deg" in text):
            # Counting DEGs for two chemicals is a summary, not a similarity test.
            intent = Analysis.DEG_SUMMARY
        elif (len(chemicals) >= 2 or "pfas" in text or "chemical" in text) and wants_similarity:
            intent, parameters = Analysis.EXPLORATORY_PAIR_LOOKUP, {}
        elif chemicals or any(
            word in text
            for word in ("deg", "differential", "summary", "response", "measured", "universe", "contrast")
        ):
            intent = Analysis.DEG_SUMMARY
        return InterpretedRequest(intent=intent, chemicals=chemicals, gene_ids=genes, requested_parameters=parameters)

    @staticmethod
    def parse_threshold(raw: str) -> float | None:
        """Return the requested threshold, or None when one was asked for but not parseable.

        A silently substituted default answers a different question than the one asked,
        so an unparsed mention is surfaced as an error rather than replaced.
        """
        match = THRESHOLD_PATTERN.search(raw)
        if match:
            value = float(match.group(1))
            return value / 100 if match.group(2) == "%" else value
        if THRESHOLD_MENTION.search(raw):
            return None
        return 0.05

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
        unavailable_columns = [
            column for column, pattern in UNSUPPORTED_COLUMN_PATTERNS.items() if pattern.search(request.question)
        ]
        if unavailable_columns:
            return PlanStatus.INCOMPATIBLE, "Unavailable column(s): " + ", ".join(unavailable_columns) + ".", "unknown_column"
        if CAUSAL_PATTERN.search(request.question) or REGULATORY_PATTERN.search(request.question):
            return PlanStatus.ABSTAINED, "These observational analyses cannot establish a causal mechanism.", "unsupported_claim"
        if RETRACTED_PATTERN.search(request.question) or ("module" in text and "correlation" in text):
            return PlanStatus.ABSTAINED, RETRACTED_WARNING, "retracted_method"
        if not interpreted.gene_ids and ("module" in text or "assignment" in text):
            known_tokens = {chemical.upper() for chemical in CHEMICALS}
            symbols = [
                token for token in GENE_SYMBOL_PATTERN.findall(request.question)
                if token.upper() not in known_tokens
            ]
            if symbols:
                return (
                    PlanStatus.INCOMPATIBLE,
                    "Gene identifier(s) are absent from the analyzed universe: " + ", ".join(symbols) + ".",
                    "unknown_gene",
                )
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
            if threshold is None:
                findings.append(ValidationFinding(
                    code="invalid_fdr_threshold", severity=FindingSeverity.ERROR,
                    message="A threshold was requested but could not be parsed; refusing to substitute the default."))
                return findings
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
