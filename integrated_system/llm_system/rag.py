from __future__ import annotations

import math
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field

from pfas_workflow.models import BiologicalRequest, PlanStatus
from pfas_workflow.planner import RuleBasedPlanner

from .corpus import Chunk, build_corpus

CHEMICALS = ("PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA")
TOKEN_RE = re.compile(r"[a-z0-9][a-z0-9_-]+")
INCOMPATIBLE = {
    "fastq": "Raw sequencing reads are unavailable.",
    "raw reads": "Raw sequencing reads are unavailable.",
    "single-cell": "The dataset is bulk RNA-seq, not single-cell RNA-seq.",
    "proteomics": "No proteomics measurements are available.",
    "metabolomics": "No metabolomics measurements are available.",
    "survival": "No survival endpoint is available.",
    "human patients": "The study organism is C. elegans, not humans.",
}
UNSUPPORTED_COLUMNS = ("age", "sex", "patient_id", "batch_id", "survival_time")
CAUSAL_TERMS = ("prove", "causes", "causal", "guarantees", "definitive mechanism")


@dataclass(frozen=True)
class Citation:
    chunk_id: str
    source: str
    heading: str
    sha256: str


@dataclass
class WorkflowProposal:
    schema_version: str = "1.0"
    status: str = "abstain"
    intent: str = "unknown"
    workflow_steps: list[str] = field(default_factory=list)
    required_inputs: list[str] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    citations: list[Citation] = field(default_factory=list)
    explanation: str = ""
    requires_expert_review: bool = True

    def as_dict(self) -> dict:
        result = asdict(self)
        return result


def _tokens(text: str) -> set[str]:
    stop = {"the", "and", "for", "with", "from", "that", "this", "into", "are", "was", "what"}
    return {token for token in TOKEN_RE.findall(text.lower()) if token not in stop}


class LexicalRetriever:
    """Transparent deterministic retrieval suitable for a small reviewed corpus."""

    def __init__(self, chunks: Iterable[Chunk] | None = None):
        self.chunks = list(chunks or build_corpus())
        self.document_frequency: dict[str, int] = {}
        for chunk in self.chunks:
            for token in _tokens(chunk.text + " " + chunk.heading):
                self.document_frequency[token] = self.document_frequency.get(token, 0) + 1

    def retrieve(self, query: str, k: int = 4, minimum_score: float = 0.08) -> list[tuple[Chunk, float]]:
        query_tokens = _tokens(query)
        ranked = []
        for chunk in self.chunks:
            chunk_tokens = _tokens(chunk.heading + " " + chunk.text)
            overlap = query_tokens & chunk_tokens
            score = sum(math.log((len(self.chunks) + 1) / (self.document_frequency[t] + 1)) + 1 for t in overlap)
            score /= math.sqrt(max(len(query_tokens), 1) * max(len(chunk_tokens), 1))
            if score >= minimum_score:
                ranked.append((chunk, score))
        return sorted(ranked, key=lambda item: (-item[1], item[0].chunk_id))[:k]


class BaseSystem:
    def _screen(self, request: str) -> WorkflowProposal | None:
        text = request.lower()
        allowed = {c.lower() for c in CHEMICALS} | {"pfas"}
        unknown = sorted(set(re.findall(r"\bpf[a-z0-9]+\b", text)) - allowed)
        if unknown:
            return WorkflowProposal(status="incompatible", explanation=f"Unsupported chemical identifier(s): {', '.join(unknown)}.")
        for phrase, reason in INCOMPATIBLE.items():
            if phrase in text:
                return WorkflowProposal(status="incompatible", explanation=reason)
        traps = [column for column in UNSUPPORTED_COLUMNS if column in text]
        if traps:
            return WorkflowProposal(status="incompatible", explanation=f"Unavailable column(s): {', '.join(traps)}.")
        if any(term in text for term in CAUSAL_TERMS):
            return WorkflowProposal(status="abstain", explanation="Observational expression and phenotype analyses cannot establish a causal mechanism.")
        if "percent-per-module" in text or ("module" in text and "correlation" in text):
            return WorkflowProposal(status="abstain", explanation="The raw module-profile correlation metric was retracted after null stress-testing.")
        if len(request.strip()) < 12 or request.lower().strip() in {"analyze biology", "analyze the data"}:
            return WorkflowProposal(status="needs_clarification", explanation="Specify a chemical, gene, endpoint, or analysis question.")
        return None

    @staticmethod
    def _intent(request: str) -> str:
        text = request.lower()
        if re.search(r"wbgene\d{8}", text):
            return "consensus_module_lookup"
        if "go" in text or "ontology" in text or "pathway" in text:
            return "go_enrichment_review"
        if "compare" in text or "shared" in text or "similar" in text:
            return "cross_chemical_comparison"
        if "deg" in text or "differential" in text or any(c.lower() in text for c in CHEMICALS):
            return "differential_expression_summary"
        return "unknown"

    @staticmethod
    def _steps(intent: str) -> tuple[list[str], list[str]]:
        mapping = {
            "consensus_module_lookup": (["Validate WormBase gene identifiers", "Join consensus module assignments and strength", "Report membership with uncertainty"], ["WormBase gene identifier", "Consensus assignment artifact"]),
            "go_enrichment_review": (["Select the corrected gene set", "Use the namespace-mapped eligible population", "Apply multiplicity correction across the declared test family", "Report terms as exploratory associations"], ["Gene set", "Mapped measured-gene population", "GO annotation version"]),
            "cross_chemical_comparison": (["Load harmonized logFC contrasts", "Apply a prespecified corrected similarity statistic", "Evaluate against an explicit null", "Report uncertainty and sensitivity across k"], ["Precomputed edgeR contrasts", "Prespecified statistic and null"]),
            "differential_expression_summary": (["Validate the shared 13,852-gene universe", "Apply per-contrast FDR threshold", "Summarize effect direction and magnitude", "State cross-contrast multiplicity limitation"], ["Precomputed edgeR DGE CSV files"]),
        }
        return mapping.get(intent, ([], []))


class PromptOnlySystem(BaseSystem):
    """Deterministic no-retrieval baseline with no source citations."""

    def propose(self, request: str) -> WorkflowProposal:
        screened = self._screen(request)
        if screened:
            return screened
        intent = self._intent(request)
        if intent == "unknown":
            return WorkflowProposal(status="needs_clarification", explanation="The request does not map to a supported workflow.")
        steps, inputs = self._steps(intent)
        return WorkflowProposal(status="accepted", intent=intent, workflow_steps=steps, required_inputs=inputs,
            assumptions=["The request concerns the existing C. elegans PFAS study."],
            warnings=["This proposal requires scientific expert review."],
            explanation="Generated from fixed workflow rules without retrieval.")


class GroundedWorkflowSystem(BaseSystem):
    """Evaluation adapter over the same planner used by the API.

    The earlier implementation duplicated request interpretation, retrieval,
    and safety rules inside this module. Keeping a second planner made an
    evaluation pass poor evidence for API behavior. This adapter preserves the
    stable evaluation schema while delegating every grounded decision to the
    production planner.
    """

    def __init__(self, retriever: LexicalRetriever | None = None):
        if retriever is not None:
            raise ValueError("custom retrievers are not supported by the production planner adapter")
        self.planner = RuleBasedPlanner()

    def propose(self, request: str) -> WorkflowProposal:
        plan = self.planner.plan(BiologicalRequest(question=request))
        statuses = {
            PlanStatus.READY_FOR_REVIEW: "accepted",
            PlanStatus.NEEDS_CLARIFICATION: "needs_clarification",
            PlanStatus.INCOMPATIBLE: "incompatible",
            PlanStatus.ABSTAINED: "abstain",
        }
        intents = {
            "deg_summary": "differential_expression_summary",
            "consensus_module_lookup": "consensus_module_lookup",
            "go_lookup": "go_enrichment_review",
            "exploratory_pair_lookup": "cross_chemical_comparison",
        }
        citations = [
            Citation(
                chunk_id=item.evidence_id,
                source=item.source,
                heading=item.supports,
                sha256=item.source_sha256,
            )
            for item in plan.citations
        ]
        return WorkflowProposal(
            status=statuses[plan.status],
            intent=intents.get(str(plan.analysis), "unknown"),
            workflow_steps=plan.steps,
            required_inputs=[],
            assumptions=["The production artifact catalog and evidence store were validated."],
            warnings=[item.message for item in plan.findings],
            citations=citations,
            explanation=plan.reason,
            requires_expert_review=plan.requires_human_review,
        )
