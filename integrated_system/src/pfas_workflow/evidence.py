from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .catalog import ArtifactCatalog
from .models import Analysis, EvidenceCitation


@dataclass(frozen=True)
class EvidenceRecord:
    evidence_id: str
    source: Path
    supports: str
    analyses: tuple[Analysis, ...]


class VerifiedEvidenceStore:
    """Allowlisted evidence index over reviewed project artifacts."""

    def __init__(self, catalog: ArtifactCatalog):
        root = catalog.root
        run_dir = catalog.run_dir
        self.catalog = catalog
        self.records = (
            EvidenceRecord(
                "dge-schema", catalog.dge_path("PFOS"),
                "DGE schema, measured-gene universe, effect sizes, and per-contrast FDR values",
                (Analysis.DEG_SUMMARY,),
            ),
            EvidenceRecord(
                "consensus-assignment", run_dir / "primary_k5_consensus_module_assignments.csv",
                "majority-voted consensus module membership",
                (Analysis.CONSENSUS_MODULE_LOOKUP,),
            ),
            EvidenceRecord(
                "consensus-strength", run_dir / "primary_k5_consensus_strength.csv",
                "per-gene consensus assignment strength",
                (Analysis.CONSENSUS_MODULE_LOOKUP,),
            ),
            EvidenceRecord(
                "go-results", run_dir / "go_enrichment" / "consensus_k5_GO_enrichment.csv",
                "stored GO enrichment records and adjusted p-values",
                (Analysis.GO_LOOKUP,),
            ),
            EvidenceRecord(
                "corrected-similarity", run_dir / "residual_metric_cross_k_robust_pairs.csv",
                "corrected residual similarities retained by the cross-k sensitivity filter",
                (Analysis.EXPLORATORY_PAIR_LOOKUP,),
            ),
            EvidenceRecord(
                "method-decisions", root / "K-spaces Analysis" / "FINDINGS.md",
                "methodological corrections, uncertainty boundaries, and retracted claims",
                tuple(Analysis),
            ),
        )

    def retrieve(self, analysis: Analysis) -> list[EvidenceCitation]:
        selected = [record for record in self.records if analysis in record.analyses]
        missing = [record.source for record in selected if not record.source.is_file()]
        if missing:
            raise FileNotFoundError("Required verified evidence is missing: " + ", ".join(map(str, missing)))
        return [EvidenceCitation(
            evidence_id=record.evidence_id,
            source=str(record.source.relative_to(self.catalog.root)),
            source_sha256=self.catalog.sha256(record.source),
            supports=record.supports,
        ) for record in selected]
