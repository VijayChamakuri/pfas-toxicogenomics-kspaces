"""Tested helpers for the merged PFAS analysis notebook."""

from .science import (
    CHEMICALS,
    build_detection_evidence,
    build_legacy_zero_masked_matrix,
    build_module_contributions,
    load_and_validate_dge,
    true_subsample_stability,
)
from .deep_analysis import (
    adjust_enrichment_family,
    chemical_module_preference,
    detection_matrix,
    exact_intersection_table,
    module_activity_pca,
    pairwise_detection_overlap,
    unique_shared_gene_sets,
)

__all__ = [
    "CHEMICALS",
    "build_detection_evidence",
    "build_legacy_zero_masked_matrix",
    "build_module_contributions",
    "load_and_validate_dge",
    "true_subsample_stability",
    "adjust_enrichment_family",
    "chemical_module_preference",
    "detection_matrix",
    "exact_intersection_table",
    "module_activity_pca",
    "pairwise_detection_overlap",
    "unique_shared_gene_sets",
]
