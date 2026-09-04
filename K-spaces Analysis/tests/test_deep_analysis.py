from pathlib import Path

import numpy as np
import pandas as pd

from pfas_master.deep_analysis import (
    adjust_enrichment_family,
    chemical_module_preference,
    detection_matrix,
    exact_intersection_table,
    module_activity_pca,
    pairwise_detection_overlap,
    unique_shared_gene_sets,
)
from pfas_master.science import CHEMICALS, load_and_validate_dge

ROOT = Path(__file__).parents[2]


def test_detection_intersections_conserve_detected_union():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    detected = detection_matrix(frames)
    intersections = exact_intersection_table(detected)
    assert intersections.gene_count.sum() == detected.any(axis=1).sum()
    assert intersections.degree.between(1, 10).all()
    assert not intersections.intersection.duplicated().any()


def test_unique_and_shared_sets_partition_each_chemical():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    detected = detection_matrix(frames)
    unique, shared, summary = unique_shared_gene_sets(detected)
    for chemical in CHEMICALS:
        assert unique[chemical].isdisjoint(shared[chemical])
        assert len(unique[chemical] | shared[chemical]) == int(detected[chemical].sum())
    assert (summary.detected_total == summary.detected_exclusively_here + summary.detected_and_shared).all()


def test_pairwise_overlap_covers_all_pairs_and_controls_family():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    detected = detection_matrix(frames)
    result = pairwise_detection_overlap(frames, detected)
    assert len(result) == 45
    assert result.q_bh_45_pairs.between(0, 1).all()
    assert result.concordant_sign_pct.between(0, 100).all()


def test_chemical_module_preference_covers_declared_50_cell_family():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    detected = detection_matrix(frames)
    genes = detected.index[detected.any(axis=1)]
    assignments = pd.Series(np.arange(len(genes)) % 5, index=genes)
    result = chemical_module_preference(detected, assignments)
    assert len(result) == 50
    assert result.q_bh_50_cells.between(0, 1).all()


def test_module_activity_pca_uses_conditions_as_observations():
    matrix = pd.DataFrame(
        np.arange(120, dtype=float).reshape(12, 10),
        index=[f"g{i}" for i in range(12)],
        columns=CHEMICALS,
    )
    assignments = pd.Series(np.arange(12) % 3, index=matrix.index)
    activity, pca = module_activity_pca(matrix, assignments)
    assert activity.shape == (3, 10)
    assert pca.shape == (10, 4)
    assert pca.index.tolist() == list(CHEMICALS)


def test_enrichment_adjustment_is_one_family():
    results = pd.DataFrame(
        {"p_raw": [0.001, 0.02, 0.5], "enrichment": ["e", "e", "p"]}
    )
    adjusted = adjust_enrichment_family(results)
    assert adjusted.p_fdr_bh_family.is_monotonic_increasing
    assert adjusted.p_fdr_bh_family.between(0, 1).all()
