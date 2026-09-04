from pathlib import Path

import numpy as np
import pandas as pd

from pfas_master.science import (
    CHEMICALS,
    build_detection_evidence,
    build_module_contributions,
    load_and_validate_dge,
    true_subsample_stability,
)

ROOT = Path(__file__).parents[2]


def test_real_dge_schema_and_universe():
    frames, manifest = load_and_validate_dge(ROOT / "Data" / "DEGs")
    assert set(frames) == set(CHEMICALS)
    assert manifest.rows.eq(13_852).all()
    assert all(len(frame) == 13_852 for frame in frames.values())


def test_detection_partition_and_thresholds():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    evidence, sensitivity = build_detection_evidence(frames)
    assert len(evidence) == 10 * 13_852
    assert sensitivity.fdr_threshold.tolist() == [0.01, 0.05, 0.10]
    multiplicity = evidence.drop_duplicates("WB_id").detection_multiplicity
    assert multiplicity.between(0, 10).all()


def test_true_stability_uses_subsample_labels():
    data = np.arange(200, dtype=float).reshape(100, 2)

    def seed_sensitive_fit(subset, seed):
        cutoff = np.quantile(subset[:, 0], 0.35 + 0.1 * seed)
        return (subset[:, 0] > cutoff).astype(int)

    result = true_subsample_stability(data, seed_sensitive_fit, seeds=[1, 2, 3], fraction=0.8)
    assert len(result) == 3
    assert result.shared_genes.min() > 0
    assert not result.ari.eq(1).all()


def test_module_contribution_keys_are_unique():
    frames, _ = load_and_validate_dge(ROOT / "Data" / "DEGs")
    genes = frames[CHEMICALS[0]].index[:100]
    assignments = pd.Series(np.arange(100) % 5, index=genes)
    strength = pd.Series(0.8, index=genes)
    gene_level, summary = build_module_contributions(frames, assignments, strength)
    assert not gene_level.duplicated(["WB_id", "chemical"]).any()
    assert not summary.duplicated(["chemical", "consensus_module"]).any()
    assert len(gene_level) == 1000
