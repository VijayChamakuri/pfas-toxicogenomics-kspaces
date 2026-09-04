from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

CHEMICALS = ("PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA")
REQUIRED_COLUMNS = ("WB_id", "logFC", "logCPM", "F", "PValue", "FDR")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_validate_dge(degs_dir: Path) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Load immutable DGE tables, normalize the ID column, and fail on schema drift."""
    frames: dict[str, pd.DataFrame] = {}
    manifest: list[dict] = []
    reference: set[str] | None = None
    for chemical in CHEMICALS:
        path = degs_dir / f"{chemical}vsControl_DGE_results.csv"
        frame = pd.read_csv(path)
        if "WB_id" not in frame and "Unnamed: 0" in frame:
            frame = frame.rename(columns={"Unnamed: 0": "WB_id"})
        missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
        if missing:
            raise ValueError(f"{path.name} missing columns: {missing}")
        frame = frame.loc[:, REQUIRED_COLUMNS].copy()
        if len(frame) != 13_852:
            raise ValueError(f"{path.name} has {len(frame)} rows, expected 13852")
        if frame.WB_id.isna().any() or frame.WB_id.duplicated().any():
            raise ValueError(f"{path.name} has missing or duplicated WB_id values")
        numeric = frame.loc[:, REQUIRED_COLUMNS[1:]].apply(pd.to_numeric, errors="raise")
        if not np.isfinite(numeric.to_numpy()).all():
            raise ValueError(f"{path.name} contains non-finite statistics")
        if not numeric.PValue.between(0, 1).all() or not numeric.FDR.between(0, 1).all():
            raise ValueError(f"{path.name} has probability values outside [0, 1]")
        universe = set(frame.WB_id)
        if reference is None:
            reference = universe
        elif universe != reference:
            raise ValueError(f"{path.name} does not share the measured gene universe")
        frames[chemical] = frame.set_index("WB_id")
        manifest.append(
            {
                "chemical": chemical,
                "path": str(path.resolve()),
                "rows": len(frame),
                "sha256": _sha256(path),
            }
        )
    return frames, pd.DataFrame(manifest)


def build_detection_evidence(
    frames: dict[str, pd.DataFrame], thresholds: Sequence[float] = (0.01, 0.05, 0.10)
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return gene-level continuous evidence and threshold sensitivity summaries."""
    genes = frames[CHEMICALS[0]].index
    rows: list[pd.DataFrame] = []
    for chemical in CHEMICALS:
        part = frames[chemical].loc[genes, ["logFC", "FDR"]].reset_index()
        part["chemical"] = chemical
        part["is_detected_fdr_0_05"] = part.FDR < 0.05
        rows.append(part)
    evidence = pd.concat(rows, ignore_index=True)
    multiplicity = evidence.groupby("WB_id").is_detected_fdr_0_05.sum().rename("detection_multiplicity")
    evidence = evidence.join(multiplicity, on="WB_id")
    evidence["detection_category"] = pd.cut(
        evidence.detection_multiplicity,
        bins=[-1, 0, 1, 9, 10],
        labels=["not_detected", "detected_exclusively_one", "detected_2_to_9", "detected_all_10"],
    ).astype(str)

    sensitivity = []
    for threshold in thresholds:
        sets = {chemical: set(frame.index[frame.FDR < threshold]) for chemical, frame in frames.items()}
        all_genes = set().union(*sets.values())
        counts = {gene: sum(gene in geneset for geneset in sets.values()) for gene in all_genes}
        sensitivity.append(
            {
                "fdr_threshold": threshold,
                "detected_union": len(all_genes),
                "detected_exclusively_one": sum(value == 1 for value in counts.values()),
                "detected_all_10": sum(value == 10 for value in counts.values()),
            }
        )
    return evidence, pd.DataFrame(sensitivity)


def build_legacy_zero_masked_matrix(
    frames: dict[str, pd.DataFrame], top_n: int = 3000
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reproduce the historical matrix for diagnosis, never for primary inference."""
    selected: set[str] = set()
    for chemical, frame in frames.items():
        significant = frame.loc[frame.FDR < 0.05].copy()
        selected.update(significant.logFC.abs().nlargest(min(top_n, len(significant))).index)
    genes = sorted(selected)
    matrix = pd.DataFrame(index=genes, columns=CHEMICALS, dtype=float)
    diagnostics = []
    for chemical in CHEMICALS:
        frame = frames[chemical].reindex(genes)
        masked = frame.logFC.where(frame.FDR < 0.05, 0.0)
        matrix[chemical] = masked
        diagnostics.append(
            {
                "chemical": chemical,
                "genes": len(genes),
                "forced_zero_fraction": float((frame.FDR >= 0.05).mean()),
                "nonzero_effects_forced_to_zero": int(((frame.FDR >= 0.05) & (frame.logFC != 0)).sum()),
                "zero_tie_size": int((masked == 0).sum()),
            }
        )
    ranked = matrix.rank(axis=0, method="average", pct=True)
    row_sd = ranked.std(axis=1, ddof=0).replace(0, np.nan)
    z_matrix = ranked.sub(ranked.mean(axis=1), axis=0).div(row_sd, axis=0).dropna()
    return z_matrix, pd.DataFrame(diagnostics)


def true_subsample_stability(
    data: np.ndarray,
    fit_predict: Callable[[np.ndarray, int], np.ndarray],
    seeds: Sequence[int],
    fraction: float = 0.8,
) -> pd.DataFrame:
    """Compare independently fitted subsamples on their shared observations."""
    if not 0.5 <= fraction < 1:
        raise ValueError("fraction must be in [0.5, 1)")
    rng = np.random.default_rng(0)
    fits: list[tuple[np.ndarray, np.ndarray, int]] = []
    size = int(len(data) * fraction)
    for seed in seeds:
        indices = np.sort(rng.choice(len(data), size=size, replace=False))
        labels = np.asarray(fit_predict(data[indices], seed))
        if labels.shape != (size,):
            raise ValueError("fit_predict returned labels with the wrong shape")
        fits.append((indices, labels, seed))
    rows = []
    for left in range(len(fits)):
        for right in range(left + 1, len(fits)):
            idx_a, labels_a, seed_a = fits[left]
            idx_b, labels_b, seed_b = fits[right]
            shared, pos_a, pos_b = np.intersect1d(idx_a, idx_b, return_indices=True)
            rows.append(
                {
                    "seed_a": seed_a,
                    "seed_b": seed_b,
                    "shared_genes": len(shared),
                    "ari": adjusted_rand_score(labels_a[pos_a], labels_b[pos_b]),
                }
            )
    return pd.DataFrame(rows)


def build_module_contributions(
    frames: dict[str, pd.DataFrame], assignments: pd.Series, strength: pd.Series
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build one tidy gene-level table and one chemical-by-module summary."""
    base = pd.DataFrame(
        {
            "WB_id": assignments.index,
            "consensus_module": assignments.astype(int).to_numpy(),
            "consensus_strength": strength.reindex(assignments.index).to_numpy(),
        }
    )
    gene_rows = []
    for chemical in CHEMICALS:
        stats = frames[chemical].reindex(assignments.index)
        part = base.copy()
        part["chemical"] = chemical
        part["logFC"] = stats.logFC.to_numpy()
        part["FDR"] = stats.FDR.to_numpy()
        part["is_detected"] = part.FDR < 0.05
        part["direction"] = np.select(
            [part.logFC > 0, part.logFC < 0], ["up", "down"], default="zero"
        )
        gene_rows.append(part)
    gene_level = pd.concat(gene_rows, ignore_index=True)
    summary = (
        gene_level.groupby(["chemical", "consensus_module"], as_index=False)
        .agg(
            module_genes=("WB_id", "size"),
            detected_genes=("is_detected", "sum"),
            up_detected=("direction", lambda s: int(((s == "up") & gene_level.loc[s.index, "is_detected"]).sum())),
            down_detected=("direction", lambda s: int(((s == "down") & gene_level.loc[s.index, "is_detected"]).sum())),
            mean_logFC=("logFC", "mean"),
            median_logFC=("logFC", "median"),
            mean_consensus_strength=("consensus_strength", "mean"),
            core_genes=("consensus_strength", lambda s: int((s >= 0.8).sum())),
        )
    )
    summary["pct_module_detected"] = 100 * summary.detected_genes / summary.module_genes
    total_detected = summary.groupby("chemical").detected_genes.transform("sum")
    summary["pct_chemical_detected_in_module"] = 100 * summary.detected_genes / total_detected
    detected_n = summary.detected_genes.clip(lower=1)
    proportion_up = summary.up_detected / detected_n
    z = 1.959963984540054
    denominator = 1 + z**2 / detected_n
    center = (proportion_up + z**2 / (2 * detected_n)) / denominator
    half_width = z * np.sqrt(
        proportion_up * (1 - proportion_up) / detected_n + z**2 / (4 * detected_n**2)
    ) / denominator
    summary["pct_up_among_detected"] = 100 * proportion_up.where(summary.detected_genes > 0)
    summary["pct_up_wilson_low"] = 100 * (center - half_width).where(summary.detected_genes > 0)
    summary["pct_up_wilson_high"] = 100 * (center + half_width).where(summary.detected_genes > 0)
    return gene_level, summary
