from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

from .science import CHEMICALS


def detection_matrix(
    frames: Mapping[str, pd.DataFrame], threshold: float = 0.05
) -> pd.DataFrame:
    """Return a gene by chemical boolean detection matrix at a declared FDR threshold."""
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    genes = frames[CHEMICALS[0]].index
    matrix = pd.DataFrame(
        {chemical: frames[chemical].reindex(genes).FDR.lt(threshold) for chemical in CHEMICALS},
        index=genes,
    )
    matrix.index.name = "WB_id"
    return matrix


def exact_intersection_table(
    detected: pd.DataFrame, include_undetected: bool = False
) -> pd.DataFrame:
    """Count exact multi-chemical detection patterns for an UpSet representation."""
    if list(detected.columns) != list(CHEMICALS):
        raise ValueError("detected columns must follow the declared chemical order")
    patterns = (
        detected.astype(bool)
        .groupby(list(CHEMICALS), observed=True)
        .size()
        .rename("gene_count")
        .reset_index()
    )
    patterns["degree"] = patterns[list(CHEMICALS)].sum(axis=1).astype(int)
    if not include_undetected:
        patterns = patterns[patterns.degree > 0]
    patterns["intersection"] = patterns.apply(
        lambda row: " & ".join(c for c in CHEMICALS if bool(row[c])) or "none", axis=1
    )
    return patterns.sort_values(["gene_count", "degree", "intersection"], ascending=[False, False, True]).reset_index(drop=True)


def unique_shared_gene_sets(
    detected: pd.DataFrame,
) -> tuple[dict[str, set[str]], dict[str, set[str]], pd.DataFrame]:
    """Partition each chemical's detected genes into exclusive and shared categories."""
    multiplicity = detected.sum(axis=1)
    unique_sets: dict[str, set[str]] = {}
    shared_sets: dict[str, set[str]] = {}
    rows = []
    for chemical in CHEMICALS:
        active = detected[chemical]
        unique = set(detected.index[active & multiplicity.eq(1)])
        shared = set(detected.index[active & multiplicity.ge(2)])
        unique_sets[chemical], shared_sets[chemical] = unique, shared
        rows.append(
            {
                "chemical": chemical,
                "detected_total": int(active.sum()),
                "detected_exclusively_here": len(unique),
                "detected_and_shared": len(shared),
            }
        )
    return unique_sets, shared_sets, pd.DataFrame(rows)


def pairwise_detection_overlap(
    frames: Mapping[str, pd.DataFrame], detected: pd.DataFrame
) -> pd.DataFrame:
    """Summarize overlap and sign concordance for all 45 chemical pairs."""
    rows = []
    for left_index, left in enumerate(CHEMICALS):
        for right in CHEMICALS[left_index + 1 :]:
            left_set = detected[left]
            right_set = detected[right]
            both = left_set & right_set
            neither = ~left_set & ~right_set
            only_left = left_set & ~right_set
            only_right = ~left_set & right_set
            odds_ratio, p_value = fisher_exact(
                [[int(both.sum()), int(only_left.sum())],
                 [int(only_right.sum()), int(neither.sum())]],
                alternative="greater",
            )
            intersection = int(both.sum())
            union = int((left_set | right_set).sum())
            minimum = min(int(left_set.sum()), int(right_set.sum()))
            overlap_ids = detected.index[both]
            signs_match = np.sign(frames[left].reindex(overlap_ids).logFC.to_numpy()) == np.sign(
                frames[right].reindex(overlap_ids).logFC.to_numpy()
            )
            rows.append(
                {
                    "chemical_a": left,
                    "chemical_b": right,
                    "overlap_genes": intersection,
                    "jaccard": intersection / union if union else np.nan,
                    "overlap_coefficient": intersection / minimum if minimum else np.nan,
                    "concordant_sign_pct": 100 * float(signs_match.mean()) if intersection else np.nan,
                    "fisher_odds_ratio": odds_ratio,
                    "p_raw": p_value,
                }
            )
    result = pd.DataFrame(rows)
    result["q_bh_45_pairs"] = multipletests(result.p_raw, method="fdr_bh")[1]
    return result.sort_values(["q_bh_45_pairs", "jaccard"], ascending=[True, False]).reset_index(drop=True)


def chemical_module_preference(
    detected: pd.DataFrame, assignments: pd.Series
) -> pd.DataFrame:
    """Test all 50 chemical-module detection enrichments in one declared family."""
    genes = assignments.index.intersection(detected.index)
    labels = assignments.reindex(genes)
    rows = []
    total_genes = len(genes)
    for chemical in CHEMICALS:
        active = detected.reindex(genes)[chemical]
        total_active = int(active.sum())
        for module in sorted(labels.unique()):
            in_module = labels.eq(module)
            observed = int((active & in_module).sum())
            module_size = int(in_module.sum())
            expected = total_active * module_size / total_genes
            only_detected = total_active - observed
            only_module = module_size - observed
            neither = total_genes - observed - only_detected - only_module
            odds_ratio, p_value = fisher_exact(
                [[observed, only_detected], [only_module, neither]],
                alternative="two-sided",
            )
            rows.append(
                {
                    "chemical": chemical,
                    "module": int(module),
                    "clustered_genes": total_genes,
                    "module_genes": module_size,
                    "detected_genes": total_active,
                    "observed": observed,
                    "expected": expected,
                    "log2_observed_expected": np.log2((observed + 0.5) / (expected + 0.5)),
                    "pearson_residual": (observed - expected) / np.sqrt(expected) if expected else np.nan,
                    "fisher_odds_ratio": odds_ratio,
                    "p_raw": p_value,
                }
            )
    result = pd.DataFrame(rows)
    result["q_bh_50_cells"] = multipletests(result.p_raw, method="fdr_bh")[1]
    return result


def plot_upset(
    intersections: pd.DataFrame, output, top_n: int = 30
) -> None:
    """Draw a deterministic UpSet plot without a hidden intersection filter."""
    shown = intersections.head(top_n).iloc[::-1].reset_index(drop=True)
    fig, (ax_bar, ax_matrix) = plt.subplots(
        1, 2, figsize=(13, max(7, 0.28 * len(shown))),
        gridspec_kw={"width_ratios": [2.2, 3.2]}, sharey=True,
    )
    y = np.arange(len(shown))
    ax_bar.barh(y, shown.gene_count, color="#4C72B0")
    ax_bar.set_xlabel("Genes in exact intersection")
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(shown.intersection, fontsize=7)
    ax_bar.grid(axis="x", alpha=0.2)
    for yi, value in zip(y, shown.gene_count, strict=True):
        ax_bar.text(value, yi, f" {value}", va="center", fontsize=7)
    for yi, row in shown.iterrows():
        active = [j for j, chemical in enumerate(CHEMICALS) if bool(row[chemical])]
        ax_matrix.scatter(range(len(CHEMICALS)), [yi] * len(CHEMICALS), s=12, color="#D0D0D0")
        if active:
            ax_matrix.plot(active, [yi] * len(active), color="#333333", linewidth=1)
            ax_matrix.scatter(active, [yi] * len(active), s=30, color="#C44E52", zorder=3)
    ax_matrix.set_xticks(range(len(CHEMICALS)))
    ax_matrix.set_xticklabels(CHEMICALS, rotation=90)
    ax_matrix.set_yticks([])
    ax_matrix.set_xlabel("Detected at FDR < 0.05")
    fig.suptitle(f"Top {min(top_n, len(shown))} exact PFAS DEG intersections")
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_detection_composition(summary: pd.DataFrame, output) -> None:
    """Plot threshold-defined exclusive and shared detections for each chemical."""
    ordered = summary.set_index("chemical").reindex(CHEMICALS)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(ordered.index, ordered.detected_exclusively_here, label="exclusive to one contrast")
    ax.bar(
        ordered.index,
        ordered.detected_and_shared,
        bottom=ordered.detected_exclusively_here,
        label="detected in at least two contrasts",
    )
    ax.set_ylabel("Genes detected at FDR < 0.05")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False)
    ax.set_title("Per-chemical DEG detection composition")
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def module_activity_pca(matrix: pd.DataFrame, assignments: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return consensus-module mean activity and a two-component condition PCA."""
    activity = matrix.assign(module=assignments.reindex(matrix.index)).groupby("module").mean()
    model = PCA(n_components=2)
    scores = model.fit_transform(activity.T)
    pca = pd.DataFrame(scores, index=activity.columns, columns=["PC1", "PC2"])
    pca["explained_PC1"] = model.explained_variance_ratio_[0]
    pca["explained_PC2"] = model.explained_variance_ratio_[1]
    return activity, pca


def plot_module_activity(activity: pd.DataFrame, pca: pd.DataFrame, output) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    image = axes[0].imshow(activity.T, aspect="auto", cmap="RdBu_r")
    axes[0].set_xticks(range(len(activity.index)), [f"M{x}" for x in activity.index])
    axes[0].set_yticks(range(len(activity.columns)), activity.columns)
    axes[0].set_title("Mean row-z-scored logFC by consensus module")
    fig.colorbar(image, ax=axes[0], fraction=0.046)
    axes[1].scatter(pca.PC1, pca.PC2, color="#4C72B0")
    for chemical, row in pca.iterrows():
        axes[1].annotate(chemical, (row.PC1, row.PC2), xytext=(4, 3), textcoords="offset points")
    axes[1].set_xlabel(f"PC1 ({100 * pca.explained_PC1.iloc[0]:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({100 * pca.explained_PC2.iloc[0]:.1f}%)")
    axes[1].set_title("Exploratory PCA of ten chemical conditions")
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def adjust_enrichment_family(results: pd.DataFrame) -> pd.DataFrame:
    """Apply one BH correction to a declared enrichment result family."""
    required = {"p_raw", "enrichment"}
    if not required.issubset(results):
        raise ValueError(f"missing enrichment columns: {sorted(required - set(results))}")
    adjusted = results.copy()
    adjusted["p_fdr_bh_family"] = multipletests(adjusted.p_raw, method="fdr_bh")[1]
    return adjusted


def ontology_edges_for_terms(go_ids: list[str], godag) -> tuple[list[str], list[tuple[str, str]]]:
    """Connect selected terms to their nearest selected ontology ancestors."""
    selected = [go_id for go_id in go_ids if go_id in godag]
    edges: list[tuple[str, str]] = []
    for child in selected:
        ancestors = set(godag[child].get_all_parents()).intersection(selected)
        if ancestors:
            parent = max(ancestors, key=lambda go_id: godag[go_id].depth)
            edges.append((parent, child))
    return selected, edges


def plot_ontology_panels(
    enrichment: pd.DataFrame,
    group_column: str,
    groups: list,
    godag,
    output,
    top_n: int = 8,
) -> None:
    """Plot ontology-backed enriched-term forests for modules or chemicals."""
    ncols = 1 if len(groups) <= 5 else 2
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, 6.5 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, group in zip(axes, groups, strict=False):
        subset = enrichment[enrichment[group_column] == group].nsmallest(top_n, "p_fdr_bh_family")
        nodes, edges = ontology_edges_for_terms(subset.GO_id.tolist(), godag)
        graph = nx.DiGraph(edges)
        graph.add_nodes_from(nodes)
        if not nodes:
            ax.text(0.5, 0.5, "No family-wide significant terms", ha="center", va="center")
            ax.set_axis_off()
            ax.set_title(str(group))
            continue
        positions = nx.kamada_kawai_layout(graph)
        labels = {
            go_id: "\n".join(textwrap.wrap(godag[go_id].name, width=28))
            for go_id in nodes
        }
        nx.draw_networkx_edges(
            graph, positions, ax=ax, arrows=True, alpha=0.45,
            edge_color="#52606D", arrowsize=14,
        )
        nx.draw_networkx_nodes(
            graph, positions, ax=ax, node_size=1500, node_color="#D9EAF2",
            edgecolors="#287094", linewidths=0.8,
        )
        nx.draw_networkx_labels(
            graph, positions, labels=labels, ax=ax, font_size=7,
            bbox={"facecolor": "white", "edgecolor": "#287094", "alpha": 0.92, "pad": 2},
        )
        ax.margins(0.30)
        ax.set_title(str(group))
        ax.set_axis_off()
    for ax in axes[len(groups):]:
        ax.set_axis_off()
    fig.suptitle("Ontology relationships among leading family-wide significant terms")
    fig.tight_layout(rect=(0, 0, 1, 0.985), pad=3.0)
    fig.savefig(output, dpi=200, bbox_inches="tight", pad_inches=0.4)
    fig.savefig(Path(output).with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.4)
    plt.close(fig)
