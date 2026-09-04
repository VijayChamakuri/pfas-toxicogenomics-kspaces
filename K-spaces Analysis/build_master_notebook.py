"""Build the unified, restart-safe PFAS scientific computing notebook."""

from __future__ import annotations

import copy
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import nbformat as nbf

ROOT = Path(__file__).resolve().parent
CANONICAL = ROOT / "PFAS_kspaces_analysis.ipynb"
OUTPUT = ROOT / "PFAS_master_analysis.ipynb"


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(text.strip())


def replace_em_dashes(value: Any) -> Any:
    """Replace project-authored em dashes throughout copied notebook content."""
    if isinstance(value, str):
        return value.replace("\N{EM DASH}", ":")
    if isinstance(value, list):
        return [replace_em_dashes(item) for item in value]
    if isinstance(value, dict):
        return {key: replace_em_dashes(item) for key, item in value.items()}
    return value


def canonical_cell(cell: nbf.NotebookNode) -> nbf.NotebookNode:
    """Copy a canonical cell and subordinate its headings to the master narrative."""
    copied = nbf.from_dict(replace_em_dashes(copy.deepcopy(dict(cell))))
    if copied.cell_type == "markdown":
        copied.source = re.sub(r"(?m)^####\s+", "##### ", copied.source)
        copied.source = re.sub(r"(?m)^###\s+", "#### ", copied.source)
        copied.source = re.sub(r"(?m)^##\s+", "### ", copied.source)
        copied.source = re.sub(r"(?m)^(#{3,6})\s+\d+\.\s+", r"\1 ", copied.source)
        copied.source = re.sub(
            r"(?m)^(#{3,6}\s+.*?)Section\s+\d+\s+", r"\1", copied.source
        )
    return copied


def append_range(
    destination: list[nbf.NotebookNode], source: nbf.NotebookNode, start: int, stop: int
) -> None:
    destination.extend(canonical_cell(cell) for cell in source.cells[start:stop])


generated_intermediate = not CANONICAL.exists()
if generated_intermediate:
    subprocess.run([sys.executable, str(ROOT / "build_notebook.py")], cwd=ROOT, check=True)

canonical = nbf.read(CANONICAL, as_version=4)
cells: list[nbf.NotebookNode] = [
    md("""
# Integrated PFAS toxicogenomics analysis

This notebook connects data provenance, reconstruction of the initial analytical approach,
methodological correction, uncertainty-aware module discovery, functional interpretation,
backend-independent workflow classification, and evidence-grounded analysis planning. The ten
canonical edgeR result tables are the only inputs to biological inference. Synthetic workflow
requests verify software behavior and never enter PFAS conclusions.

Decisions remain traceable to source files, executed cells, versioned artifacts, and explicit
status fields. Human scientific review is required before an accepted plan is executed or an
exploratory association is reported as a conclusion.
"""),
    md("""
## 1. Research Context and Analytical Objectives

The study asks which transcriptional responses are shared across ten PFAS exposures and which
are condition-associated in *C. elegans*. The analytical contract uses the measured 13,852-gene
universe, preserves continuous effect sizes, separates detection thresholds from magnitudes,
and quantifies optimizer and resampling uncertainty separately. The multi-fit consensus is the
recommended module assignment. Raw five-module percentage correlations are not valid evidence
of chemical similarity. Corrected residual pair tests control BH FDR across 45 pairs separately
at each tested k, but remain exploratory because clustering uncertainty and external validation
are unresolved.
"""),
    md("""
## 2. Reproducible Environment and Data Provenance

The setup resolves project paths, records the runtime, and verifies all ten canonical DGE tables
before analysis begins. This is the single execution notebook. Maintained helper modules, input data,
and generated outputs stay mapped to their project locations so the notebook does not duplicate or
silently diverge from tested source code.

| Role | Mapped location relative to `K-spaces Analysis` |
|---|---|
| Ten canonical edgeR inputs | `../Data/DEGs/{CHEM}vsControl_DGE_results.csv` |
| Scientific helper code | `src/pfas_master/` |
| Workflow and RAG code | `../integrated_system/src/pfas_workflow/` and `../integrated_system/llm_system/` |
| Tables and gene lists | `output/` |
| Figures, trees, and UpSet plots | `output/figures/` |

Launch Jupyter from the repository root or the `K-spaces Analysis` directory. If a different
working directory is required, set `PFAS_KSPACES_DIR` to the repository's analysis directory.
"""),
    code("""
from pathlib import Path
import hashlib, json, os, platform, subprocess, sys
import pandas as pd

configured_root = os.environ.get('PFAS_KSPACES_DIR')
candidates = [
    Path(configured_root).expanduser() if configured_root else None,
    Path.cwd().resolve(),
    Path.cwd().resolve() / 'K-spaces Analysis',
    Path.cwd().resolve().parent / 'K-spaces Analysis',
]
ROOT = next((path.resolve() for path in candidates if path and (path / 'build_notebook.py').exists()), None)
if ROOT is None:
    raise RuntimeError(
        'Set PFAS_KSPACES_DIR to the K-spaces Analysis folder that contains build_notebook.py'
    )
PROJECT = ROOT.parent
SYSTEM_ROOT = PROJECT / 'integrated_system'
for required in (ROOT / 'src', SYSTEM_ROOT, SYSTEM_ROOT / 'src'):
    if str(required) not in sys.path:
        sys.path.insert(0, str(required))

environment = {
    'python': platform.python_version(),
    'analysis_root': str(ROOT),
    'system_root': str(SYSTEM_ROOT),
    'canonical_builder_sha256': hashlib.sha256((ROOT / 'build_notebook.py').read_bytes()).hexdigest(),
}
display(pd.Series(environment, name='value').to_frame())
"""),
]

# Preserve maintained setup and provenance code in its original execution order.
append_range(cells, canonical, 1, 7)
cells.append(md("""
## 3. Data Integrity and Quality Assessment

Quality control distinguishes the complete measured universe from threshold-defined detection
sets. Continuous log-fold changes feed module discovery. The initial zero-masked matrix is
reconstructed only to expose ties introduced when non-significance was treated as zero effect.
"""))
append_range(cells, canonical, 7, 9)
cells.extend([
    code("""
from pfas_master.science import (
    build_detection_evidence,
    build_legacy_zero_masked_matrix,
    load_and_validate_dge,
)
dge_frames, provenance = load_and_validate_dge(PROJECT / 'Data' / 'DEGs')
assert provenance['rows'].eq(13852).all()
evidence, threshold_sensitivity = build_detection_evidence(dge_frames)
legacy_matrix, legacy_diagnostics = build_legacy_zero_masked_matrix(dge_frames)
display(provenance)
display(threshold_sensitivity)
display(legacy_diagnostics)
"""),
    md("""
### Exact DEG intersections and per-chemical detection composition

The original UpSet and unique-versus-shared analyses asked useful questions, but their categories
depend on the FDR threshold. The corrected version preserves all continuous effects in the primary
matrix and uses FDR below 0.05 only to describe exact detection intersections. “Exclusive” means
detected in exactly one of these ten contrasts; it does not mean biologically absent elsewhere.
"""),
    code("""
from pfas_master.deep_analysis import (
    detection_matrix,
    exact_intersection_table,
    plot_detection_composition,
    plot_upset,
    pairwise_detection_overlap,
    unique_shared_gene_sets,
)

deep_output = ROOT / 'output' / 'deep_analysis'
deep_figures = FIG_DIR / 'deep_analysis'
deep_output.mkdir(parents=True, exist_ok=True)
deep_figures.mkdir(parents=True, exist_ok=True)

detected_matrix = detection_matrix(dge_frames, threshold=FDR_THRESHOLD)
exact_intersections = exact_intersection_table(detected_matrix)
unique_gene_sets, shared_gene_sets, detection_composition = unique_shared_gene_sets(detected_matrix)
pairwise_overlap = pairwise_detection_overlap(dge_frames, detected_matrix)

up_detected_matrix = pd.DataFrame({
    chemical: detected_matrix[chemical] & dge_frames[chemical].logFC.gt(0)
    for chemical in CHEMICALS
})
down_detected_matrix = pd.DataFrame({
    chemical: detected_matrix[chemical] & dge_frames[chemical].logFC.lt(0)
    for chemical in CHEMICALS
})
up_intersections = exact_intersection_table(up_detected_matrix)
down_intersections = exact_intersection_table(down_detected_matrix)

exact_intersections.to_csv(deep_output / 'exact_DEG_intersections.csv', index=False)
up_intersections.to_csv(deep_output / 'exact_upregulated_DEG_intersections.csv', index=False)
down_intersections.to_csv(deep_output / 'exact_downregulated_DEG_intersections.csv', index=False)
detection_composition.to_csv(deep_output / 'per_chemical_unique_shared_counts.csv', index=False)
pairwise_overlap.to_csv(deep_output / 'pairwise_DEG_overlap_and_sign_concordance.csv', index=False)
detected_matrix.astype(int).to_csv(deep_output / 'gene_by_chemical_detection_matrix.csv')

exclusive_rows = []
for owner, genes in unique_gene_sets.items():
    for gene in sorted(genes):
        row = {'WB_id': gene, 'exclusive_detection_contrast': owner}
        for chemical in CHEMICALS:
            row[f'logFC_{chemical}'] = dge_frames[chemical].loc[gene, 'logFC']
            row[f'FDR_{chemical}'] = dge_frames[chemical].loc[gene, 'FDR']
        exclusive_rows.append(row)
exclusive_continuous_evidence = pd.DataFrame(exclusive_rows)
exclusive_continuous_evidence.to_csv(
    deep_output / 'exclusive_detection_genes_with_all_continuous_effects.csv', index=False
)
plot_upset(exact_intersections, deep_figures / 'upset_exact_DEG_intersections.png', top_n=30)
plot_upset(up_intersections, deep_figures / 'upset_exact_upregulated_intersections.png', top_n=30)
plot_upset(down_intersections, deep_figures / 'upset_exact_downregulated_intersections.png', top_n=30)
plot_detection_composition(detection_composition, deep_figures / 'per_chemical_DEG_composition.png')

display(detection_composition)
display(pairwise_overlap)
display(exact_intersections.head(30))
display(Image(filename=str(deep_figures / 'upset_exact_DEG_intersections.png')))
display(Image(filename=str(deep_figures / 'upset_exact_upregulated_intersections.png')))
display(Image(filename=str(deep_figures / 'upset_exact_downregulated_intersections.png')))
display(Image(filename=str(deep_figures / 'per_chemical_DEG_composition.png')))
"""),
    md("""
## 4. Reconstruction of the Original Analytical Workflow

The original model-selection and stability steps are retained for auditability. Their outputs
are diagnostic records, not final scientific findings. The initial resampling code discarded
subset labels and measured repeated full-data optimizer variation instead of gene-subsample
stability. The revised strategy keeps shared-gene labels and reports optimizer and sampling
variation separately.
"""),
])
append_range(cells, canonical, 9, 13)
cells.append(md("""
## 5. Methodological Corrections and Validation Strategy

The revised workflow uses continuous effects, independent repeated fits, true gene subsamples,
null simulations, consensus voting, chemical-level jackknife checks, and permutation tests.
Threshold categories describe detection evidence and are not chemical-specific mechanisms.
"""))
append_range(cells, canonical, 13, 18)
cells.append(md("""
### Audit-only single-fit reconstruction

The cell below recreates the primary single fit only because later corrected analyses require its
assignment for direct diagnostics. Raw module-percentage similarity tables, rankings, heatmaps,
and structure-activity figures are intentionally excluded. They are retained only in the internal
archive and must not be cited as scientific results.
"""))
append_range(cells, canonical, 19, 20)
cells.append(md("""
## 6. PFAS Response Patterns and Statistical Evidence

Module profiles, sensitivity variants, consensus assignments, null-adjusted comparisons, and
structure-activity probes are interpreted together. The ten chemicals are conditions rather
than biological replicates, so panel-level patterns do not establish generalization or causality.
"""))
append_range(cells, canonical, 47, 56)
cells.append(code("""
def module_profile_from_assignment(assign: pd.Series, matrix_df: pd.DataFrame, k: int):
    mat_df = matrix_df.loc[assign.index].copy()
    mat_df['module'] = assign
    mean_value = mat_df.groupby('module')[list(matrix_df.columns)].mean().round(3)
    pct_rows = []
    for chemical in CHEMICALS:
        detected = merged.index[merged[f'FDR_{chemical}'] < FDR_THRESHOLD]
        detected_assignments = assign.loc[assign.index.intersection(detected)]
        counts = detected_assignments.value_counts().reindex(range(k), fill_value=0)
        row = (counts / len(detected) * 100).round(1)
        row.name = chemical
        pct_rows.append(row)
    pct_table = pd.DataFrame(pct_rows)
    pct_table.index.name = 'chemical'
    return mean_value, pct_table

meanval_consensus, pct_consensus = module_profile_from_assignment(
    assign_consensus, matrix_primary, K
)
core_mask = consensus_strength_s >= 0.8
assign_core = assign_consensus[core_mask]
print(f'Consensus profile prepared for {len(assign_consensus)} genes; '
      f'{core_mask.sum()} have at least 80% assignment agreement.')
"""))
cells.extend([
    md("""
### Per-chemical and per-module gene contributions

Every chemical-by-module combination is regenerated from the consensus assignment. The detailed
gene table retains logFC, FDR, direction, detection status, module, and consensus strength. The
summary reports counts and proportions without treating thresholded non-detection as zero effect.
Separate CSV files make all 50 combinations inspectable outside the notebook.
"""),
    code("""
from pfas_master.science import build_module_contributions
from pfas_master.deep_analysis import (
    chemical_module_preference,
    module_activity_pca,
    plot_module_activity,
)

gene_contributions, contribution_summary = build_module_contributions(
    dge_frames, assign_consensus, consensus_strength_s
)
gene_contributions.to_csv(deep_output / 'gene_contributions_all_chemicals_modules.csv', index=False)
contribution_summary.to_csv(deep_output / 'summary_gene_contributions.csv', index=False)
module_preference = chemical_module_preference(detected_matrix, assign_consensus)
module_preference.to_csv(deep_output / 'chemical_module_preference_50_tests.csv', index=False)

gene_list_dir = deep_output / 'gene_lists_by_chemical_module'
gene_list_dir.mkdir(parents=True, exist_ok=True)
for (chemical, module), group in gene_contributions.groupby(['chemical', 'consensus_module']):
    group.sort_values(['is_detected', 'FDR', 'logFC'], ascending=[False, True, False]).to_csv(
        gene_list_dir / f'{chemical}_module_{module}_genes.csv', index=False
    )

fig, axes = plt.subplots(1, 2, figsize=(15, 5))
summary_pivot = contribution_summary.pivot(
    index='chemical', columns='consensus_module', values='pct_module_detected'
).reindex(CHEMICALS)
image = axes[0].imshow(summary_pivot, aspect='auto', cmap='YlOrRd')
axes[0].set_xticks(range(summary_pivot.shape[1]), [f'M{x}' for x in summary_pivot.columns])
axes[0].set_yticks(range(len(summary_pivot)), summary_pivot.index)
axes[0].set_title('Percent of each module detected per chemical')
fig.colorbar(image, ax=axes[0], label='Percent of module genes')

up = contribution_summary.pivot(index='chemical', columns='consensus_module', values='up_detected').reindex(CHEMICALS)
down = contribution_summary.pivot(index='chemical', columns='consensus_module', values='down_detected').reindex(CHEMICALS)
x = np.arange(len(CHEMICALS)); width = 0.075
for j, module in enumerate(up.columns):
    offset = (j - (len(up.columns) - 1) / 2) * width
    axes[1].bar(x + offset, up[module], width, label=f'M{module} up')
    axes[1].bar(x + offset, -down[module], width, alpha=0.45)
axes[1].axhline(0, color='black', linewidth=0.7)
axes[1].set_xticks(x, CHEMICALS, rotation=45)
axes[1].set_ylabel('Detected genes: up positive, down negative')
axes[1].set_title('Direction-resolved contribution by chemical and module')
axes[1].legend(ncol=3, fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(deep_figures / 'per_chemical_per_module_contributions.png', dpi=200, bbox_inches='tight')
plt.close(fig)

preference_heat = module_preference.pivot(
    index='chemical', columns='module', values='log2_observed_expected'
).reindex(CHEMICALS)
fig, ax = plt.subplots(figsize=(8, 6))
limit = np.nanmax(np.abs(preference_heat.to_numpy()))
image = ax.imshow(preference_heat, aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit)
ax.set_xticks(range(preference_heat.shape[1]), [f'M{x}' for x in preference_heat.columns])
ax.set_yticks(range(len(preference_heat)), preference_heat.index)
ax.set_title('Module-size-adjusted detection preference')
for row in range(preference_heat.shape[0]):
    for column in range(preference_heat.shape[1]):
        record = module_preference[
            module_preference.chemical.eq(preference_heat.index[row])
            & module_preference.module.eq(preference_heat.columns[column])
        ].iloc[0]
        marker = '*' if record.q_bh_50_cells < 0.05 else ''
        ax.text(column, row, f'{preference_heat.iloc[row, column]:.2f}{marker}', ha='center', va='center', fontsize=8)
fig.colorbar(image, ax=ax, label='log2((observed + 0.5)/(expected + 0.5))')
fig.tight_layout()
fig.savefig(deep_figures / 'chemical_module_preference_heatmap.png', dpi=200, bbox_inches='tight')
plt.close(fig)

module_activity, activity_pca = module_activity_pca(matrix_primary, assign_consensus)
module_activity.to_csv(deep_output / 'consensus_module_activity.csv')
activity_pca.to_csv(deep_output / 'chemical_module_activity_PCA.csv')
plot_module_activity(module_activity, activity_pca, deep_figures / 'module_activity_and_PCA.png')

display(contribution_summary)
display(module_preference)
display(Image(filename=str(deep_figures / 'per_chemical_per_module_contributions.png')))
display(Image(filename=str(deep_figures / 'chemical_module_preference_heatmap.png')))
display(Image(filename=str(deep_figures / 'module_activity_and_PCA.png')))
"""),
])
cells.append(md("""
### Quarantined raw similarity diagnostic

The raw percent-per-module correlation was stress-tested and rejected because unequal module
sizes create a high null correlation. Its pair rankings, figures, and candidate tables are not
displayed in this notebook. The corrected analysis below subtracts the module-size expectation,
tests against label-shuffle nulls, and treats surviving pairs as exploratory.
"""))
append_range(cells, canonical, 66, 75)

# Consensus sensitivity and corrected residual analyses are maintained in dependency order.
append_range(cells, canonical, 82, 84)
cells.extend([
    code("""
cross_strategy_consensus = []
for name, other_assign in [
    ('percentile_rank_sensitivity', assign_consensus_sens),
    ('topN1000', assign_consensus_top1000),
]:
    shared = assign_consensus.index.intersection(other_assign.index)
    cross_strategy_consensus.append({
        'variant': name,
        'n_shared_genes': len(shared),
        'ARI_consensus_vs_consensus': round(
            adjusted_rand_score(assign_consensus.loc[shared], other_assign.loc[shared]), 3
        ),
    })
cross_strategy_consensus_df = pd.DataFrame(cross_strategy_consensus)
cross_strategy_consensus_df.to_csv(
    ROOT / 'output' / 'kspaces_runs' / 'cross_strategy_ARI_consensus.csv', index=False
)
display(cross_strategy_consensus_df)
"""),
    md("""
Consensus-to-consensus ARI isolates input-strategy disagreement from single-fit optimizer noise.
The values are sensitivity evidence for the partition, not chemical-similarity rankings.
"""),
])
cells.append(code("""
from itertools import combinations
from scipy.stats import spearmanr

ATTRS = pd.DataFrame([
    ('PFOSA', 'Long', 'Sulfonamide', 'No'), ('PFBSA', 'Short', 'Sulfonamide', 'No'),
    ('PFNA', 'Long', 'Carboxylic acid', 'No'), ('PFOS', 'Long', 'Sulfonic acid', 'No'),
    ('PFOA', 'Long', 'Carboxylic acid', 'No'), ('GenX', 'Short', 'Carboxylic acid', 'Yes'),
    ('PFEESA', 'Short', 'Sulfonic acid', 'Yes'), ('PFPeA', 'Short', 'Carboxylic acid', 'No'),
    ('PFBA', 'Short', 'Carboxylic acid', 'No'), ('PFBS', 'Short', 'Sulfonic acid', 'No'),
], columns=['chemical', 'chain_length', 'functional_group', 'ether']).set_index('chemical')
pairs_all = list(combinations(CHEMICALS, 2))
candidate_pairs = [('GenX', 'PFBS'), ('PFOS', 'PFNA'), ('PFOS', 'PFOA'), ('PFOSA', 'PFPeA')]
"""))
append_range(cells, canonical, 86, 108)
cells.append(md("""
## 7. Functional and Pathway-Level Interpretation

Functional summaries connect to consensus assignments and explicitly declared eligible populations.
Identifier mapping loss and ontology version are result provenance. Consensus-module enrichment uses
the clustered DEG-union genes mapped in each namespace. Chemical-exclusive and chemical-shared studies
use the measured, namespace-mapped universe, while each chemical-by-module study uses the mapped genes
inside that module. Each declared family receives one BH correction across all of its tested terms.
Terms remain exploratory associations because annotation and module uncertainty are not fully propagated.
"""))
append_range(cells, canonical, 36, 41)
cells.append(code("""
col_order = [
    'PFOS', 'PFNA', 'PFOA', 'GenX', 'PFBS',
    'PFEESA', 'PFBA', 'PFBSA', 'PFOSA', 'PFPeA',
]
import plotly.graph_objects as go

module_colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#64B5CD']

def hex_to_rgba(hex_color, alpha):
    value = hex_color.lstrip('#')
    red, green, blue = int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)
    return f'rgba({red},{green},{blue},{alpha})'
"""))
append_range(cells, canonical, 75, 82)
cells.extend([
    md("""
### Corrected enrichment across chemicals, modules, and detection categories

Three prespecified enrichment families restore the original analytical depth: 50 chemical-by-module
detected-gene sets, ten chemical-exclusive sets, and ten chemical-shared sets. Each family uses
the namespace-mapped measured population and one BH correction across every module, chemical,
namespace, and term tested within that family. Empty or underpowered sets are recorded rather than
silently dropped. The resulting ontology panels use actual GO parent relationships; no manual
GO-slim category assignment or artificial Sankey padding is used.
"""),
    code("""
from pfas_master.deep_analysis import adjust_enrichment_family, plot_ontology_panels

MIN_GO_STUDY_GENES = 10

def run_enrichment_family(studies, family_name, population_mode, top_n=15):
    rows = []
    skipped = []
    with contextlib.redirect_stdout(io.StringIO()):
        for namespace, associations in ns2assoc.items():
            engine_cache = {}
            for study_id, study_genes in studies.items():
                if population_mode == 'measured':
                    population_key = 'measured'
                    population = measured_background_by_ns[namespace]
                elif population_mode == 'within_module':
                    module = int(str(study_id).split('|M')[1])
                    population_key = f'module_{module}'
                    population = sorted(
                        set(assign_consensus_aligned.index[assign_consensus_aligned == module])
                        .intersection(associations)
                    )
                else:
                    raise ValueError(f'unknown population mode: {population_mode}')
                if population_key not in engine_cache:
                    engine_cache[population_key] = GOEnrichmentStudy(
                        population, associations, godag, propagate_counts=True,
                        alpha=FDR_ALPHA, methods=['fdr_bh'],
                    )
                engine = engine_cache[population_key]
                eligible = sorted(set(study_genes).intersection(population))
                if len(eligible) < MIN_GO_STUDY_GENES:
                    skipped.append({
                        'family': family_name, 'study_id': str(study_id),
                        'namespace': namespace, 'input_genes': len(set(study_genes)),
                        'mapped_eligible_genes': len(eligible),
                        'reason': f'fewer than {MIN_GO_STUDY_GENES} mapped eligible genes',
                    })
                    continue
                for result in engine.run_study(eligible):
                    rows.append({
                        'family': family_name, 'study_id': str(study_id),
                        'namespace': namespace, 'GO_id': result.GO, 'term': result.name,
                        'enrichment': result.enrichment, 'p_raw': result.p_uncorrected,
                        'study_count': result.study_count, 'study_n': result.study_n,
                        'pop_count': result.pop_count, 'pop_n': result.pop_n,
                    })
    tested = adjust_enrichment_family(pd.DataFrame(rows))
    significant = tested[
        tested.enrichment.eq('e') & tested.p_fdr_bh_family.lt(FDR_ALPHA)
    ].sort_values(['study_id', 'namespace', 'p_fdr_bh_family'])
    displayed = significant.groupby(['study_id', 'namespace'], group_keys=False).head(top_n)
    tested.to_csv(deep_output / f'GO_{family_name}_all_tests.csv', index=False)
    displayed.to_csv(deep_output / f'GO_{family_name}_significant_top{top_n}.csv', index=False)
    skipped_df = pd.DataFrame(skipped, columns=[
        'family', 'study_id', 'namespace', 'input_genes',
        'mapped_eligible_genes', 'reason',
    ])
    skipped_df.to_csv(deep_output / f'GO_{family_name}_skipped.csv', index=False)
    return tested, displayed, skipped_df

chemical_module_studies = {}
for chemical in CHEMICALS:
    detected_genes = set(detected_matrix.index[detected_matrix[chemical]])
    for module in sorted(assign_consensus_aligned.unique()):
        module_genes = set(assign_consensus_aligned.index[assign_consensus_aligned == module])
        chemical_module_studies[f'{chemical}|M{module}'] = detected_genes.intersection(module_genes)

cm_go_tests, cm_go, cm_go_skipped = run_enrichment_family(
    chemical_module_studies, 'chemical_by_module', 'within_module'
)
unique_go_tests, unique_go, unique_go_skipped = run_enrichment_family(
    unique_gene_sets, 'chemical_exclusive', 'measured'
)
shared_go_tests, shared_go, shared_go_skipped = run_enrichment_family(
    shared_gene_sets, 'chemical_shared', 'measured'
)

cm_go_dir = deep_output / 'GO_per_chemical_per_module'
cm_go_dir.mkdir(parents=True, exist_ok=True)
for study_id, group in cm_go.groupby('study_id'):
    group.to_csv(cm_go_dir / f'{study_id.replace("|", "_")}.csv', index=False)

cm_counts = (
    cm_go.groupby('study_id').size().reindex(chemical_module_studies, fill_value=0)
    .rename('significant_terms').reset_index()
)
cm_counts[['chemical', 'module']] = cm_counts.study_id.str.split('|', expand=True)
cm_heat = cm_counts.pivot(index='chemical', columns='module', values='significant_terms').reindex(CHEMICALS)
fig, ax = plt.subplots(figsize=(9, 6))
image = ax.imshow(cm_heat, aspect='auto', cmap='viridis')
ax.set_xticks(range(cm_heat.shape[1]), cm_heat.columns)
ax.set_yticks(range(len(cm_heat)), cm_heat.index)
ax.set_title('Family-wide significant GO terms per chemical and consensus module')
for row in range(cm_heat.shape[0]):
    for column in range(cm_heat.shape[1]):
        ax.text(column, row, int(cm_heat.iloc[row, column]), ha='center', va='center', color='white')
fig.colorbar(image, ax=ax, label='Significant enriched terms')
fig.tight_layout()
fig.savefig(deep_figures / 'GO_per_chemical_per_module_heatmap.png', dpi=200, bbox_inches='tight')
plt.close(fig)

plot_ontology_panels(
    go_df_consensus, 'module', sorted(go_df_consensus.module.unique()), godag,
    deep_figures / 'ontology_trees_consensus_modules.png', top_n=8,
)
plot_ontology_panels(
    unique_go.assign(chemical=unique_go.study_id), 'chemical', list(CHEMICALS), godag,
    deep_figures / 'ontology_trees_exclusive_by_chemical.png', top_n=6,
)
plot_ontology_panels(
    shared_go.assign(chemical=shared_go.study_id), 'chemical', list(CHEMICALS), godag,
    deep_figures / 'ontology_trees_shared_by_chemical.png', top_n=6,
)

def gene_set_sankey(enrichment, title, stem):
    best = enrichment.sort_values('p_fdr_bh_family').groupby('study_id', as_index=False).first()
    chemicals = list(CHEMICALS)
    terms = list(dict.fromkeys(best.term.tolist()))
    nodes = chemicals + terms
    node_index = {name: index for index, name in enumerate(nodes)}
    figure = go.Figure(go.Sankey(
        node={'label': nodes, 'pad': 15, 'thickness': 14},
        link={
            'source': [node_index[row.study_id] for row in best.itertuples()],
            'target': [node_index[row.term] for row in best.itertuples()],
            'value': [max(int(row.study_count), 1) for row in best.itertuples()],
        },
    ))
    figure.update_layout(title_text=title, width=1200, height=650)
    figure.write_html(deep_figures / f'{stem}.html')
    figure.write_image(deep_figures / f'{stem}.png', scale=2)
    return figure

exclusive_sankey = gene_set_sankey(
    unique_go, 'PFAS-exclusive detected genes to leading corrected GO term',
    'sankey_exclusive_DEGs_to_GO',
)
shared_sankey = gene_set_sankey(
    shared_go, 'PFAS-shared detected genes to leading corrected GO term',
    'sankey_shared_DEGs_to_GO',
)

display(pd.DataFrame({
    'family': ['chemical_by_module', 'chemical_exclusive', 'chemical_shared'],
    'tested_terms': [len(cm_go_tests), len(unique_go_tests), len(shared_go_tests)],
    'displayed_significant_terms': [len(cm_go), len(unique_go), len(shared_go)],
    'skipped_study_namespaces': [len(cm_go_skipped), len(unique_go_skipped), len(shared_go_skipped)],
}))
display(cm_counts)
display(Image(filename=str(deep_figures / 'GO_per_chemical_per_module_heatmap.png')))
display(Image(filename=str(deep_figures / 'ontology_trees_consensus_modules.png')))
display(Image(filename=str(deep_figures / 'ontology_trees_exclusive_by_chemical.png')))
display(Image(filename=str(deep_figures / 'ontology_trees_shared_by_chemical.png')))
display(Image(filename=str(deep_figures / 'sankey_exclusive_DEGs_to_GO.png')))
display(Image(filename=str(deep_figures / 'sankey_shared_DEGs_to_GO.png')))
"""),
    md("""
### Reactome provenance boundary

Historical Reactome images exist, but the repository does not contain a versioned Reactome gene-set
collection, complete producing tables, or sufficient provenance to recompute them consistently.
They remain in the legacy audit and are not presented as corrected results. Adding Reactome requires
a pinned, legally usable gene-set release and a declared multiple-testing family; it is not replaced
with live web enrichment during a supposedly reproducible notebook run.
"""),
])
cells.extend([
    md("""
### Reassessment of the historical glutathione-transferase label

The earlier single-fit analysis highlighted a top GO label containing “glutathione transferase.”
Because the corrected continuous-effect primary matrix can change both module membership and
family-wide significant terms, the notebook tests whether that label is recovered instead of
assuming a fixed module number. Absence is a non-replication result, not an execution error.
"""),
    code("""
term_column = 'top_term_single_fit (Section 8)'
target_hits = compare_df.index[
    compare_df[term_column].str.contains('glutathione transferase', case=False, na=False)
]
historical_term_check = {
    'query': 'glutathione transferase',
    'recovered_among_corrected_top_terms': bool(len(target_hits)),
    'matching_modules': [int(value) for value in target_hits],
    'interpretation': (
        'Proceed to membership-overlap review before comparing labels.'
        if len(target_hits)
        else 'The historical top-term claim does not reproduce under the corrected primary matrix and GO family.'
    ),
}
display(pd.Series(historical_term_check, name='historical_term_reassessment').to_frame())
"""),
])
cells.append(md("""
### Scientific synthesis and interpretation boundaries

Consensus assignments replace single-fit labels for interpretation. Corrected residual analyses
show condition-associated structure but do not support the original raw-correlation rankings.
Functional terms and cross-condition pairs remain exploratory where clustering or annotation
uncertainty has not been fully propagated. Raw counts and complete sample metadata
are unavailable, so reproducibility begins with the canonical edgeR contrast tables.
"""))
cells.extend([
    md("""
## 8. Predictive Modeling and Stability Assessment

A shared modeling contract verifies workflow-request classification through independent neural
network implementations. Both backends use identical features, group-disjoint partitions,
preprocessing, capacity, seeds, early stopping, checkpoints, metrics, and prediction records.
This synthetic dataset tests implementation parity; it is isolated from PFAS data and supports
no biological claim.
"""),
    code("""
configured_system_python = os.environ.get('PFAS_SYSTEM_PYTHON')
venv_system_python = SYSTEM_ROOT / '.venv' / 'bin' / 'python'
system_python = Path(configured_system_python) if configured_system_python else venv_system_python
if not system_python.exists():
    system_python = Path(sys.executable)
parity_dir = SYSTEM_ROOT / 'artifacts' / 'framework_parity_notebook'
subprocess.run(
    [str(system_python), '-m', 'ml_frameworks', '--smoke', '--framework', 'both',
     '--output', str(parity_dir)],
    cwd=SYSTEM_ROOT,
    env={**os.environ, 'PYTHONPATH': str(SYSTEM_ROOT)},
    check=True,
)
parity_path = parity_dir / 'parity_report.json'
parity = json.loads(parity_path.read_text(encoding='utf-8'))
backend_rows = []
for backend, result in parity['results'].items():
    backend_rows.append({
        'backend': backend,
        'version': result['framework_version'],
        'seed': result['seed'],
        'parameters': result['parameter_count'],
        'epochs': result['epochs_completed'],
        'held_out_accuracy': result['metrics']['accuracy'],
        'held_out_macro_f1': result['metrics']['macro_f1'],
        'error_count': len(result['errors']),
        'checkpoint': result['checkpoint'],
    })
display(pd.Series(parity['dataset'], name='shared_contract').to_frame())
display(pd.DataFrame(backend_rows))
print('Differences on 16 synthetic held-out records test implementation consistency, not backend superiority.')
"""),
    md("""
## 9. Evidence-Grounded Analysis Planning

The planning entry point interprets a request, checks compatibility, retrieves reviewed evidence,
emits a structured workflow, validates parameters and scientific rules, detects unsupported
claims, and requires human approval. Missing metadata, nonexistent entities, incompatible
operations, insufficient evidence, and unsafe conclusions trigger clarification or abstention.

Model adaptation is restricted to request routing and structured-plan support. Its synthetic,
group-partitioned records do not alter retrieval, scientific validation, or human review.
"""),
    code("""
from pfas_workflow.catalog import ArtifactCatalog
from pfas_workflow.models import BiologicalRequest
from pfas_workflow.planner import RuleBasedPlanner
from llm_system.corpus import build_corpus

catalog = ArtifactCatalog(PROJECT)
planner = RuleBasedPlanner(catalog)
supported_plan = planner.plan(BiologicalRequest(
    question='Compare GenX and PFEESA using the corrected exploratory pair analysis'
))
unsupported_plan = planner.plan(BiologicalRequest(
    question='Prove that PFOS causes human liver cancer from the worm DGE tables'
))
display(pd.Series(supported_plan.model_dump(mode='json'), name='review_ready_plan').to_frame())
display(pd.Series(unsupported_plan.model_dump(mode='json'), name='abstained_plan').to_frame())

corpus_chunks = build_corpus()
evidence_inventory = pd.DataFrame([
    {
        'source_id': chunk.source_id,
        'source_path': chunk.relative_path,
        'source_sha256': chunk.source_sha256,
        'chunk_id': chunk.chunk_id,
        'heading': chunk.heading,
    }
    for chunk in corpus_chunks
])
display(
    evidence_inventory.groupby(['source_id', 'source_path', 'source_sha256'])
    .size().rename('reviewed_chunks').reset_index()
)
display(pd.DataFrame([citation.model_dump() for citation in supported_plan.citations]))

# This is a deterministic validation fixture, not a real scientific approval.
# Live execution stops here until an external reviewer records a signed decision.
validation_fixture = {
    'fixture_type': 'deterministic plan-validation example',
    'plan_digest': supported_plan.review_digest(),
    'schema_valid': True,
    'execution_status': 'not executed: external scientific approval required',
    'approval_claimed': False,
}
display(pd.Series(validation_fixture, name='validation_fixture').to_frame())

adaptation_path = SYSTEM_ROOT / 'llm_system' / 'training_output' / 'metrics.json'
run_adaptation = os.environ.get('PFAS_RUN_MODEL_ADAPTATION', '0') == '1'
if run_adaptation:
    subprocess.run(
        [str(system_python), '-m', 'llm_system.finetune', '--max-steps', '20'],
        cwd=SYSTEM_ROOT,
        env={**os.environ, 'PYTHONPATH': str(SYSTEM_ROOT)},
        check=True,
    )
if adaptation_path.exists():
    adaptation_metrics = json.loads(adaptation_path.read_text())
    adaptation_metrics['execution_mode'] = (
        'regenerated in this run' if run_adaptation else
        'previous bounded artifact; set PFAS_RUN_MODEL_ADAPTATION=1 to regenerate'
    )
    display(pd.Series(adaptation_metrics, name='adaptation_run').to_frame())
else:
    print('No adaptation artifact is present. Set PFAS_RUN_MODEL_ADAPTATION=1 to run the bounded path.')
"""),
    md("""
## 10. Workflow Reliability and Scientific Quality Control

Reliability cases cover supported, ambiguous, incomplete, incompatible, statistically invalid,
and hallucinated requests. Deterministic checks assess routing, workflow validity, grounding,
citations, unsupported claims, and abstention. Automated checks do not replace scientific review.
"""),
    code("""
from llm_system.evaluate import evaluate

quality_report = evaluate(SYSTEM_ROOT / 'llm_system' / 'eval_cases.json')
display(pd.DataFrame({
    system: result['metrics'] for system, result in quality_report['systems'].items()
}).T)
failed_cases = []
for system, result in quality_report['systems'].items():
    for case in result['cases']:
        failed = [name for name, passed in case['checks'].items() if not passed]
        if failed:
            failed_cases.append({
                'system': system, 'case': case['id'],
                'failed_checks': ', '.join(failed), 'decision': case['result']['status'],
            })
display(pd.DataFrame(failed_cases))
print((SYSTEM_ROOT / 'llm_system' / 'HUMAN_REVIEW_RUBRIC.md').read_text())
"""),
    md("""
## 11. Integrated Results and Biological Interpretation

The scientific result is the consensus module analysis derived from canonical DGE tables, with
detection thresholds, optimizer sensitivity, sampling stability, and null behavior beside each
interpretation. Planning can locate and propose approved operations, but cannot upgrade an
exploratory association into a confirmed mechanism. Modeling backends and the adapted router
verify software behavior only.

The internal audit preserves the breadth of the initial work. Every archived file has a
disposition and corrected-method note. A note is traceability, not proof that a replacement
result passed scientific validation.
"""),
    code("""
coverage = pd.read_csv(ROOT / 'output' / 'master' / 'legacy_output_coverage.csv')
assert coverage['legacy_path'].is_unique
assert coverage[['category', 'disposition', 'corrected_counterpart']].notna().all().all()
original_methods = {
    'detection-set summary': 'Thresholded DEG set summaries',
    'module contribution table': 'Single-fit module contribution profiles',
    'low-dimensional visualization': 'Raw module-percentage correlation geometry',
    'GO enrichment': 'Module-wise enrichment with legacy background and correction family',
    'GO presentation grouping': 'Post hoc functional label grouping',
    'sankey': 'Chemical-to-module-to-term flow display',
    'eigengene': 'Module eigengenes correlated across ten chemical conditions',
    'notebook': 'Initial executable analysis record',
    'generated output': 'Derived tables and figures from the initial workflow',
    'metadata': 'Archive and platform metadata',
}
identified_issues = {
    'detection-set summary': 'Thresholding was mixed with effect-size representation',
    'module contribution table': 'A single stochastic fit did not represent assignment uncertainty',
    'low-dimensional visualization': 'Unequal module sizes induced a high correlation null',
    'GO enrichment': 'Background, annotation version, and correction family were not fully controlled',
    'GO presentation grouping': 'Labels were assigned after inspection and need explicit provenance',
    'sankey': 'Presentation padding could imply unsupported flow magnitudes',
    'eigengene': 'Circular construction and ten-condition inference limit confirmatory use',
    'notebook': 'Resampling discarded subset labels and mixed exploratory with supported results',
    'generated output': 'Scientific status depended on the upstream method that generated each artifact',
    'metadata': 'Files do not constitute scientific outputs',
}
output_audit = (
    coverage.groupby(['category', 'disposition', 'corrected_counterpart'], dropna=False)
    .agg(original_outputs=('legacy_path', 'count'), example=('legacy_path', 'first'))
    .reset_index()
    .rename(columns={
        'category': 'Original output family',
        'disposition': 'Scientific status',
        'corrected_counterpart': 'Identified issue and corrected method',
        'example': 'Example original output',
    })
)
output_audit['Original method'] = output_audit['Original output family'].map(original_methods)
output_audit['Identified issue'] = output_audit['Original output family'].map(identified_issues)
output_audit = output_audit.rename(
    columns={'Identified issue and corrected method': 'Corrected method'}
)
output_audit['Final location'] = 'Executed master notebook and output/master audit artifacts'
display(output_audit[[
    'Original output family', 'Original method', 'Identified issue', 'Corrected method',
    'Final location', 'Scientific status', 'original_outputs', 'Example original output',
]])
print({'archived_entries': len(coverage), 'audited_output_families': len(output_audit)})
"""),
    md("""
## 12. Reproducibility and Execution Summary

This matrix records what can be regenerated, current evidence, and scope boundaries.
Infrastructure supports local containers and an AWS target, but this notebook neither provisions
resources nor incurs charges. Fast checks cover schemas, scientific helpers, both model backends,
retrieval, unsupported-claim handling, and notebook construction. Expensive full-data fitting
and model adaptation remain explicit jobs.
"""),
    code("""
reproducibility = pd.DataFrame([
    {'component': 'Canonical DGE inputs', 'evidence': 'SHA-256 and 13,852-row assertions',
     'status': 'verified during execution', 'scope_boundary': 'Raw counts and sample metadata unavailable'},
    {'component': 'Consensus modules', 'evidence': 'Thirty aligned fixed-seed fits plus per-gene agreement',
     'status': 'computational replay and inferential stability reported separately',
     'scope_boundary': 'No external biological replication'},
    {'component': 'Backend consistency', 'evidence': str(parity_path.relative_to(PROJECT)),
     'status': 'both backend runs regenerated in this execution', 'scope_boundary': 'Synthetic requests only'},
    {'component': 'Grounded planning', 'evidence': 'Deterministic evaluation in Section 10',
     'status': 'acceptance and abstention evaluated', 'scope_boundary': 'Human approval required'},
    {'component': 'Adapted router', 'evidence': 'Manifest, split hashes, model card, and metrics',
     'status': 'bounded artifact loaded when present', 'scope_boundary': 'Loss is not workflow validity'},
    {'component': 'Packaged execution', 'evidence': 'Dockerfile and CI configuration inspected outside this kernel',
     'status': 'configuration present; runtime evidence is environment-specific',
     'scope_boundary': 'Deployment requires approval and local Docker requires a running daemon'},
])
display(reproducibility)
print('Regenerate: python3 build_master_notebook.py')
print('Execute: jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=kspaces_venv PFAS_master_analysis.ipynb')
"""),
])

notebook = nbf.v4.new_notebook(cells=cells)
for index, cell in enumerate(notebook.cells):
    cell["id"] = hashlib.sha1(
        f"{index}:{cell.cell_type}:{cell.source}".encode()
    ).hexdigest()[:12]
notebook.metadata.kernelspec = canonical.metadata.get("kernelspec", {})
notebook.metadata.language_info = canonical.metadata.get("language_info", {})
nbf.write(notebook, OUTPUT)
if generated_intermediate:
    CANONICAL.unlink()
print(f"Wrote {OUTPUT} with {len(cells)} cells")
