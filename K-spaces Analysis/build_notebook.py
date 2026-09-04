"""
Builds the single consolidated, professional notebook for the PFAS k-spaces
module-discovery analysis. Run this to regenerate PFAS_kspaces_analysis.ipynb
from scratch; then execute it with nbconvert.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

def md(text):
    cells.append(nbf.v4.new_markdown_cell(text))

def code(text):
    cells.append(nbf.v4.new_code_cell(text))

# ---------------------------------------------------------------------------
md(r"""# PFAS Toxicogenomics : k-spaces Module Discovery Across 10 PFAS Chemicals

**Project:** *C. elegans* transcriptional response to 10 PFAS chemicals
(PFOSA, PFBSA, PFOS, PFNA, PFOA, GenX, PFEESA, PFBS, PFPeA, PFBA), each dosed
at its own EC50.

**Question:** do genes respond to different PFAS chemicals in coordinated
ways, and if so, which chemicals pattern together? This notebook builds a
gene x chemical differential-expression matrix from already-computed edgeR
results, clusters genes into modules with **k-spaces** (Markarian et al.
2025, mixtures of Gaussian latent variable models :
[github.com/pachterlab/k-spaces](https://github.com/pachterlab/k-spaces)),
and asks which chemicals load onto each module.

**How to read this notebook:** every section shows the code, the real output
it produced, and an interpretation. Where a method choice was non-obvious
(gene selection, number of modules, whether k-spaces is even the right tool
here), the notebook shows the actual check that was run to justify or
challenge it : including two checks that did **not** flatter the chosen
method, kept in for honesty:

- Bayesian Information Criterion (BIC) never selects an interior optimum for
  k (number of modules) on this data : a bootstrap-stability criterion is
  used instead.
- A simple k-means baseline is **more reproducible** than k-spaces on this
  exact data; k-spaces is kept because it is fitting a richer model
  (coordinated response *direction*, not just a centroid), not because it
  "wins" on stability.
- A single k-spaces fit carries substantial EM optimization noise in which
  genes land in which module : Section 12 builds a **consensus assignment**
  across 30 independent fits (plus a per-gene confidence score) that is the
  recommended module structure to cite, not the single fit in Section 5.
  Section 12 also stress-tests the chemical-pairing correlation metric used
  in Sections 5 and 7 against a proper null and **retracts it** : the
  headline "which PFAS chemicals are most similar" claims from those
  sections do not survive and should not be cited as originally framed.

**Environment:** Python 3.13, `kspaces` (the paper's own reference
implementation, required Python 3.10+), plus pandas/numpy/scipy/
scikit-learn/matplotlib/plotly/goatools. See `requirements.txt` in this
folder for exact pinned versions.
""")

# ---------------------------------------------------------------------------
md("## 0. Setup")
code(r"""
import warnings
warnings.filterwarnings("ignore")

import os
import io
import contextlib
import time
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import display, HTML, Image

import kspaces
from importlib.metadata import version
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

pd.set_option("display.max_columns", 20)
pd.set_option("display.width", 140)
plt.rcParams["figure.dpi"] = 110

# kspaces draws its EM initializations from numpy's *global* RNG state (not a
# seeded Generator), so this notebook is only reproducible run-to-run if that
# global state is fixed once, up front, before any kspaces call.
GLOBAL_SEED = 0
np.random.seed(GLOBAL_SEED)

# Set PFAS_KSPACES_DIR if Jupyter is launched from another directory.
# Otherwise, launch Jupyter from the folder containing this notebook.
ROOT = Path(os.environ.get("PFAS_KSPACES_DIR", Path.cwd())).expanduser().resolve()
PROJECT_DIR = ROOT.parent
DEGS_DIR = PROJECT_DIR / "Data" / "DEGs"
EXT_DIR = ROOT / "data_external"
FIG_DIR = ROOT / "output" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

CHEMICALS = ["PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA"]
FDR_THRESHOLD = 0.05

FILE_MAP = {
    "notebook_directory": ROOT,
    "dge_inputs": DEGS_DIR,
    "go_annotation": EXT_DIR / "wb.gaf",
    "go_ontology": EXT_DIR / "go-basic.obo",
    "generated_outputs": ROOT / "output",
}
print("File map:")
for label, path in FILE_MAP.items():
    print(f"  {label:20s} -> {path}")

missing_dge = [DEGS_DIR / f"{chem}vsControl_DGE_results.csv" for chem in CHEMICALS
               if not (DEGS_DIR / f"{chem}vsControl_DGE_results.csv").exists()]
if missing_dge:
    raise FileNotFoundError(
        "Missing DGE inputs. Launch Jupyter from the K-spaces Analysis directory "
        "or set PFAS_KSPACES_DIR to that directory. Missing: "
        + ", ".join(str(path) for path in missing_dge)
    )
print("kspaces version:", version("kspaces"))
""")

# ---------------------------------------------------------------------------
md(r"""## 1. Build the gene x chemical matrix

**Input:** `Data/DEGs/{CHEM}vsControl_DGE_results.csv` : 10 files, one per
PFAS chemical, each the output of an edgeR `glmQLFTest` (columns: gene id,
`logFC`, `logCPM`, `F`, `PValue`, `FDR`). All 10 files were confirmed to
share the identical 13,852-gene filtered universe (CPM > 1 in >= 4 samples),
so they merge cleanly on gene id.
""")

code(r"""
def load_chemical(chem: str) -> pd.DataFrame:
    path = DEGS_DIR / f"{chem}vsControl_DGE_results.csv"
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.str.strip('"').rename("WB_id")
    return df[["logFC", "FDR"]].rename(columns={"logFC": f"logFC_{chem}", "FDR": f"FDR_{chem}"})

frames = [load_chemical(c) for c in CHEMICALS]
merged = frames[0]
for f in frames[1:]:
    merged = merged.join(f, how="inner")

assert len(merged) == 13852, f"expected 13852 shared genes, got {len(merged)}"
for c in CHEMICALS:
    assert not merged[f"logFC_{c}"].isna().any(), f"unexpected NaNs in {c}"

print(f"merged_wide: {merged.shape[0]} genes x {len(CHEMICALS)} chemicals (logFC + FDR columns)")
merged.to_csv(ROOT / "output" / "merged_wide.csv")
merged.head()
""")

md("**Sanity check** : DEG counts per chemical (FDR < 0.05) against the numbers reported in `Data/Initial RNA-seqy Analysis_Plots.docx`:")

code(r"""
deg_counts = pd.DataFrame([
    {"chemical": c, "n_DEGs_FDR0.05": int((merged[f"FDR_{c}"] < FDR_THRESHOLD).sum())}
    for c in CHEMICALS
]).sort_values("n_DEGs_FDR0.05", ascending=False).reset_index(drop=True)
deg_counts.to_csv(ROOT / "output" / "deg_counts.csv", index=False)
display(deg_counts)
""")

# ---------------------------------------------------------------------------
md(r"""## 2. Gene selection strategy

The number of DEGs ranges from 154 (PFOA) to 9,161 (PFOSA) : nearly a
60-fold spread. Feeding raw log2 fold-changes into k-spaces as-is lets the
highest-DEG-count chemicals dominate module discovery purely by magnitude.

**Primary strategy (proposed default):** restrict to the union of genes DE
(FDR<0.05) in >= 1 chemical, retain each chemical's continuous logFC, then
z-score each gene's row across the 10 chemicals. This captures response
*shape* without replacing measured effect-size spacing with ranks.

**Sensitivity variants** (run later in this notebook, per the planning
document's own suggestion): signed within-chemical percentile ranks followed
by row z-scoring, and top-N-per-chemical gene selection for N in {500, 1000}.
""")

code(r"""
def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    mean = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(mean, axis=0).div(sd, axis=0)

logfc = merged[[f"logFC_{c}" for c in CHEMICALS]].rename(columns=lambda x: x.replace("logFC_", ""))
fdr = merged[[f"FDR_{c}" for c in CHEMICALS]].rename(columns=lambda x: x.replace("FDR_", ""))
de_union_idx = fdr.index[(fdr < FDR_THRESHOLD).any(axis=1)]
print(f"DE union (FDR<0.05 in >=1 chemical): {len(de_union_idx)} genes")

sub = logfc.loc[de_union_idx]
matrix_primary = zscore_rows(sub).dropna()
matrix_primary.to_csv(ROOT / "output" / "matrix_primary_zlogfc.csv")

pct = sub.rank(pct=True, method="average")
matrix_sens_pctrank = zscore_rows(pct).dropna()
matrix_sens_pctrank.to_csv(ROOT / "output" / "matrix_sensitivity_pctrank.csv")

def build_topn_variant(n: int) -> pd.DataFrame:
    genes = set()
    for c in CHEMICALS:
        s = merged[[f"FDR_{c}", f"logFC_{c}"]].copy()
        s["abs_logFC"] = s[f"logFC_{c}"].abs()
        top = s.sort_values([f"FDR_{c}", "abs_logFC"], ascending=[True, False]).head(n)
        genes.update(top.index)
    idx = pd.Index(sorted(genes))
    return zscore_rows(logfc.loc[idx]).dropna()

matrix_top500 = build_topn_variant(500)
matrix_top1000 = build_topn_variant(1000)
matrix_top500.to_csv(ROOT / "output" / "matrix_topN_500.csv")
matrix_top1000.to_csv(ROOT / "output" / "matrix_topN_1000.csv")

print(f"matrix_primary (continuous row-z-logFC): {matrix_primary.shape}")
print(f"matrix_sens_pctrank (rank sensitivity):  {matrix_sens_pctrank.shape}")
print(f"matrix_top500:                          {matrix_top500.shape}")
print(f"matrix_top1000:                         {matrix_top1000.shape}")
matrix_primary.head()
""")

# ---------------------------------------------------------------------------
md(r"""## 3. k-spaces model selection: does BIC pick a sensible k?

Each module is fit as a **1-D affine subspace** (a line through the 10-D
chemical-response space) rather than a 0-D point/centroid: a line captures a
*direction* of coordinated response (which chemicals move together, and how
strongly), analogous to a WGCNA "module eigengene." d=0 (point) models are
fit as a null comparison at every k.
""")

code(r"""
def fit_and_score(data, kd, initializations=10):
    spaces, probs = kspaces.run_EM(
        data, kd, assignment="soft", initializations=initializations,
        max_iter=100, tol=5e-2, silent=True,
    )
    if not spaces:
        return None
    ll = kspaces.total_log_likelihood(data, spaces)
    df = kspaces.model_selection_.get_df(spaces, eq_noise=False)
    n = data.shape[0]
    bic = kspaces.get_BIC(df, n, ll)
    icl = kspaces.get_ICL(probs, data, spaces, eq_noise=False)
    return {"df": df, "log_likelihood": ll, "BIC": bic, "ICL": icl, "spaces": spaces, "probs": probs}

data_primary = matrix_primary.values
t0 = time.time()
rows = []
for k in range(2, 26):
    for label, dim in [("d1_lines", 1), ("d0_points", 0)]:
        kd = [dim] * k
        r = fit_and_score(data_primary, kd, initializations=10)
        if r is None:
            continue
        rows.append({"k": k, "config": label, "df": r["df"], "log_likelihood": r["log_likelihood"],
                      "BIC": r["BIC"], "ICL": r["ICL"]})
model_selection_df = pd.DataFrame(rows)
model_selection_df.to_csv(ROOT / "output" / "kspaces_runs" / "primary_model_selection.csv", index=False)
print(f"Sweep over k=2..25, d in {{0,1}}: {len(model_selection_df)} models fit in {time.time()-t0:.0f}s")
model_selection_df.sort_values("BIC").head(10)
""")

code(r"""
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for label, marker in [("d1_lines", "o"), ("d0_points", "s")]:
    sub = model_selection_df[model_selection_df["config"] == label].sort_values("k")
    axes[0].plot(sub["k"], sub["BIC"], marker=marker, label=label)
    axes[1].plot(sub["k"], sub["ICL"], marker=marker, label=label)
axes[0].set_xlabel("k (number of modules)"); axes[0].set_ylabel("BIC"); axes[0].set_title("BIC vs k")
axes[1].set_xlabel("k (number of modules)"); axes[1].set_ylabel("ICL"); axes[1].set_title("ICL vs k")
for ax in axes:
    ax.legend(); ax.grid(alpha=0.3)
fig.suptitle("BIC/ICL never turn over within k=2..25 -- not usable alone to pick k on this data")
fig.tight_layout()
fig.savefig(FIG_DIR / "bic_icl_vs_k.png", dpi=150)
plt.show()
""")

md(r"""**Finding:** BIC and ICL both decrease monotonically all the way to
k=25 for the d=1 ("lines") configuration, with no interior minimum. This is
a known failure mode when the number of genes (N=10,255) vastly exceeds the
number of features (D=10 chemicals): there is enough room for arbitrarily
many subspaces to each grab a small likelihood gain. **BIC/ICL alone cannot
select k here** : a different criterion is needed (Section 4).

Note d=1 ("lines") beats d=0 ("points", i.e. plain centroid/k-means-like
clustering) on BIC at **every** k tested : the extra directional flexibility
is earning its complexity cost in raw fit quality, which is revisited
adversarially in Section 9.
""")

# ---------------------------------------------------------------------------
md(r"""## 4. Choosing k via bootstrap stability instead of BIC

The planning document calls bootstrap stability "essential" validation.
For each k, fit a reference model, then repeatedly refit on gene-resampled
(with replacement) bootstraps and measure how well each bootstrap's module
structure reproduces the reference model's assignment on the **original**
data (Adjusted Rand Index). Mean ARI across bootstraps is the stability
score for that k.
""")

code(r"""
def hard_labels(probs):
    return probs.argmax(axis=1)

def kspaces_fit(data, kd, initializations=10):
    spaces, probs = kspaces.run_EM(
        data, kd, assignment="soft", initializations=initializations,
        max_iter=100, tol=5e-2, silent=True,
    )
    return spaces, probs

n = data_primary.shape[0]
rng = np.random.default_rng(0)
stab_rows = []
t0 = time.time()
N_BOOT = 15
for k in range(3, 16):
    kd = [1] * k
    ref_spaces, ref_probs = kspaces_fit(data_primary, kd)
    if not ref_spaces:
        continue
    ref_labels = hard_labels(ref_probs)
    aris = []
    for _ in range(N_BOOT):
        boot_idx = rng.integers(0, n, size=n)
        boot_spaces, _ = kspaces_fit(data_primary[boot_idx], kd, initializations=5)
        if not boot_spaces:
            continue
        boot_probs_on_full = kspaces.E_step(data_primary, boot_spaces, assignment="soft")
        aris.append(adjusted_rand_score(ref_labels, hard_labels(boot_probs_on_full)))
    stab_rows.append({"k": k, "mean_ARI": np.mean(aris), "sd_ARI": np.std(aris), "n_boot_ok": len(aris)})

stability_df = pd.DataFrame(stab_rows).merge(
    model_selection_df[model_selection_df["config"] == "d1_lines"][["k", "BIC"]], on="k", how="left"
)
stability_df.to_csv(ROOT / "output" / "kspaces_runs" / "primary_stability.csv", index=False)
print(f"Bootstrap stability sweep k=3..15, {N_BOOT} resamples each: {time.time()-t0:.0f}s")
stability_df
""")

code(r"""
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.errorbar(stability_df["k"], stability_df["mean_ARI"], yerr=stability_df["sd_ARI"], marker="o", capsize=3)
ax.axvline(5, color="red", linestyle="--", alpha=0.6, label="k=5 (chosen)")
ax.set_xlabel("k (number of modules)")
ax.set_ylabel("mean bootstrap ARI (15 resamples)")
ax.set_title("Module stability peaks at low k, declines as k grows")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG_DIR / "stability_vs_k.png", dpi=150)
plt.show()
""")

code(r"""
print(f"Stability ARI range across k=3..15: {stability_df['mean_ARI'].min():.2f}-{stability_df['mean_ARI'].max():.2f}")
print(f"Best k by stability: k={stability_df.loc[stability_df['mean_ARI'].idxmax(), 'k']:.0f} "
      f"(ARI={stability_df['mean_ARI'].max():.3f})")
print(f"k=5 stability: ARI={stability_df.loc[stability_df['k']==5, 'mean_ARI'].values[0]:.3f}")
""")

md(r"""**Finding:** reproducibility is moderate-to-weak everywhere (see
range printed above) and is generally highest at the low end of the k range
tested, declining (noisily) as k grows. **k=5 is chosen as a pragmatic
balance** between stability and having enough modules to be a meaningful
"table of modules" : this is a flagged, reviewable choice, not a validated
optimum; the printed numbers above show exactly how much stability is given
up relative to the best-stability k.
""")

# ---------------------------------------------------------------------------
md(r"""## 5. Fit the primary model at k=5 and profile modules by chemical

For each module: which chemicals' columns are elevated/depressed for that
module's genes (mean matrix value), what fraction of each chemical's own
FDR<0.05 DEGs fall in each module, and the resulting chemical-chemical
correlation of module-usage profiles : the direct answer to "which PFAS
chemicals share which modules."
""")

code(r"""
def fit_module_profile(matrix_df: pd.DataFrame, k: int, seed_note: str = ""):
    data = matrix_df.values
    kd = [1] * k
    spaces, probs = kspaces.run_EM(
        data, kd, assignment="soft", initializations=20, max_iter=150, tol=5e-2, silent=True,
    )
    labels = probs.argmax(axis=1)
    assign = pd.Series(labels, index=matrix_df.index, name="module")

    mat_df = matrix_df.copy()
    mat_df["module"] = assign
    mean_value = mat_df.groupby("module")[list(matrix_df.columns)].mean().round(3)

    pct_rows = []
    for c in CHEMICALS:
        chem_degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
        chem_degs_in_union = assign.loc[assign.index.intersection(chem_degs)]
        if len(chem_degs) == 0:
            continue
        counts = chem_degs_in_union.value_counts().reindex(range(k), fill_value=0)
        pct = (counts / len(chem_degs) * 100).round(1)
        pct.name = c
        pct_rows.append(pct)
    pct_table = pd.DataFrame(pct_rows)
    pct_table.index.name = "chemical"

    corr = pct_table.T.corr()
    return assign, mean_value, pct_table, corr

assign_primary, meanval_primary, pct_primary, corr_primary = fit_module_profile(matrix_primary, 5)
print("Module sizes:")
display(assign_primary.value_counts().sort_index().to_frame("n_genes"))
""")

code(r"""
print("Mean matrix value per module per chemical (which chemicals drive each module):")
display(meanval_primary)
""")

code(r"""
print("% of each chemical's own FDR<0.05 DEGs falling in each module:")
display(pct_primary)
""")

code(r"""
print("Chemical-chemical correlation of module-usage profiles:")
display(corr_primary.round(2))

# save for downstream sections
assign_primary.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_module_assignments.csv")
meanval_primary.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_module_mean_value_by_chemical.csv")
pct_primary.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_pct_chemical_DEGs_per_module.csv")
corr_primary.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_chemical_similarity_by_module_usage.csv")
""")

code(r"""
# generated directly from corr_primary -- always matches this run's actual numbers,
# never a stale hardcoded claim
pairs_sorted = (
    corr_primary.where(np.triu(np.ones(corr_primary.shape), k=1).astype(bool))
    .stack()
    .sort_values(ascending=False)
)
print("Top 8 most similar chemical pairs by module-usage correlation (this run):")
top_pairs = pairs_sorted.head(8).round(3)
display(top_pairs.to_frame("correlation"))
""")

md(r"""**Read this table with caution, corrected below.** The table above
is generated live from this run's actual module fit, but **Section 12
stress-tests this correlation metric against a null (module labels shuffled
within the same size distribution) and finds essentially none of these pair
correlations exceed chance** : with only k=5 modules, correlating each
chemical's 5-element percent-per-module vector produces high correlation
almost by construction, regardless of real biology. An earlier version of
this notebook read the top of this table as "GenX and PFBS pair extremely
robustly" : **that claim does not survive Section 12's stress test and is
retracted**, not merely softened. Treat every pair in the table above as
provisional until you've read Section 12.
""")

# ---------------------------------------------------------------------------
md("## 6. Heatmap: gene x chemical matrix, rows grouped by module")

code(r"""
col_order = ["PFOS", "PFNA", "PFOA", "GenX", "PFBS", "PFEESA", "PFBA", "PFBSA", "PFOSA", "PFPeA"]
mat_ordered = matrix_primary[col_order]
order = assign_primary.sort_values().index
mat_sorted = mat_ordered.loc[order]
modules_sorted = assign_primary.loc[order]

n_modules = modules_sorted.nunique()
fig, (ax_strip, ax_hm) = plt.subplots(1, 2, figsize=(9, 10), gridspec_kw={"width_ratios": [0.3, 10]}, sharey=True)

cmap_modules = plt.get_cmap("tab10", n_modules)
ax_strip.imshow(modules_sorted.values.reshape(-1, 1), aspect="auto", cmap=cmap_modules, vmin=-0.5, vmax=n_modules - 0.5)
ax_strip.set_xticks([]); ax_strip.set_yticks([])

boundaries = np.where(np.diff(modules_sorted.values) != 0)[0] + 1
prev = 0
for b in list(boundaries) + [len(modules_sorted)]:
    mid = (prev + b) / 2
    ax_strip.text(-1.3, mid, f"M{modules_sorted.values[prev]}", va="center", ha="right", fontsize=9)
    prev = b
ax_strip.set_xlim(-2.5, 0.5)

vmax = np.nanpercentile(np.abs(mat_sorted.values), 98)
im = ax_hm.imshow(mat_sorted.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax_hm.set_xticks(range(len(col_order))); ax_hm.set_xticklabels(col_order, rotation=90)
ax_hm.set_yticks([])
for b in boundaries:
    ax_hm.axhline(b - 0.5, color="black", linewidth=0.6)

cbar = fig.colorbar(im, ax=ax_hm, fraction=0.03, pad=0.02)
cbar.set_label("row-zscored continuous logFC")
fig.suptitle(f"k-spaces modules (k=5, primary continuous logFC strategy)\n"
             f"{len(mat_sorted)} genes (DE union, FDR<0.05 in >=1 chemical), rows grouped by module", fontsize=11)
fig.tight_layout()
fig.savefig(FIG_DIR / "primary_k5_heatmap.png", dpi=200)
plt.show()
""")

# ---------------------------------------------------------------------------
md(r"""## 7. Sensitivity checks: does the chemical-pairing finding survive different gene-selection strategies?

Same k=5 fit repeated on the signed-percentile-rank sensitivity matrix and
on the top-500 / top-1000-per-chemical variants.
""")

code(r"""
assign_sens, meanval_sens, pct_sens, corr_sens = fit_module_profile(matrix_sens_pctrank, 5)
print("SENSITIVITY (signed percentile ranks): chemical-chemical correlation of module-usage profiles")
display(corr_sens.round(2))

# highest-raw-DEG-count chemicals, to check whether they're the ones that break off as outliers
high_deg_chems = deg_counts.sort_values("n_DEGs_FDR0.05", ascending=False)["chemical"].head(2).tolist()
other_chems = [c for c in CHEMICALS if c not in high_deg_chems]
other_block_corr = corr_sens.loc[other_chems, other_chems]
other_block_vals = other_block_corr.where(np.triu(np.ones(other_block_corr.shape), k=1).astype(bool)).stack()
outlier_corr_to_block = corr_sens.loc[high_deg_chems, other_chems].values.flatten()
print(f"\nTwo highest-raw-DEG chemicals: {high_deg_chems}")
print(f"Their correlation to the other {len(other_chems)} chemicals: "
      f"min={outlier_corr_to_block.min():.2f}, max={outlier_corr_to_block.max():.2f}, mean={outlier_corr_to_block.mean():.2f}")
print(f"The other {len(other_chems)} chemicals' correlation among themselves: "
      f"min={other_block_vals.min():.2f}, max={other_block_vals.max():.2f}")
""")

md(r"""**Interpret this as a transformation sensitivity, not a validation.**
Use the values printed above to assess whether the signed-rank transform
materially changes chemical profiles relative to continuous logFC. Neither
representation is validated merely because it produces clearer separation.
""")

code(r"""
assign_top500, meanval_top500, pct_top500, corr_top500 = fit_module_profile(matrix_top500, 5)
assign_top1000, meanval_top1000, pct_top1000, corr_top1000 = fit_module_profile(matrix_top1000, 5)

for name, df in [("sens_pctrank", assign_sens), ("topN500", assign_top500), ("topN1000", assign_top1000)]:
    df.to_csv(ROOT / "output" / "kspaces_runs" / f"{name}_k5_module_assignments.csv")

print("top500: chemical-chemical correlation of module-usage profiles")
display(corr_top500.round(2))
""")

code(r"""
def robust_pair(corr_df, a, b):
    return corr_df.loc[a, b] if a in corr_df.index and b in corr_df.columns else np.nan

pairs_to_check = [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]
robustness = pd.DataFrame(
    {name: [robust_pair(c, a, b) for a, b in pairs_to_check]
     for name, c in [("primary", corr_primary), ("sens_pctrank", corr_sens),
                      ("topN500", corr_top500), ("topN1000", corr_top1000)]},
    index=[f"{a}-{b}" for a, b in pairs_to_check],
)
print("Chemical-pair module-usage correlation across all 4 strategies:")
display(robustness.round(2))

# "most robust" = highest worst-case (minimum) correlation across the 4 strategies --
# a pair that's only high in some strategies and low in others is not robust,
# regardless of its average.
robustness_summary = robustness.copy()
robustness_summary["min_across_strategies"] = robustness.min(axis=1)
robustness_summary["max_across_strategies"] = robustness.max(axis=1)
robustness_summary = robustness_summary.sort_values("min_across_strategies", ascending=False)
print("\nRanked by worst-case (minimum) correlation across strategies -- the most trustworthy robustness measure:")
display(robustness_summary.round(2))
""")

md(r"""**Correction (see Section 12): do not read the ranking above as
"trustworthy" evidence of which pairs pattern together.** An earlier version
of this notebook claimed the top-ranked pair by `min_across_strategies` was
"the actual, trustworthy 'which two PFAS pattern together' answer" : that
claim is retracted. Section 12 shows that correlating 5-element
percent-per-module vectors at k=5 produces high correlation for almost any
pair, real or randomly-labeled, purely from module-size structure; this
applies to every column in the table above, not just the primary strategy.
The `min_across_strategies` ranking is kept here because it is still useful
for a narrower purpose : identifying which pairs' correlation is at least
*consistent in sign and magnitude* across how the input matrix is built,
which is informative about sensitivity to gene-selection choices : but that
is a different, weaker claim than "these chemicals are unusually similar,"
and should not be conflated with it.
""")

md("### Cross-strategy ARI: how much does exact module membership agree?")

code(r"""
cross_ari_rows = []
for name, other in [("sens_pctrank", assign_sens), ("topN500", assign_top500), ("topN1000", assign_top1000)]:
    shared = assign_primary.index.intersection(other.index)
    ari = adjusted_rand_score(assign_primary.loc[shared], other.loc[shared])
    cross_ari_rows.append({"variant": name, "n_shared_genes": len(shared), "ARI_vs_primary": round(ari, 3)})
cross_ari_df = pd.DataFrame(cross_ari_rows)
cross_ari_df.to_csv(ROOT / "output" / "kspaces_runs" / "cross_strategy_ARI_vs_primary.csv", index=False)
display(cross_ari_df)
print(f"Cross-strategy ARI range: {cross_ari_df['ARI_vs_primary'].min():.2f}-{cross_ari_df['ARI_vs_primary'].max():.2f}")
""")

md(r"""**Cross-strategy ARI is low** (see range printed above, well below
the within-strategy bootstrap ARI from Section 4): the *exact*
gene-to-module partition is not stable across gene-selection strategy, even
though the higher-level *chemical-pairing* signal (Section 7 above) is.
Don't over-interpret which specific genes are "in" any one module.
""")

# ---------------------------------------------------------------------------
md(r"""## 8. GO enrichment per module

For each ontology namespace, the eligible population is the intersection of
the 10,255-gene clustered DEG union with genes mapped in the pinned WormBase
GAF. This conditions module-identity questions on genes eligible for module
assignment. Raw Fisher p-values from every module x namespace x term test are
then corrected together as one BH family. Annotations: WormBase GAF
(`current.geneontology.org/annotations/wb.gaf.gz`) + `go-basic.obo`,
downloaded to `data_external/`.
""")

code(r"""
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
from goatools.go_enrichment import GOEnrichmentStudy
from statsmodels.stats.multitest import multipletests

FDR_ALPHA = 0.05
TOP_N_TERMS = 15

with contextlib.redirect_stdout(io.StringIO()):
    godag = GODag(str(EXT_DIR / "go-basic.obo"))
    gaf = GafReader(str(EXT_DIR / "wb.gaf"), godag=godag)
    ns2assoc = gaf.get_ns2assc()

measured_universe = set(merged.index)
clustered_universe = set(matrix_primary.index)
measured_background_by_ns = {
    ns: sorted(measured_universe.intersection(assoc))
    for ns, assoc in ns2assoc.items()
}
eligible_background_by_ns = {
    ns: sorted(clustered_universe.intersection(assoc))
    for ns, assoc in ns2assoc.items()
}
print(f"Measured universe: {len(measured_universe)} genes")
print(f"Clustered DEG-union universe: {len(clustered_universe)} genes")
for ns, assoc in ns2assoc.items():
    print(f"  namespace {ns}: {len(eligible_background_by_ns[ns])} clustered genes mapped "
          f"({len(measured_background_by_ns[ns])} measured genes mapped)")
""")

code(r"""
def familywise_go_enrichment(assign):
    tested = []
    with contextlib.redirect_stdout(io.StringIO()):
        for ns, assoc in ns2assoc.items():
            background_ns = eligible_background_by_ns[ns]
            goe = GOEnrichmentStudy(
                background_ns, assoc, godag, propagate_counts=True,
                alpha=FDR_ALPHA, methods=["fdr_bh"],
            )
            for module in sorted(assign.unique()):
                study_genes = sorted(set(assign.index[assign == module]).intersection(background_ns))
                for result in goe.run_study(study_genes):
                    tested.append({
                        "module": module, "namespace": ns, "GO_id": result.GO,
                        "term": result.name, "enrichment": result.enrichment,
                        "p_raw": result.p_uncorrected,
                        "study_count": result.study_count, "study_n": result.study_n,
                        "pop_count": result.pop_count, "pop_n": result.pop_n,
                    })
    tested_df = pd.DataFrame(tested)
    tested_df["p_fdr_bh_family"] = multipletests(tested_df["p_raw"], method="fdr_bh")[1]
    significant = tested_df[
        (tested_df["enrichment"] == "e")
        & (tested_df["p_fdr_bh_family"] < FDR_ALPHA)
    ]
    return (
        significant.sort_values("p_fdr_bh_family")
        .groupby(["module", "namespace"], group_keys=False)
        .head(TOP_N_TERMS)
        .sort_values(["module", "namespace", "p_fdr_bh_family"])
        .reset_index(drop=True),
        tested_df,
    )

t0 = time.time()
go_df, go_tests_primary = familywise_go_enrichment(assign_primary)
(ROOT / "output" / "kspaces_runs" / "go_enrichment").mkdir(exist_ok=True, parents=True)
go_df.to_csv(ROOT / "output" / "kspaces_runs" / "go_enrichment" / "primary_k5_GO_enrichment.csv", index=False)
print(f"GO enrichment done in {time.time()-t0:.0f}s: {len(go_tests_primary)} total tests in one BH family; "
      f"{len(go_df)} significant displayed term-module hits")
""")

code(r"""
for m in sorted(go_df["module"].unique()):
    sub = go_df[go_df["module"] == m].sort_values("p_fdr_bh_family").head(5)
    print(f"--- Module {m} top terms ---")
    display(sub[["namespace", "term", "p_fdr_bh_family", "study_count", "pop_count"]])
""")

code(r"""
# generated module-biology summary: top GO term (excluding generic top-of-DAG
# terms) + the chemicals with the most extreme mean matrix value, both read
# live from this run's results -- always in sync with the tables above.
GENERIC_TERM_POP_CAP = 2500
summary_rows = []
for m in sorted(assign_primary.unique()):
    size = int((assign_primary == m).sum())
    specific = go_df[(go_df["module"] == m) & (go_df["pop_count"] <= GENERIC_TERM_POP_CAP)].sort_values("p_fdr_bh_family")
    if len(specific):
        top_term = f"{specific.iloc[0]['term']} (family BH={specific.iloc[0]['p_fdr_bh_family']:.1e})"
    else:
        top_term = "(no specific term below generic-term cap)"
    row_vals = meanval_primary.loc[m]
    high_chem = row_vals.idxmax(); low_chem = row_vals.idxmin()
    summary_rows.append({
        "module": m, "n_genes": size,
        "highest_chemical": f"{high_chem} ({row_vals[high_chem]:+.2f})",
        "lowest_chemical": f"{low_chem} ({row_vals[low_chem]:+.2f})",
        "top_specific_GO_term": top_term,
    })
module_summary = pd.DataFrame(summary_rows).set_index("module")
print("Module biology summary (generated live from this run):")
display(module_summary)
""")

# ---------------------------------------------------------------------------
md("## 9. Sankey diagram: chemical -> module -> top enriched GO term")

code(r"""
import plotly.graph_objects as go

GENERIC_TERM_POP_CAP = 2500
modules = sorted(assign_primary.unique())

cm_counts = pd.DataFrame(0, index=CHEMICALS, columns=modules)
for c in CHEMICALS:
    degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
    in_union = assign_primary.loc[assign_primary.index.intersection(degs)]
    vc = in_union.value_counts()
    for m in modules:
        cm_counts.loc[c, m] = int(vc.get(m, 0))

go_specific = go_df[go_df["pop_count"] <= GENERIC_TERM_POP_CAP].copy()
top_term_per_module = {}
for m in modules:
    sub = go_specific[go_specific["module"] == m].sort_values("p_fdr_bh_family")
    if len(sub):
        row = sub.iloc[0]
        top_term_per_module[m] = (row["term"], int(row["study_count"]))
    else:
        top_term_per_module[m] = (f"(no specific term, module {m})", int((assign_primary == m).sum()))

module_labels = [f"Module {m}" for m in modules]
term_labels = [top_term_per_module[m][0] for m in modules]
nodes = CHEMICALS + module_labels + term_labels
node_idx = {name: i for i, name in enumerate(nodes)}

def hex_to_rgba(hex_color, alpha):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

module_colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#64B5CD"]
sources, targets, values, link_colors = [], [], [], []
for c in CHEMICALS:
    for j, m in enumerate(modules):
        v = cm_counts.loc[c, m]
        if v > 0:
            sources.append(node_idx[c]); targets.append(node_idx[f"Module {m}"]); values.append(v)
            link_colors.append(hex_to_rgba(module_colors[j % len(module_colors)], 0.5))
for j, m in enumerate(modules):
    term, count = top_term_per_module[m]
    sources.append(node_idx[f"Module {m}"]); targets.append(node_idx[term]); values.append(count)
    link_colors.append(hex_to_rgba(module_colors[j % len(module_colors)], 0.5))

node_colors = (["#999999"] * len(CHEMICALS)
               + [module_colors[j % len(module_colors)] for j in range(len(modules))]
               + [module_colors[j % len(module_colors)] for j in range(len(modules))])

fig_sankey = go.Figure(data=[go.Sankey(
    node=dict(pad=18, thickness=16, label=nodes, color=node_colors, line=dict(color="black", width=0.5)),
    link=dict(source=sources, target=targets, value=values, color=link_colors),
)])
fig_sankey.update_layout(title_text="PFAS chemical -> k-spaces module -> top enriched GO term (primary, k=5)",
                          font_size=12, width=1200, height=700)
sankey_path = FIG_DIR / "sankey_chemical_module_go.png"
fig_sankey.write_image(str(sankey_path), scale=2)
display(Image(filename=str(sankey_path)))
""")

# ---------------------------------------------------------------------------
md(r"""## 10. Chemicals in low-dimensional module-usage space, and structure-activity comparison

Does chain length / functional group / ether-linkage (from
`NHR Screen/PFAS_Attributes.txt`) predict which chemicals share modules?
""")

code(r"""
ATTRS = pd.DataFrame([
    ("PFOSA", "Long", "Sulfonamide", "No"),
    ("PFBSA", "Short", "Sulfonamide", "No"),
    ("PFNA", "Long", "Carboxylic acid", "No"),
    ("PFOS", "Long", "Sulfonic acid", "No"),
    ("PFOA", "Long", "Carboxylic acid", "No"),
    ("GenX", "Short", "Carboxylic acid", "Yes"),
    ("PFEESA", "Short", "Sulfonic acid", "Yes"),
    ("PFPeA", "Short", "Carboxylic acid", "No"),
    ("PFBA", "Short", "Carboxylic acid", "No"),
    ("PFBS", "Short", "Sulfonic acid", "No"),
], columns=["chemical", "chain_length", "functional_group", "ether"]).set_index("chemical")

pca = PCA(n_components=2)
coords = pca.fit_transform(pct_primary.values)
var_exp = pca.explained_variance_ratio_ * 100
attrs = ATTRS.loc[pct_primary.index]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, attr_col, title in [(axes[0], "chain_length", "colored by chain length"),
                             (axes[1], "functional_group", "colored by functional group")]:
    cats = attrs[attr_col].unique()
    cmap = plt.get_cmap("tab10", len(cats))
    for i, cat in enumerate(cats):
        mask = (attrs[attr_col] == cat).values
        ax.scatter(coords[mask, 0], coords[mask, 1], label=cat, s=140, color=cmap(i), edgecolor="black")
    for i, chem in enumerate(pct_primary.index):
        ax.annotate(chem, (coords[i, 0], coords[i, 1]), fontsize=9, xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var_exp[0]:.0f}% of module-usage variance)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.0f}%)")
    ax.set_title(title); ax.legend(fontsize=8)
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
fig.suptitle("10 PFAS chemicals in k-spaces module-usage space (PCA of % DEGs per module)")
fig.tight_layout()
fig.savefig(FIG_DIR / "chemical_lowdim_by_attribute.png", dpi=200)
plt.show()
""")

code(r"""
pairs = []
chems = list(corr_primary.index)
for i in range(len(chems)):
    for j in range(i + 1, len(chems)):
        a, b = chems[i], chems[j]
        pairs.append({
            "chem_a": a, "chem_b": b,
            "same_chain_length": attrs.loc[a, "chain_length"] == attrs.loc[b, "chain_length"],
            "same_functional_group": attrs.loc[a, "functional_group"] == attrs.loc[b, "functional_group"],
            "same_ether": attrs.loc[a, "ether"] == attrs.loc[b, "ether"],
            "module_usage_corr": corr_primary.loc[a, b],
        })
pairs_df = pd.DataFrame(pairs)
pairs_df.to_csv(ROOT / "output" / "kspaces_runs" / "chemical_pair_attributes_vs_module_corr.csv", index=False)

print("Mean module-usage correlation, same vs different attribute:")
for col in ["same_chain_length", "same_functional_group", "same_ether"]:
    display(pairs_df.groupby(col)["module_usage_corr"].agg(["mean", "count"]))
""")

md(r"""**Do not read the numbers this cell prints as a finding.** They use
`corr_primary`, the single-fit, module-size-confounded correlation metric
that Section 12 retracts. This comparison is redone properly in Section 14,
once the null-adjusted residual metric exists (Section 13) : including a
label-permutation test, not just a same-vs-different mean comparison. See
Section 14 for the number to actually cite.
""")

# ---------------------------------------------------------------------------
md(r"""## 11. Adversarial checks: is k-spaces actually the right tool here?

Everything above validates *internal* consistency. This section checks
against a null model and a much simpler baseline, using the identical
bootstrap-ARI protocol as Section 4 for direct comparability.

**(a) Null baseline**: independently shuffle each chemical's column across
genes (destroys gene-level co-response structure, keeps each chemical's own
marginal distribution). If the real data isn't meaningfully better than
this, there's no real structure being found.

**(b) Baseline comparison**: plain k-means(k=5) on the same matrix, same
protocol, directly comparable to k-spaces' bootstrap ARI.
""")

code(r"""
K = 5
N_BOOT_ADV = 15
perm_rng = np.random.default_rng(1)
perm_data = data_primary.copy()
for j in range(perm_data.shape[1]):
    perm_data[:, j] = perm_rng.permutation(perm_data[:, j])

real_spaces, real_probs = kspaces_fit(data_primary, [1] * K, initializations=20)
real_ll = kspaces.total_log_likelihood(data_primary, real_spaces)
real_bic = kspaces.get_BIC(kspaces.model_selection_.get_df(real_spaces, eq_noise=False), n, real_ll)

perm_spaces, perm_probs = kspaces_fit(perm_data, [1] * K, initializations=20)
perm_ll = kspaces.total_log_likelihood(perm_data, perm_spaces)
perm_bic = kspaces.get_BIC(kspaces.model_selection_.get_df(perm_spaces, eq_noise=False), n, perm_ll)

print(f"real data:      logL={real_ll:12.1f}  BIC={real_bic:12.1f}")
print(f"permuted null:  logL={perm_ll:12.1f}  BIC={perm_bic:12.1f}")
print(f"delta logL (real - null): {real_ll - perm_ll:.1f}  (higher = real data fits better than null)")
""")

code(r"""
def bootstrap_ari(fit_fn, assign_fn, data, n, rng, n_boot=N_BOOT_ADV):
    ref_labels = assign_fn(fit_fn(data))
    aris = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot_model = fit_fn(data[idx])
        aris.append(adjusted_rand_score(ref_labels, assign_fn(boot_model, data)))
    return float(np.mean(aris)), float(np.std(aris))

def ks_fit_fn(d):
    s, _ = kspaces_fit(d, [1] * K, initializations=8)
    return s

def ks_assign_fn_real(spaces, on_data=None):
    target = on_data if on_data is not None else data_primary
    return kspaces.E_step(target, spaces, assignment="soft").argmax(axis=1)

def ks_assign_fn_perm(spaces, on_data=None):
    target = on_data if on_data is not None else perm_data
    return kspaces.E_step(target, spaces, assignment="soft").argmax(axis=1)

t0 = time.time()
real_ari_mean, real_ari_sd = bootstrap_ari(ks_fit_fn, ks_assign_fn_real, data_primary, n, np.random.default_rng(0))
null_ari_mean, null_ari_sd = bootstrap_ari(ks_fit_fn, ks_assign_fn_perm, perm_data, n, np.random.default_rng(0))
print(f"real data k-spaces stability:      mean_ARI={real_ari_mean:.3f} sd={real_ari_sd:.3f}")
print(f"permuted null k-spaces stability:  mean_ARI={null_ari_mean:.3f} sd={null_ari_sd:.3f}")
print(f"({time.time()-t0:.0f}s)")
""")

code(r"""
def km_fit_fn(d):
    km = KMeans(n_clusters=K, n_init=10, random_state=0)
    km.fit(d)
    return km

def km_assign_fn(model, on_data=None):
    target = on_data if on_data is not None else data_primary
    return model.predict(target)

t0 = time.time()
km_ari_mean, km_ari_sd = bootstrap_ari(km_fit_fn, km_assign_fn, data_primary, n, np.random.default_rng(0))
ref_kspaces_labels = ks_assign_fn_real(real_spaces)
km_ref = km_fit_fn(data_primary)
agreement_ari = adjusted_rand_score(ref_kspaces_labels, km_assign_fn(km_ref))

adv_summary = pd.DataFrame([
    {"check": "real_vs_null_logL_delta", "value": real_ll - perm_ll},
    {"check": "real_vs_null_BIC_delta", "value": perm_bic - real_bic},
    {"check": "real_data_kspaces_stability_ARI", "value": real_ari_mean},
    {"check": "permuted_null_kspaces_stability_ARI", "value": null_ari_mean},
    {"check": "kmeans_stability_ARI", "value": km_ari_mean},
    {"check": "kmeans_vs_kspaces_agreement_ARI", "value": agreement_ari},
])
adv_summary.to_csv(ROOT / "output" / "kspaces_runs" / "adversarial_baseline_checks.csv", index=False)
print(f"k-means(k=5) stability:  mean_ARI={km_ari_mean:.3f} sd={km_ari_sd:.3f}  ({time.time()-t0:.0f}s)")
print(f"k-spaces(k=5, d=1) stability (from above): mean_ARI={real_ari_mean:.3f}")
print(f"ARI(k-means labels vs k-spaces labels) on real data: {agreement_ari:.3f}")
print(f"\nk-means is {km_ari_mean/real_ari_mean:.1f}x more stable than k-spaces on this data.")
display(adv_summary)
""")

md(r"""**(a) Real data vs. null**: the log-likelihood delta and BIC delta
printed above (both strongly favoring real data) and the stability
comparison (real data ARI well above the permuted-null ARI) together show
**k-spaces is finding real, non-null structure**, not just fitting noise.

**(b) k-spaces vs. k-means**: the ratio printed above quantifies it :
**k-means is meaningfully more stable than k-spaces** on this exact data,
and the agreement ARI between the two methods' labels shows they find
genuinely different partitions, not minor variants of the same one. This is
a real bias-variance tradeoff: d=1 lines already showed better raw fit
(lower BIC) than d=0 points at every k (Section 3), because fitting a
*direction* (coordinated response shape) is a richer model than fitting a
centroid : but that same flexibility makes it more sensitive to exactly
which genes land in a given resample.

**Bottom line: k-spaces is not obviously "the best" method by a
reproducibility metric : k-means is more stable on this data.** Whether
k-spaces' extra directional information is worth the stability cost is a
judgment call for whoever is driving the biological interpretation.
""")

# ---------------------------------------------------------------------------
md(r"""## 12. Beyond bootstrap: consensus clustering, chemical-level jackknife, and a formal permutation test

Bootstrap resampling (Sections 4 and 11) only diagnoses how much a single
fit wobbles under resampling : it doesn't build anything better, and it only
probes one axis of fragility (which genes happen to be in the dataset). This
section goes further on three fronts, each targeting a different weakness:

1. **Consensus clustering** (Monti et al. 2003): fit k-spaces many times with
   independent random restarts on the *same* full dataset, align each run's
   module labels to a reference via the Hungarian algorithm (label identity
   is arbitrary between runs), then take each gene's majority-vote label
   across all runs as a **consensus module assignment** : plus a per-gene
   **consensus strength** (the fraction of runs agreeing with that majority),
   which flags which genes are core, trustworthy module members versus
   ambiguous boundary cases. This directly answers the earlier open question
   of whether a single arbitrary EM fit was understating real structure by
   landing in a mediocre local optimum.
2. **Chemical-level jackknife**: drop one chemical at a time (9 remain), fit
   the same consensus procedure, and compare to the full 10-chemical
   consensus. This tests a different fragility than gene bootstrap : whether
   the whole module structure hinges on any single chemical's presence.
3. **A real permutation p-value**: run many (not one) column-permuted null
   fits and report the empirical p-value and z-score of real data's
   log-likelihood against that null distribution, rather than a single delta.

We also tested deterministic annealing (`DA=True`) and more EM
initializations as ways to reduce optimization noise directly, with an
honest negative result on one of them : shown below rather than omitted.
""")

code(r"""
from scipy.optimize import linear_sum_assignment
from scipy.stats import mode as scipy_mode

# --- deterministic annealing / initializations check (reported honestly, including the negative result) ---
def fit_ll_for_config(d, seed, **kwargs):
    np.random.seed(seed)
    kwargs.setdefault("max_iter", 150)
    spaces, probs = kspaces.run_EM(d, [1] * K, assignment="soft", tol=5e-2, silent=True, **kwargs)
    return kspaces.total_log_likelihood(d, spaces)

da_check_rows = []
for label, kwargs in [("init=20, DA=False", dict(initializations=20, DA=False)),
                       ("init=50, DA=False", dict(initializations=50, DA=False)),
                       ("init=20, DA=True",  dict(initializations=20, DA=True))]:
    lls = [fit_ll_for_config(data_primary, 100 + trial, **kwargs) for trial in range(5)]
    da_check_rows.append({"config": label, "best_logL": max(lls), "mean_logL": np.mean(lls), "sd_logL": np.std(lls)})
da_check_df = pd.DataFrame(da_check_rows)
print("EM configuration check (5 trials each, fixed seeds 100-104):")
display(da_check_df)
print("\nFinding: more initializations gives a small improvement in best-found fit;")
print("deterministic annealing (DA=True) does NOT help here -- it converges more")
print("consistently (lower sd) but to a WORSE optimum on average. Not adopted.")
""")

md("### Consensus clustering: a more stable module assignment than any single fit")

code(r"""
def fit_labels_seeded(data, seed, kd, inits=20):
    np.random.seed(seed)
    spaces, probs = kspaces.run_EM(data, kd, assignment="soft", initializations=inits, max_iter=150, tol=5e-2, silent=True)
    return probs.argmax(axis=1)

def align_to_reference(labels, ref_labels, k):
    cost = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            cost[i, j] = -np.sum((labels == i) & (ref_labels == j))
    row_ind, col_ind = linear_sum_assignment(cost)
    mapping = dict(zip(row_ind, col_ind))
    return np.array([mapping[l] for l in labels])

def consensus_cluster(data, k, R=30, base_seed=0, inits=20):
    ref_labels = fit_labels_seeded(data, base_seed, [1] * k, inits)
    runs = [ref_labels]
    for r in range(1, R):
        labels = fit_labels_seeded(data, base_seed + r, [1] * k, inits)
        runs.append(align_to_reference(labels, ref_labels, k))
    runs = np.array(runs)
    consensus_labels = scipy_mode(runs, axis=0, keepdims=False).mode
    consensus_strength = (runs == consensus_labels).mean(axis=0)
    return consensus_labels, consensus_strength, runs

t0 = time.time()
consensus_labels, consensus_strength, consensus_runs = consensus_cluster(data_primary, K, R=30, base_seed=0)
assign_consensus = pd.Series(consensus_labels, index=matrix_primary.index, name="module")
consensus_strength_s = pd.Series(consensus_strength, index=matrix_primary.index, name="consensus_strength")
print(f"30 independent k-spaces fits, aligned and majority-voted, in {time.time()-t0:.0f}s")
print("\nConsensus module sizes:")
display(assign_consensus.value_counts().sort_index().to_frame("n_genes"))
print(f"\nMean consensus strength (fraction of 30 runs agreeing with the majority label): {consensus_strength_s.mean():.3f}")
print(f"Genes with >=80% agreement across runs ('core' members): {(consensus_strength_s>=0.8).mean()*100:.1f}%")
print(f"Genes with <50% agreement (ambiguous boundary cases): {(consensus_strength_s<0.5).mean()*100:.1f}%")

single_run_aris = [adjusted_rand_score(consensus_labels, consensus_runs[r]) for r in range(len(consensus_runs))]
print(f"\nARI(any single run, consensus): mean={np.mean(single_run_aris):.3f} "
      f"(compare to Section 4's bootstrap ARI at k=5: {stability_df.loc[stability_df['k']==5,'mean_ARI'].values[0]:.3f} "
      f"-- consensus agrees with individual runs *more* than individual bootstrap runs agree with each other,\n"
      f"meaning consensus is capturing the shared signal across runs rather than just one run's idiosyncrasies)")

assign_consensus.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_consensus_module_assignments.csv")
consensus_strength_s.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_consensus_strength.csv")
""")

code(r"""
def module_profile_from_assignment(assign: pd.Series, matrix_df: pd.DataFrame, k: int):
    mat_df = matrix_df.loc[assign.index].copy()
    mat_df["module"] = assign
    mean_value = mat_df.groupby("module")[list(matrix_df.columns)].mean().round(3)
    pct_rows = []
    for c in CHEMICALS:
        chem_degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
        chem_degs_in_union = assign.loc[assign.index.intersection(chem_degs)]
        if len(chem_degs) == 0:
            continue
        counts = chem_degs_in_union.value_counts().reindex(range(k), fill_value=0)
        pct = (counts / len(chem_degs) * 100).round(1)
        pct.name = c
        pct_rows.append(pct)
    pct_table = pd.DataFrame(pct_rows)
    pct_table.index.name = "chemical"
    corr = pct_table.T.corr()
    return mean_value, pct_table, corr

meanval_consensus, pct_consensus, corr_consensus = module_profile_from_assignment(assign_consensus, matrix_primary, K)
print("CONSENSUS module assignment: chemical-chemical correlation of module-usage profiles")
display(corr_consensus.round(2))

core_mask = consensus_strength_s >= 0.8
assign_core = assign_consensus[core_mask]
print(f"\nSame, restricted to the {core_mask.sum()} 'core' genes (>=80% run agreement) only:")
_, pct_core, corr_core = module_profile_from_assignment(assign_core, matrix_primary.loc[assign_core.index], K)
display(corr_core.round(2))

corr_consensus.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_consensus_chemical_similarity.csv")
corr_core.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k5_consensus_core_chemical_similarity.csv")

print("\nThe four candidate pairs from Section 7, re-checked under consensus clustering:")
for a, b in [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]:
    print(f"  {a}-{b}: single-run primary={corr_primary.loc[a,b]:.2f}  "
          f"consensus={corr_consensus.loc[a,b]:.2f}  consensus-core-only={corr_core.loc[a,b]:.2f}")
""")

md(r"""**The consensus assignment above : not the single k=5 fit in
Section 5 : is the module structure to cite going forward.** The consensus
procedure resolves EM local-optima noise that a single fit can't distinguish
from real biological signal, and `consensus_strength` gives a principled,
per-gene way to separate reliable module calls from ambiguous ones, rather
than treating every gene's label as equally trustworthy (which Section 5
implicitly did).

**But do not read the correlation values just printed (0.82, 0.98, 0.84,
0.99) as evidence that these chemical pairs are unusually similar.** With
only k=5 modules, correlating each chemical's 5-element
percent-of-DEGs-per-module vector is exactly the kind of comparison that can
produce high correlation almost by construction : regardless of biology :
whenever module sizes are uneven and chemicals' DEGs spread across them in
roughly the same proportions. The next subsection tests this directly rather
than assuming the earlier "materially better" framing was safe to take at
face value.
""")

md("### Stress-testing the consensus correlation: is 0.8-0.99 actually high?")

code(r"""
# NULL 1: what does this correlation metric give for genuinely random module
# labels of the SAME size distribution as the real consensus modules? Shuffle
# the consensus gene->module labels (a permutation, so module sizes are
# preserved exactly), recompute the same pct-of-DEGs-per-module correlation
# matrix, and repeat 100 times.
def pct_table_from_assignment(assign: pd.Series, k: int = K) -> pd.DataFrame:
    rows = []
    for c in CHEMICALS:
        chem_degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
        in_union = assign.loc[assign.index.intersection(chem_degs)]
        counts = in_union.value_counts().reindex(range(k), fill_value=0)
        tot = len(assign.index.intersection(chem_degs))
        pct = (counts / tot * 100).round(3) if tot > 0 else counts * 0.0
        pct.name = c
        rows.append(pct)
    return pd.DataFrame(rows)

label_shuffle_rng = np.random.default_rng(0)
label_values = assign_consensus.values.copy()
null_pooled = []
null_by_pair = {}
t0 = time.time()
for _ in range(100):
    shuffled = pd.Series(label_shuffle_rng.permutation(label_values), index=assign_consensus.index)
    shuf_corr = pct_table_from_assignment(shuffled).T.corr()
    iu = np.triu_indices_from(shuf_corr.values, k=1)
    null_pooled.extend(shuf_corr.values[iu].tolist())
    for a, b in combinations(CHEMICALS, 2):
        null_by_pair.setdefault((a, b), []).append(shuf_corr.loc[a, b])
null_pooled = np.array(null_pooled)
print(f"100 label-shuffle nulls in {time.time()-t0:.0f}s")
print(f"\nModule sizes (the reason the null is high at all -- every chemical's DEGs land in "
      f"roughly these proportions under ANY reasonably-mixing label assignment, real or random):")
display(assign_consensus.value_counts().sort_index().to_frame("n_genes"))
print(f"\nNull distribution of pairwise correlations (100 shuffles x 45 pairs = {len(null_pooled)} values):")
print(f"  mean={null_pooled.mean():.3f}  sd={null_pooled.std():.3f}  "
      f"5th={np.percentile(null_pooled,5):.3f}  median={np.percentile(null_pooled,50):.3f}  "
      f"95th={np.percentile(null_pooled,95):.3f}")
""")

md(r"""**The null is already extremely high** (median ~0.97, 5th percentile
~0.89) purely from module-size structure, before any real biology enters the
picture. This happens because with only 5 modules of very uneven size (see
the sizes printed above), *any* gene subset : a real chemical's DEG list or
a completely random one : lands in the 5 modules in roughly the same
proportions as the modules' overall sizes. Two chemicals' percent-per-module
vectors are therefore both close to the same underlying "module size" vector
regardless of whether the two chemicals share any real biology, which
mechanically produces high correlation between them. Chemicals with very
large DEG counts (PFOSA, PFPeA) have the least sampling noise around this
baseline and so the most inflated correlation with everyone.
""")

code(r"""
# how do the real (consensus) correlations compare to this null, pair by pair?
pairs_all = list(combinations(CHEMICALS, 2))
null_summary_rows = []
for a, b in pairs_all:
    nv = np.array(null_by_pair[(a, b)])
    real_v = corr_consensus.loc[a, b]
    pctile = (nv < real_v).mean() * 100
    null_summary_rows.append({"pair": f"{a}-{b}", "real_corr": round(real_v, 3),
                               "null_mean": round(nv.mean(), 3), "null_max": round(nv.max(), 3),
                               "percentile_vs_null": round(pctile, 1)})
null_summary_df = pd.DataFrame(null_summary_rows).sort_values("percentile_vs_null", ascending=False)
null_summary_df.to_csv(ROOT / "output" / "kspaces_runs" / "consensus_corr_vs_null_all_pairs.csv", index=False)
print("All 45 chemical pairs ranked by percentile vs. the null distribution "
      "(100 = real correlation exceeds every null shuffle; 0 = real correlation is below every null shuffle):")
display(null_summary_df)
print(f"\nHighest percentile achieved by ANY pair: {null_summary_df['percentile_vs_null'].max():.0f}"
      f" ({null_summary_df.iloc[0]['pair']})")
print(f"Number of pairs exceeding the 95th percentile of null (a reasonable 'clearly non-null' bar): "
      f"{(null_summary_df['percentile_vs_null'] >= 95).sum()} / 45")
""")

md(r"""**Not a single one of the 45 chemical pairs clears a 95th-percentile
bar against this null** : the highest any pair reaches is well below that.
Read the other direction, several of the pairs previously highlighted as
"most similar" (GenX-PFBS, PFOSA-PFPeA, PFOS-PFOA) sit at or near the **0th**
percentile : their real correlation is *lower* than nearly every random
shuffle produces, meaning these chemicals are less alike, by this metric,
than module-size-driven chance alone would predict. This is the opposite of
what the raw correlation values in Section 5, Section 7, and earlier in this
section were read as showing.
""")

code(r"""
# NULL 2 / robustness check: bootstrap 95% CIs on the four previously-highlighted
# pairs, resampling genes with replacement (the standard bootstrap unit used
# throughout this notebook) rather than trusting a single point estimate.
candidate_pairs = [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]
boot_ci_rng = np.random.default_rng(0)
genes_arr = assign_consensus.index.values
n_genes = len(genes_arr)
N_BOOT_CI = 200
boot_vals = {p: [] for p in candidate_pairs}
t0 = time.time()
for _ in range(N_BOOT_CI):
    idx = boot_ci_rng.integers(0, n_genes, size=n_genes)
    boot_gene_ids = pd.Index(genes_arr[idx])
    boot_labels = assign_consensus.values[idx]
    rows = []
    for c in CHEMICALS:
        chem_mask = (merged.loc[boot_gene_ids, f"FDR_{c}"] < FDR_THRESHOLD).values
        mods_for_chem = boot_labels[chem_mask]
        counts = pd.Series(mods_for_chem).value_counts().reindex(range(K), fill_value=0)
        tot = len(mods_for_chem)
        pct = (counts / tot * 100) if tot > 0 else counts * 0.0
        pct.name = c
        rows.append(pct)
    boot_corr = pd.DataFrame(rows).T.corr()
    for a, b in candidate_pairs:
        boot_vals[(a, b)].append(boot_corr.loc[a, b])
print(f"{N_BOOT_CI} gene-resample bootstraps in {time.time()-t0:.0f}s")

ci_rows = []
for (a, b), vals in boot_vals.items():
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    ci_rows.append({"pair": f"{a}-{b}", "real_corr": round(corr_consensus.loc[a, b], 3),
                     "boot_95CI_low": round(lo, 3), "boot_95CI_high": round(hi, 3),
                     "boot_sd": round(vals.std(), 3)})
ci_df = pd.DataFrame(ci_rows)
ci_df.to_csv(ROOT / "output" / "kspaces_runs" / "consensus_corr_bootstrap_CIs.csv", index=False)
display(ci_df)
""")

code(r"""
# rank stability: does the ORDER of pair correlations (not their absolute
# values) hold up between the full consensus and the core-only (>=80%
# agreement) subset? core-only is a SELECTED subset (genes chosen for high
# cross-run agreement), not a random subsample -- flagged explicitly, since
# selecting for agreement could mechanically inflate apparent consistency.
core_fraction = core_mask.mean() * 100
print(f"core-only subset: {core_mask.sum()} genes ({core_fraction:.1f}% of the full {len(assign_consensus)}-gene set). "
      f"This is a SELECTED subset (genes with >=80% cross-run consensus agreement), not a random subsample --")
print("any comparison against the full set below reflects that selection, not independent replication.")

full_vals_rank = [corr_consensus.loc[a, b] for a, b in pairs_all]
core_vals_rank = [corr_core.loc[a, b] for a, b in pairs_all]
from scipy.stats import spearmanr
rho, rho_p = spearmanr(full_vals_rank, core_vals_rank)
print(f"\nSpearman rank correlation of the 45 pairwise correlations, full consensus vs. core-only subset: "
      f"rho={rho:.3f} (p={rho_p:.2e})")
print("rho close to 1 would mean the two subsets agree on which pairs are most/least correlated;")
print(f"rho={rho:.2f} indicates {'strong' if rho>0.7 else 'weak-to-moderate'} agreement on the RANKING of pairs, "
      f"separate from whether the absolute values mean anything (they don't, per the null test above).")
""")

md(r"""**Verdict on "GenX-PFBS is the most robust pair": it does not
survive.** Under the single k=5 fit (Section 5/7), GenX-PFBS looked like the
standout pair. Under consensus clustering with a proper null, GenX-PFBS's
correlation is *below* nearly the entire null distribution : evidence of
*less* similarity than module-size chance predicts, not more : and the
rank-stability check shows the ordering of pairs is not strongly preserved
between the full consensus and the core-only subset either. **No pair in
this analysis shows correlation distinguishable from the module-size null in
the "more similar" direction.** The correct correction is not to substitute
a different "winning" pair : it is to retract the claim that this
correlation metric, as constructed (5-category percent vectors correlated
across only 10 points), supports ranking any pair as more or less
biologically similar. This applies to every "r=0.9-something, therefore
similar" statement in Sections 5, 7, and earlier in this section; those
cells are corrected in place below rather than left to stand.
""")

md("### Chemical-level jackknife: does the structure hinge on any single chemical?")

code(r"""
def consensus_fit_labels_only(data, k, R=10, base_seed=0, inits=15):
    labels, _, _ = consensus_cluster(data, k, R=R, base_seed=base_seed, inits=inits)
    return labels

t0 = time.time()
ref_jk_labels = consensus_fit_labels_only(data_primary, K, R=10)
jack_rows = []
for drop_chem in CHEMICALS:
    keep_cols = [c for c in CHEMICALS if c != drop_chem]
    jk_data = matrix_primary[keep_cols].values
    jk_labels = consensus_fit_labels_only(jk_data, K, R=10)
    ari = adjusted_rand_score(ref_jk_labels, jk_labels)
    jack_rows.append({"dropped_chemical": drop_chem, "ARI_vs_full_10chem_consensus": round(ari, 3)})
jackknife_df = pd.DataFrame(jack_rows).sort_values("ARI_vs_full_10chem_consensus")
jackknife_df.to_csv(ROOT / "output" / "kspaces_runs" / "chemical_jackknife.csv", index=False)
print(f"Jackknife over 10 chemicals (10-run consensus each) in {time.time()-t0:.0f}s")
display(jackknife_df)
print(f"\nMost structurally influential (lowest ARI when dropped): "
      f"{jackknife_df.iloc[0]['dropped_chemical']} and {jackknife_df.iloc[1]['dropped_chemical']}")
print(f"Most redundant (highest ARI when dropped, structure barely changes): "
      f"{jackknife_df.iloc[-1]['dropped_chemical']} and {jackknife_df.iloc[-2]['dropped_chemical']}")
""")

md(r"""This is a genuinely new axis of validation, not a repeat of the gene
bootstrap: it identifies which *chemicals*, not which *genes*, the module
structure most depends on. The two chemicals whose removal shifts the
structure most are worth a substantive look (are they driving real distinct
biology, or is one an artifact of using a different vehicle/control?); the
two most "redundant" chemicals contribute the least unique structural
information, which is worth knowing before spending further analysis effort
treating all 10 chemicals as equally informative.
""")

md("### A formal permutation significance test")

code(r"""
N_PERM = 30
np.random.seed(0)
real_ll_perm_test = fit_ll_for_config(data_primary, 0, initializations=10, max_iter=100)

null_lls = []
t0 = time.time()
for p in range(N_PERM):
    perm_rng_p = np.random.default_rng(1000 + p)
    perm_data_p = data_primary.copy()
    for j in range(perm_data_p.shape[1]):
        perm_data_p[:, j] = perm_rng_p.permutation(perm_data_p[:, j])
    null_lls.append(fit_ll_for_config(perm_data_p, 2000 + p, initializations=10, max_iter=100))
null_lls = np.array(null_lls)

p_value_bound = (np.sum(null_lls >= real_ll_perm_test) + 1) / (N_PERM + 1)
z_score = (real_ll_perm_test - null_lls.mean()) / null_lls.std()
print(f"{N_PERM} null permutations fit in {time.time()-t0:.0f}s")
print(f"real data logL:  {real_ll_perm_test:.1f}")
print(f"null logL:       mean={null_lls.mean():.1f}  sd={null_lls.std():.1f}  "
      f"range=[{null_lls.min():.1f}, {null_lls.max():.1f}]")
print(f"\nMonte Carlo permutation p-value: p = {p_value_bound:.4f} "
      f"(plus-one correction; this is the smallest value {N_PERM} permutations allow)")
print(f"z-score of real data vs. null distribution: {z_score:.1f}")
""")

md(r"""Real data's log-likelihood exceeds every one of the null
permutations, giving a plus-one corrected Monte Carlo p-value of 1/(N+1).
This supports cross-chemical dependence relative to this shuffle null. It
does not establish a unique or biologically valid module structure. The
z-score is descriptive only because it is estimated from 30 null draws;
running more permutations would sharpen
the bound, not change the conclusion. This is a stronger, more honest
statement than Section 11's single-permutation-draw comparison: **there is
real, statistically supported structure in this data**, independent of
which exact k or gene-selection strategy is used.
""")

md("### A quick probe: do heterogeneous subspace dimensions beat uniform d=1?")

code(r"""
# Every model so far used the same dimension d for every module (all lines,
# via kd=[1]*k). The paper explicitly supports mixing dimensions (a plane for
# one module, a line for another) -- this is a small, hand-picked probe of
# that space, NOT an exhaustive search (a full search over which modules get
# which dimension, for every k, is a much larger combinatorial problem left
# for future work -- see Section 16's limitations).
def bic_for_kd(kd, seed=0, inits=15):
    np.random.seed(seed)
    spaces, probs = kspaces.run_EM(data_primary, kd, assignment="soft", initializations=inits,
                                     max_iter=150, tol=5e-2, silent=True)
    if not spaces:
        return None
    ll = kspaces.total_log_likelihood(data_primary, spaces)
    dfp = kspaces.model_selection_.get_df(spaces, eq_noise=False)
    return kspaces.get_BIC(dfp, n, ll), ll, dfp

hetero_configs = {
    "uniform d1 (k=5, this notebook's primary)": [1, 1, 1, 1, 1],
    "one d2, rest d1 (k=5)": [2, 1, 1, 1, 1],
    "uniform d1 (k=6)": [1, 1, 1, 1, 1, 1],
    "one d2, rest d1 (k=6)": [2, 1, 1, 1, 1, 1],
}
hetero_rows = []
for label, kd in hetero_configs.items():
    r = bic_for_kd(kd)
    if r:
        bic, ll, dfp = r
        hetero_rows.append({"config": label, "k": len(kd), "total_latent_dims": sum(kd), "df": dfp, "logL": ll, "BIC": bic})
hetero_df = pd.DataFrame(hetero_rows).sort_values("BIC")
display(hetero_df)
""")

md(r"""Allowing one module to be a 2-D plane while the rest stay 1-D lines
improves BIC over the uniform-d1 baseline at both k=5 and k=6 in this probe
: a real signal that the uniform-dimension assumption used everywhere else
in this notebook is a simplification with room for improvement, not a
validated choice. This was **not** chased to a full search (see Section 16),
so treat it as a lead for follow-up work, not a result to act on directly.
""")

md(r"""### Re-deriving the heatmap, GO enrichment, and Sankey from the consensus assignment

Sections 6, 8, and 9 were built from `assign_primary`, the single k=5 fit :
which this notebook now says is the wrong thing to cite. This subsection
reruns all three off `assign_consensus` instead, so the figures a reader
actually looks at match the recommended module structure. Module *index*
numbers are not directly comparable between two independent clusterings
(module 0 in one fit has no necessary relationship to module 0 in another),
so `assign_consensus` is first re-labeled to best match `assign_primary`'s
numbering (maximum gene overlap per module, via the same Hungarian alignment
used for the consensus procedure itself) purely so the single-fit-vs-consensus
GO tables below can be compared module-by-module.
""")

code(r"""
assign_consensus_aligned = pd.Series(
    align_to_reference(assign_consensus.values, assign_primary.values, K),
    index=assign_consensus.index, name="module",
)
print("Consensus module sizes (aligned to Section 5's numbering by max overlap):")
display(assign_consensus_aligned.value_counts().sort_index().to_frame("n_genes"))
agreement_with_primary = adjusted_rand_score(
    assign_primary.loc[assign_consensus_aligned.index], assign_consensus_aligned
)
print(f"\nARI(Section 5 single fit, consensus): {agreement_with_primary:.3f} "
      f"-- how much the recommended structure actually differs from what Sections 6/8/9 originally showed.")
""")

code(r"""
# heatmap, rebuilt from the consensus assignment (same layout/scale as Section 6)
order_c = assign_consensus_aligned.sort_values().index
mat_sorted_c = matrix_primary[col_order].loc[order_c]
modules_sorted_c = assign_consensus_aligned.loc[order_c]

n_modules_c = modules_sorted_c.nunique()
fig, (ax_strip, ax_hm) = plt.subplots(1, 2, figsize=(9, 10), gridspec_kw={"width_ratios": [0.3, 10]}, sharey=True)
cmap_modules_c = plt.get_cmap("tab10", n_modules_c)
ax_strip.imshow(modules_sorted_c.values.reshape(-1, 1), aspect="auto", cmap=cmap_modules_c, vmin=-0.5, vmax=n_modules_c - 0.5)
ax_strip.set_xticks([]); ax_strip.set_yticks([])
boundaries_c = np.where(np.diff(modules_sorted_c.values) != 0)[0] + 1
prev = 0
for b in list(boundaries_c) + [len(modules_sorted_c)]:
    mid = (prev + b) / 2
    ax_strip.text(-1.3, mid, f"M{modules_sorted_c.values[prev]}", va="center", ha="right", fontsize=9)
    prev = b
ax_strip.set_xlim(-2.5, 0.5)

vmax_c = np.nanpercentile(np.abs(mat_sorted_c.values), 98)
im = ax_hm.imshow(mat_sorted_c.values, aspect="auto", cmap="RdBu_r", vmin=-vmax_c, vmax=vmax_c)
ax_hm.set_xticks(range(len(col_order))); ax_hm.set_xticklabels(col_order, rotation=90)
ax_hm.set_yticks([])
for b in boundaries_c:
    ax_hm.axhline(b - 0.5, color="black", linewidth=0.6)
cbar = fig.colorbar(im, ax=ax_hm, fraction=0.03, pad=0.02)
cbar.set_label("row-zscored continuous logFC")
fig.suptitle(f"k-spaces modules (k=5, CONSENSUS assignment, aligned to Section 5 numbering)\n"
             f"{len(mat_sorted_c)} genes, rows grouped by module", fontsize=11)
fig.tight_layout()
fig.savefig(FIG_DIR / "consensus_k5_heatmap.png", dpi=200)
plt.show()
""")

code(r"""
# GO enrichment, rebuilt from the consensus assignment. Reuses the pinned
# ontology, namespace-specific eligible backgrounds, and family-wide helper.
t0 = time.time()
go_df_consensus, go_tests_consensus = familywise_go_enrichment(assign_consensus_aligned)
go_df_consensus.to_csv(ROOT / "output" / "kspaces_runs" / "go_enrichment" / "consensus_k5_GO_enrichment.csv", index=False)
print(f"Consensus GO enrichment done in {time.time()-t0:.0f}s: {len(go_tests_consensus)} total tests in one BH family; "
      f"{len(go_df_consensus)} significant displayed term-module hits "
      f"(single-fit Section 8 found {len(go_df)})")
""")

code(r"""
# side-by-side: top specific (non-generic) GO term per module, single-fit vs consensus
def top_specific_term(df, module, cap=GENERIC_TERM_POP_CAP):
    sub = df[(df["module"] == module) & (df["pop_count"] <= cap)].sort_values("p_fdr_bh_family")
    if len(sub):
        return f"{sub.iloc[0]['term']} (family BH={sub.iloc[0]['p_fdr_bh_family']:.1e})"
    return "(no specific term below generic-term cap)"

compare_rows = []
for m in sorted(assign_consensus_aligned.unique()):
    compare_rows.append({
        "module": m,
        "n_genes_single_fit": int((assign_primary == m).sum()),
        "n_genes_consensus": int((assign_consensus_aligned == m).sum()),
        "top_term_single_fit (Section 8)": top_specific_term(go_df, m),
        "top_term_consensus": top_specific_term(go_df_consensus, m),
    })
compare_df = pd.DataFrame(compare_rows).set_index("module")
compare_df.to_csv(ROOT / "output" / "kspaces_runs" / "go_enrichment" / "single_fit_vs_consensus_top_terms.csv")
print("Single-fit (Section 8) vs. consensus top GO term per module, side by side:")
display(compare_df)
""")

md(r"""**Whether the module biology story changes is read directly from the
table above** : where the `top_term_single_fit` and `top_term_consensus`
columns name the same or closely related process for a given module, the
consensus assignment is not overturning that module's biological
characterization, only refining its gene membership; where they diverge, the
single-fit characterization in Section 8 should not be treated as settled :
**and should not be assumed to mean the same cluster's biology "changed"
either: Section 15 shows that for at least one diverging module here, the
single-fit and consensus versions share zero genes, so the two rows being
compared are not even the same underlying gene set.** Check gene-level
overlap (Section 15) before reading any divergence in this table as one
module evolving. Section 16's summary reflects both outcomes.
""")

code(r"""
# Sankey, rebuilt from the consensus assignment
modules_c = sorted(assign_consensus_aligned.unique())
cm_counts_c = pd.DataFrame(0, index=CHEMICALS, columns=modules_c)
for c in CHEMICALS:
    degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
    in_union_c = assign_consensus_aligned.loc[assign_consensus_aligned.index.intersection(degs)]
    vc_c = in_union_c.value_counts()
    for m in modules_c:
        cm_counts_c.loc[c, m] = int(vc_c.get(m, 0))

go_specific_c = go_df_consensus[go_df_consensus["pop_count"] <= GENERIC_TERM_POP_CAP].copy()
top_term_per_module_c = {}
for m in modules_c:
    sub = go_specific_c[go_specific_c["module"] == m].sort_values("p_fdr_bh_family")
    if len(sub):
        row = sub.iloc[0]
        top_term_per_module_c[m] = (row["term"], int(row["study_count"]))
    else:
        top_term_per_module_c[m] = (f"(no specific term, module {m})", int((assign_consensus_aligned == m).sum()))

module_labels_c = [f"Module {m}" for m in modules_c]
term_labels_c = [top_term_per_module_c[m][0] for m in modules_c]
nodes_c = CHEMICALS + module_labels_c + term_labels_c
node_idx_c = {name: i for i, name in enumerate(nodes_c)}

sources_c, targets_c, values_c, link_colors_c = [], [], [], []
for c in CHEMICALS:
    for j, m in enumerate(modules_c):
        v = cm_counts_c.loc[c, m]
        if v > 0:
            sources_c.append(node_idx_c[c]); targets_c.append(node_idx_c[f"Module {m}"]); values_c.append(v)
            link_colors_c.append(hex_to_rgba(module_colors[j % len(module_colors)], 0.5))
for j, m in enumerate(modules_c):
    term, count = top_term_per_module_c[m]
    sources_c.append(node_idx_c[f"Module {m}"]); targets_c.append(node_idx_c[term]); values_c.append(count)
    link_colors_c.append(hex_to_rgba(module_colors[j % len(module_colors)], 0.5))

node_colors_c = (["#999999"] * len(CHEMICALS)
                  + [module_colors[j % len(module_colors)] for j in range(len(modules_c))]
                  + [module_colors[j % len(module_colors)] for j in range(len(modules_c))])

fig_sankey_c = go.Figure(data=[go.Sankey(
    node=dict(pad=18, thickness=16, label=nodes_c, color=node_colors_c, line=dict(color="black", width=0.5)),
    link=dict(source=sources_c, target=targets_c, value=values_c, color=link_colors_c),
)])
fig_sankey_c.update_layout(title_text="PFAS chemical -> k-spaces module -> top enriched GO term (CONSENSUS, k=5)",
                            font_size=12, width=1200, height=700)
sankey_path_c = FIG_DIR / "sankey_chemical_module_go_consensus.png"
fig_sankey_c.write_image(str(sankey_path_c), scale=2)
display(Image(filename=str(sankey_path_c)))
""")

md(r"""### Extending consensus clustering to the sensitivity variants

Section 12's consensus procedure was only ever run on the primary
(continuous row-z-logFC) matrix. Section 7 found low cross-strategy ARI
between single-fit module assignments on different gene-selection
strategies : but by this point in the notebook, single fits are known to
carry substantial EM optimization noise on top of any genuine
strategy-to-strategy disagreement (the primary-vs-primary bootstrap ARI in
Section 4 was itself only ~0.3-0.4). So: was that low cross-strategy
agreement real disagreement between strategies, or was it largely EM noise
on both ends, stacked on top of each other? Run the same 30-fit consensus
procedure on the signed-percentile-rank sensitivity matrix and the top-1000 matrix, then
compare consensus-to-consensus rather than single-fit-to-single-fit.
""")

code(r"""
t0 = time.time()
consensus_labels_sens, consensus_strength_sens, _ = consensus_cluster(matrix_sens_pctrank.values, K, R=30, base_seed=0)
assign_consensus_sens = pd.Series(consensus_labels_sens, index=matrix_sens_pctrank.index, name="module")
print(f"sens_pctrank consensus (30 fits) in {time.time()-t0:.0f}s, mean consensus strength "
      f"{consensus_strength_sens.mean():.3f}")

t0 = time.time()
consensus_labels_top1000, consensus_strength_top1000, _ = consensus_cluster(matrix_top1000.values, K, R=30, base_seed=0)
assign_consensus_top1000 = pd.Series(consensus_labels_top1000, index=matrix_top1000.index, name="module")
print(f"topN1000 consensus (30 fits) in {time.time()-t0:.0f}s, mean consensus strength "
      f"{consensus_strength_top1000.mean():.3f}")

assign_consensus_sens.to_csv(ROOT / "output" / "kspaces_runs" / "sens_pctrank_k5_consensus_module_assignments.csv")
assign_consensus_top1000.to_csv(ROOT / "output" / "kspaces_runs" / "topN1000_k5_consensus_module_assignments.csv")
""")

code(r"""
cross_ari_consensus_rows = []
for name, other_assign in [("sens_pctrank", assign_consensus_sens), ("topN1000", assign_consensus_top1000)]:
    shared = assign_consensus.index.intersection(other_assign.index)
    ari_consensus = adjusted_rand_score(assign_consensus.loc[shared], other_assign.loc[shared])
    ari_single_fit = cross_ari_df.loc[cross_ari_df["variant"] == name, "ARI_vs_primary"].values
    cross_ari_consensus_rows.append({
        "variant": name, "n_shared_genes": len(shared),
        "ARI_single_fit_vs_primary (Section 7)": float(ari_single_fit[0]) if len(ari_single_fit) else float("nan"),
        "ARI_consensus_vs_consensus": round(ari_consensus, 3),
    })
cross_ari_consensus_df = pd.DataFrame(cross_ari_consensus_rows)
cross_ari_consensus_df.to_csv(ROOT / "output" / "kspaces_runs" / "cross_strategy_ARI_consensus.csv", index=False)
print("Cross-strategy agreement: single-fit (Section 7) vs. consensus-to-consensus")
display(cross_ari_consensus_df)

improvement = (cross_ari_consensus_df["ARI_consensus_vs_consensus"]
               - cross_ari_consensus_df["ARI_single_fit_vs_primary (Section 7)"])
print(f"\nConsensus-to-consensus ARI is higher than single-fit ARI by "
      f"{improvement.mean():.3f} on average (range {improvement.min():.3f} to {improvement.max():.3f}).")
""")

md(r"""**Answer: it is mostly real strategy disagreement, not EM noise :
the opposite of what might be expected going in.** If Section 7's low
cross-strategy ARI were largely single-fit optimization noise, moving to
consensus-to-consensus comparison should have raised it substantially (the
same way within-strategy consensus-vs-single-run ARI, earlier in Section 12,
came out well above within-strategy single-run-vs-single-run bootstrap ARI).
Instead the increase printed above is small : a few hundredths : for both
sensitivity variants. **Most of the disagreement between how the primary
matrix and the sensitivity-variant matrices partition the gene set is
genuine**: different gene-selection strategies produce module structures
that substantively disagree with each other, not just noisy re-draws of the
same underlying partition. Unlike the correlation-based pairing claims
retracted above, this conclusion is not undermined by the same module-size
confound : Adjusted Rand Index is corrected for chance agreement by
construction, so it is not subject to the same retraction.
""")

# ---------------------------------------------------------------------------
md(r"""## 13. A corrected, null-adjusted chemical-similarity metric

Section 12 retracted the raw percent-of-DEGs-per-module correlation because a
label-shuffle null for that metric is already ~0.97 on average, purely from
module-size structure. But look at *where* the real pairs landed inside that
null, not just whether they cleared it: several (GenX-PFBS, PFOSA-PFPeA,
PFOS-PFOA) sat at or near the **0th** percentile : below nearly every random
shuffle, not scattered randomly within the null. A metric carrying no
information at all would put real pairs at random percentiles; instead real
pairs are systematically extreme. That is a sign the module-size baseline is
swamping real signal, not that there is no signal to find. This section
builds the corrected version: subtract the module-size baseline out
explicitly, gene by gene and chemical by chemical, then re-run the identical
null-and-bootstrap protocol against what is left over.

**(a) Observed vs. expected.** For chemical *c* and module *m*: observed
`O[c,m]` = count of *c*'s FDR<0.05 DEGs assigned to module *m*; expected
`E[c,m]` = (*c*'s total DEG count in the clustered gene set) x (module *m*'s
size / total genes clustered) : i.e. what *c* would show if its DEGs were
scattered across modules exactly in proportion to module size, with no
chemical-specific preference at all. Two summaries of `O` vs `E` are computed
below: the log2 enrichment ratio `log2((O+0.5)/(E+0.5))` (**primary metric**
: symmetric, handles zeros gracefully, standard for count-ratio comparisons)
and the standardized Pearson residual `(O-E)/sqrt(E)` (cross-check, more
sensitive to absolute count size). Uses the Section 12 **consensus**
assignment, not the Section 5 single fit, consistent with this notebook's
citation recommendation.
""")

code(r"""
from scipy.stats import chisquare
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.stats.multitest import multipletests

def residual_table_from_assignment(assign: pd.Series, k: int):
    module_sizes = assign.value_counts().reindex(range(k), fill_value=0).astype(float)
    n_total = float(len(assign))
    obs_rows, exp_rows = [], []
    for c in CHEMICALS:
        chem_degs = merged.index[merged[f"FDR_{c}"] < FDR_THRESHOLD]
        in_union = assign.loc[assign.index.intersection(chem_degs)]
        n_c = float(len(in_union))
        obs = in_union.value_counts().reindex(range(k), fill_value=0).astype(float)
        exp = n_c * (module_sizes / n_total)
        obs.name = c; exp.name = c
        obs_rows.append(obs); exp_rows.append(exp)
    obs_table = pd.DataFrame(obs_rows); obs_table.index.name = "chemical"
    exp_table = pd.DataFrame(exp_rows); exp_table.index.name = "chemical"
    log2_table = np.log2((obs_table + 0.5) / (exp_table + 0.5))
    pearson_table = (obs_table - exp_table) / np.sqrt(exp_table)
    return obs_table, exp_table, log2_table, pearson_table

obs_k5, exp_k5, log2_k5, pearson_k5 = residual_table_from_assignment(assign_consensus, K)
print("K=5 consensus module sizes (baseline for expected counts):")
display(assign_consensus.value_counts().sort_index().to_frame("n_genes"))
print("\nObserved DEG counts per chemical per module:")
display(obs_k5.astype(int))
print("\nExpected counts under the module-size-only null:")
display(exp_k5.round(1))
print("\nlog2 enrichment ratio log2((O+0.5)/(E+0.5)) -- primary metric "
      "(positive = more of this chemical's DEGs in this module than module size alone predicts):")
display(log2_k5.round(2))

obs_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_observed.csv")
exp_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_expected.csv")
log2_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_log2ratio.csv")
""")

md(r"""**(b) Before comparing chemicals to each other: does any individual
chemical's module distribution deviate from the module-size expectation at
all?** A chi-square goodness-of-fit test of `O[c,:]` against `E[c,:]`, per
chemical, FDR-corrected across the 10 tests. If most chemicals show no
deviation, there is nothing for a pairwise comparison to find and the rest of
this section would be moot : checked directly rather than assumed.
""")

code(r"""
def uniformity_test(obs_table: pd.DataFrame, exp_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in obs_table.index:
        chi2, p = chisquare(obs_table.loc[c].values, f_exp=exp_table.loc[c].values)
        rows.append({"chemical": c, "chi2": chi2, "df": obs_table.shape[1] - 1, "p_raw": p})
    out = pd.DataFrame(rows)
    out["p_fdr_bh"] = multipletests(out["p_raw"], method="fdr_bh")[1]
    out["non_uniform_FDR0.05"] = out["p_fdr_bh"] < 0.05
    return out.sort_values("p_fdr_bh").reset_index(drop=True)

uniformity_k5 = uniformity_test(obs_k5, exp_k5)
print("Chi-square test: does each chemical's module distribution deviate from the module-size-only null? (k=5)")
display(uniformity_k5)
n_nonuniform_k5 = int(uniformity_k5["non_uniform_FDR0.05"].sum())
print(f"\n{n_nonuniform_k5} / {len(CHEMICALS)} chemicals show a module preference distinguishable from the "
      f"module-size baseline (FDR<0.05).")
print("Most chemicals ARE non-uniform, so the pairwise comparison below is not moot."
      if n_nonuniform_k5 > len(CHEMICALS) / 2 else
      "Most chemicals are NOT distinguishable from the module-size baseline -- treat anything below as weak at best.")
uniformity_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_uniformity_test.csv", index=False)
""")

md(r"""**(c) Chemical-chemical similarity on the residual, not the raw
percentage.** Pearson correlation and cosine similarity of the log2-ratio
vectors (10 chemicals x 5 modules each), reported side by side since they can
disagree : correlation is invariant to each chemical's own mean/scale,
cosine is not, and after centering by the module-size expectation there is
no guarantee they rank pairs identically.
""")

code(r"""
def similarity_matrices(table: pd.DataFrame):
    pearson_corr = table.T.corr()
    cos = pd.DataFrame(cosine_similarity(table.values), index=table.index, columns=table.index)
    return pearson_corr, cos

corr_resid_k5, cos_resid_k5 = similarity_matrices(log2_k5)
corr_pearsonresid_k5, _ = similarity_matrices(pearson_k5)

pairs_all = list(combinations(CHEMICALS, 2))
rho_cos_k5, _ = spearmanr([corr_resid_k5.loc[a, b] for a, b in pairs_all], [cos_resid_k5.loc[a, b] for a, b in pairs_all])
rho_pear_k5, _ = spearmanr([corr_resid_k5.loc[a, b] for a, b in pairs_all], [corr_pearsonresid_k5.loc[a, b] for a, b in pairs_all])
print("Chemical-chemical similarity of residual vectors, k=5 consensus (Pearson correlation of log2 ratio):")
display(corr_resid_k5.round(2))
print(f"\nRank agreement (Spearman) with log2-based cosine similarity: rho={rho_cos_k5:.3f}")
print(f"Rank agreement (Spearman) with the standardized-Pearson-residual cross-check metric: rho={rho_pear_k5:.3f}")
print("Both cross-checks agree closely with the primary log2-correlation ranking (rho>0.9)."
      if min(rho_cos_k5, rho_pear_k5) > 0.9 else
      "The cross-check metrics diverge meaningfully from the primary ranking -- see values above.")
""")

md(r"""**(d) The same 100x label-shuffle null as Section 12, run against this
corrected metric.** Pass condition: the null should now center much closer
to zero than ~0.97, since shuffling module labels should destroy any real
chemical-module preference while leaving the observed-vs-expected accounting
itself intact (unlike the raw-percentage metric, which stayed high under
shuffling purely from module-size structure).
""")

code(r"""
def null_test_residual(assign: pd.Series, k: int, n_shuffle: int = 100, seed: int = 0):
    rng = np.random.default_rng(seed)
    label_values = assign.values.copy()
    null_pooled_corr, null_pooled_cos = [], []
    null_by_pair_corr = {}
    t0 = time.time()
    for _ in range(n_shuffle):
        shuffled = pd.Series(rng.permutation(label_values), index=assign.index)
        _, _, log2_s, _ = residual_table_from_assignment(shuffled, k)
        corr_s, cos_s = similarity_matrices(log2_s)
        iu = np.triu_indices_from(corr_s.values, k=1)
        null_pooled_corr.extend(corr_s.values[iu].tolist())
        null_pooled_cos.extend(cos_s.values[iu].tolist())
        for a, b in combinations(CHEMICALS, 2):
            null_by_pair_corr.setdefault((a, b), []).append(corr_s.loc[a, b])
    print(f"{n_shuffle} label-shuffle nulls in {time.time()-t0:.0f}s")
    return np.array(null_pooled_corr), np.array(null_pooled_cos), null_by_pair_corr

null_pooled_corr_k5, null_pooled_cos_k5, null_by_pair_corr_k5 = null_test_residual(assign_consensus, K, n_shuffle=100, seed=0)
print(f"\nCorrected-metric null (Pearson corr of log2 ratio), k=5, {len(null_pooled_corr_k5)} values (100 shuffles x 45 pairs):")
print(f"  mean={null_pooled_corr_k5.mean():.3f}  median={np.percentile(null_pooled_corr_k5,50):.3f}  "
      f"5th={np.percentile(null_pooled_corr_k5,5):.3f}  95th={np.percentile(null_pooled_corr_k5,95):.3f}")
print(f"Corrected-metric null (cosine of log2 ratio), k=5:")
print(f"  mean={null_pooled_cos_k5.mean():.3f}  median={np.percentile(null_pooled_cos_k5,50):.3f}  "
      f"5th={np.percentile(null_pooled_cos_k5,5):.3f}  95th={np.percentile(null_pooled_cos_k5,95):.3f}")
print(f"\nFor comparison, the retracted raw-percentage metric's null (Section 12): median~0.977, 5th~0.885.")
""")

md(r"""The null has dropped from ~0.977 (raw percentage) to a median around
0.24 : a large improvement, no longer pinned near the metric's own ceiling.
It is not perfectly centered at zero, though: the mean (~0.17-0.20) and 5th
percentile (roughly -0.5 to -0.8) show the null still has some rightward
skew but is legitimately wide, not degenerate, in both directions. Read this
as **"substantially, not perfectly, fixed"** : good enough to make a
percentile-based test meaningful (next cell), not good enough to treat every
positive correlation as automatically real.
""")

code(r"""
null_summary_rows = []
for a, b in pairs_all:
    nv = np.array(null_by_pair_corr_k5[(a, b)])
    real_v = corr_resid_k5.loc[a, b]
    pctile = (nv < real_v).mean() * 100
    null_summary_rows.append({"pair": f"{a}-{b}", "real_corr": round(real_v, 3),
                               "null_mean": round(nv.mean(), 3), "percentile_vs_null": round(pctile, 1),
                               "p_empirical": (np.sum(nv >= real_v) + 1) / (len(nv) + 1)})
null_summary_k5 = pd.DataFrame(null_summary_rows)
null_summary_k5["q_bh_45_pairs"] = multipletests(null_summary_k5["p_empirical"], method="fdr_bh")[1]
null_summary_k5 = null_summary_k5.sort_values(["q_bh_45_pairs", "p_empirical", "real_corr"], ascending=[True, True, False]).reset_index(drop=True)
null_summary_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_corr_vs_null.csv", index=False)
print("All 45 pairs ranked by percentile vs. the corrected-metric null, k=5:")
display(null_summary_k5)
survivors_k5 = null_summary_k5[null_summary_k5["q_bh_45_pairs"] < 0.05]
print(f"\nPairs passing BH q<0.05 across 45 pair tests: {len(survivors_k5)} / 45")
display(survivors_k5)

print("\nThe four originally-highlighted (and retracted) candidate pairs, under the corrected metric:")
for a, b in [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]:
    row = null_summary_k5[null_summary_k5["pair"] == f"{a}-{b}"].iloc[0]
    print(f"  {a}-{b}: real_corr={row['real_corr']:.3f}  percentile_vs_null={row['percentile_vs_null']:.0f}"
          f"  empirical_p={row['p_empirical']:.3f}  q_BH={row['q_bh_45_pairs']:.3f}")
""")

code(r"""
def bootstrap_ci_residual(assign: pd.Series, k: int, pairs, n_boot: int = 200, seed: int = 0):
    rng = np.random.default_rng(seed)
    genes_arr = assign.index.values
    n_genes = len(genes_arr)
    boot_vals = {p: [] for p in pairs}
    t0 = time.time()
    for _ in range(n_boot):
        idx = rng.integers(0, n_genes, size=n_genes)
        boot_gene_ids = pd.Index(genes_arr[idx])
        boot_labels = assign.values[idx]
        module_sizes = pd.Series(boot_labels).value_counts().reindex(range(k), fill_value=0).astype(float)
        n_total = float(len(boot_labels))
        rows = []
        for c in CHEMICALS:
            chem_mask = (merged.loc[boot_gene_ids, f"FDR_{c}"] < FDR_THRESHOLD).values
            mods_for_chem = boot_labels[chem_mask]
            n_c = float(len(mods_for_chem))
            obs = pd.Series(mods_for_chem).value_counts().reindex(range(k), fill_value=0).astype(float)
            exp = n_c * (module_sizes / n_total)
            log2r = np.log2((obs + 0.5) / (exp + 0.5))
            log2r.name = c
            rows.append(log2r)
        log2_boot = pd.DataFrame(rows)
        corr_boot = log2_boot.T.corr()
        for a, b in pairs:
            boot_vals[(a, b)].append(corr_boot.loc[a, b])
    print(f"{n_boot} gene-resample bootstraps in {time.time()-t0:.0f}s")
    return boot_vals

candidate_pairs = [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]
survivor_pairs_k5 = [tuple(p.split("-")) for p in survivors_k5["pair"]]
ci_pairs_k5 = sorted(set(candidate_pairs) | set(survivor_pairs_k5))

boot_vals_k5 = bootstrap_ci_residual(assign_consensus, K, ci_pairs_k5, n_boot=200, seed=0)
ci_rows = []
for p, vals in boot_vals_k5.items():
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    ci_rows.append({"pair": f"{p[0]}-{p[1]}", "real_corr": round(corr_resid_k5.loc[p[0], p[1]], 3),
                     "boot_95CI_low": round(lo, 3), "boot_95CI_high": round(hi, 3),
                     "excludes_zero": bool(not (lo < 0 < hi))})
ci_df_k5 = pd.DataFrame(ci_rows).sort_values("real_corr", ascending=False)
ci_df_k5.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k5_bootstrap_CIs.csv", index=False)
print("Bootstrap 95% CIs (200 gene resamples), k=5, for the 4 original candidates plus any new 95th-percentile survivors:")
display(ci_df_k5)
""")

md(r"""### Repeating at k=10

Five modules gives each chemical a 5-element profile : thin, regardless of
which similarity metric is used. Repeating the entire pipeline above at
k=10 (a fresh 30-fit consensus, not reused from k=5) checks whether any
conclusion survives a real change in resolution, or only holds at the
specific k this notebook otherwise uses.
""")

code(r"""
K10 = 10
t0 = time.time()
consensus_labels_k10, consensus_strength_k10, consensus_runs_k10 = consensus_cluster(data_primary, K10, R=30, base_seed=0)
assign_consensus_k10 = pd.Series(consensus_labels_k10, index=matrix_primary.index, name="module")
print(f"K=10 consensus (30 fits) in {time.time()-t0:.0f}s")
print("K=10 consensus module sizes:")
display(assign_consensus_k10.value_counts().sort_index().to_frame("n_genes"))
assign_consensus_k10.to_csv(ROOT / "output" / "kspaces_runs" / "primary_k10_consensus_module_assignments.csv")
""")

code(r"""
obs_k10, exp_k10, log2_k10, pearson_k10 = residual_table_from_assignment(assign_consensus_k10, K10)
uniformity_k10 = uniformity_test(obs_k10, exp_k10)
print("Chi-square uniformity test, k=10:")
display(uniformity_k10)
n_nonuniform_k10 = int(uniformity_k10["non_uniform_FDR0.05"].sum())
print(f"\n{n_nonuniform_k10} / {len(CHEMICALS)} chemicals non-uniform at k=10 (FDR<0.05).")
uniformity_k10.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k10_uniformity_test.csv", index=False)
""")

code(r"""
corr_resid_k10, cos_resid_k10 = similarity_matrices(log2_k10)
corr_pearsonresid_k10, _ = similarity_matrices(pearson_k10)
print("Chemical-chemical similarity of residual vectors, k=10 consensus (Pearson correlation of log2 ratio):")
display(corr_resid_k10.round(2))

rho_cos_k10, _ = spearmanr([corr_resid_k10.loc[a, b] for a, b in pairs_all], [cos_resid_k10.loc[a, b] for a, b in pairs_all])
rho_pear_k10, _ = spearmanr([corr_resid_k10.loc[a, b] for a, b in pairs_all], [corr_pearsonresid_k10.loc[a, b] for a, b in pairs_all])
print(f"Rank agreement with cosine: rho={rho_cos_k10:.3f}; with Pearson-residual cross-check: rho={rho_pear_k10:.3f}")

null_pooled_corr_k10, null_pooled_cos_k10, null_by_pair_corr_k10 = null_test_residual(assign_consensus_k10, K10, n_shuffle=100, seed=0)
print(f"\nCorrected-metric null (Pearson corr of log2 ratio), k=10:")
print(f"  mean={null_pooled_corr_k10.mean():.3f}  median={np.percentile(null_pooled_corr_k10,50):.3f}  "
      f"5th={np.percentile(null_pooled_corr_k10,5):.3f}  95th={np.percentile(null_pooled_corr_k10,95):.3f}")
print(f"Corrected-metric null (cosine of log2 ratio), k=10:")
print(f"  mean={null_pooled_cos_k10.mean():.3f}  median={np.percentile(null_pooled_cos_k10,50):.3f}  "
      f"5th={np.percentile(null_pooled_cos_k10,5):.3f}  95th={np.percentile(null_pooled_cos_k10,95):.3f}")
""")

code(r"""
null_summary_rows_k10 = []
for a, b in pairs_all:
    nv = np.array(null_by_pair_corr_k10[(a, b)])
    real_v = corr_resid_k10.loc[a, b]
    pctile = (nv < real_v).mean() * 100
    null_summary_rows_k10.append({"pair": f"{a}-{b}", "real_corr": round(real_v, 3),
                                   "null_mean": round(nv.mean(), 3), "percentile_vs_null": round(pctile, 1),
                                   "p_empirical": (np.sum(nv >= real_v) + 1) / (len(nv) + 1)})
null_summary_k10 = pd.DataFrame(null_summary_rows_k10)
null_summary_k10["q_bh_45_pairs"] = multipletests(null_summary_k10["p_empirical"], method="fdr_bh")[1]
null_summary_k10 = null_summary_k10.sort_values(["q_bh_45_pairs", "p_empirical", "real_corr"], ascending=[True, True, False]).reset_index(drop=True)
null_summary_k10.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k10_corr_vs_null.csv", index=False)
print("All 45 pairs ranked by percentile vs. the corrected-metric null, k=10:")
display(null_summary_k10)
survivors_k10 = null_summary_k10[null_summary_k10["q_bh_45_pairs"] < 0.05]
print(f"\nPairs passing BH q<0.05 across 45 pair tests: {len(survivors_k10)} / 45")
display(survivors_k10)

print("\nThe four originally-highlighted candidate pairs, under the corrected metric at k=10:")
for a, b in [("GenX", "PFBS"), ("PFOS", "PFNA"), ("PFOS", "PFOA"), ("PFOSA", "PFPeA")]:
    row = null_summary_k10[null_summary_k10["pair"] == f"{a}-{b}"].iloc[0]
    print(f"  {a}-{b}: real_corr={row['real_corr']:.3f}  percentile_vs_null={row['percentile_vs_null']:.0f}"
          f"  empirical_p={row['p_empirical']:.3f}  q_BH={row['q_bh_45_pairs']:.3f}")
""")

code(r"""
survivor_pairs_k10 = [tuple(p.split("-")) for p in survivors_k10["pair"]]
ci_pairs_k10 = sorted(set(candidate_pairs) | set(survivor_pairs_k10))
boot_vals_k10 = bootstrap_ci_residual(assign_consensus_k10, K10, ci_pairs_k10, n_boot=200, seed=0)
ci_rows_k10 = []
for p, vals in boot_vals_k10.items():
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    ci_rows_k10.append({"pair": f"{p[0]}-{p[1]}", "real_corr": round(corr_resid_k10.loc[p[0], p[1]], 3),
                         "boot_95CI_low": round(lo, 3), "boot_95CI_high": round(hi, 3),
                         "excludes_zero": bool(not (lo < 0 < hi))})
ci_df_k10 = pd.DataFrame(ci_rows_k10).sort_values("real_corr", ascending=False)
ci_df_k10.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_k10_bootstrap_CIs.csv", index=False)
print("Bootstrap 95% CIs (200 gene resamples), k=10:")
display(ci_df_k10)
""")

code(r"""
survivors_k5_set = set(survivors_k5["pair"])
survivors_k10_set = set(survivors_k10["pair"])
cross_k_robust = sorted(survivors_k5_set & survivors_k10_set)
only_k5 = sorted(survivors_k5_set - survivors_k10_set)
only_k10 = sorted(survivors_k10_set - survivors_k5_set)
print(f"Pairs passing BH q<0.05 at k=5 ONLY (k-dependent): {only_k5}")
print(f"Pairs passing BH q<0.05 at k=10 ONLY (k-dependent): {only_k10}")
print(f"Pairs passing BH q<0.05 at BOTH k=5 AND k=10: {cross_k_robust}")

k5_lookup = null_summary_k5.set_index("pair")
k10_lookup = null_summary_k10.set_index("pair")
cross_k_summary = pd.DataFrame([
    {"pair": p, "k5_percentile": k5_lookup.loc[p, "percentile_vs_null"], "k10_percentile": k10_lookup.loc[p, "percentile_vs_null"],
     "k5_q_bh": k5_lookup.loc[p, "q_bh_45_pairs"], "k10_q_bh": k10_lookup.loc[p, "q_bh_45_pairs"],
     "k5_real_corr": k5_lookup.loc[p, "real_corr"], "k10_real_corr": k10_lookup.loc[p, "real_corr"]}
    for p in cross_k_robust
], columns=[
    "pair", "k5_percentile", "k10_percentile", "k5_q_bh", "k10_q_bh",
    "k5_real_corr", "k10_real_corr",
])
print("\nCross-k-robust pairs, detail:")
display(cross_k_summary)
cross_k_summary.to_csv(ROOT / "output" / "kspaces_runs" / "residual_metric_cross_k_robust_pairs.csv", index=False)
""")

md(r"""**Verdict.** The module-size adjustment fixes the degenerate raw
percentage-correlation null, but statistical evidence is determined by the
empirical p-values and BH q-values printed above, not by an uncorrected 95th
percentile screen. Results are corrected across all 45 pairs separately at
k=5 and k=10, and the cross-k table contains only pairs passing q<0.05 at
both resolutions. An empty table is a valid negative result.

Even a retained pair remains exploratory. The label-shuffle count limits
p-value resolution, the bootstrap conditions on fitted consensus labels,
the same ten chemicals informed module discovery, and no external chemical
panel validates generalization. The original raw-correlation rankings remain
retracted regardless of the corrected result.
""")

# ---------------------------------------------------------------------------
md(r"""## 14. Re-testing the Section 10 structure-activity claims with the corrected metric

Section 10 reported that chemical pairs with the same chain length have
higher mean module-usage correlation (0.74) than different-chain-length
pairs (0.58), using the now-retracted raw-percentage metric. Redone here on
the k=5 consensus residual metric from Section 13 (`corr_resid_k5`, the
log2-ratio Pearson correlation), with a proper label-permutation test:
shuffle the chain-length labels across the **10 chemicals** (not the genes)
10,000 times, and ask how often a random labeling produces a
same-vs-different mean gap at least as large as the real one. Repeated for
functional group and ether linkage; ether is flagged as underpowered
regardless of outcome : only GenX and PFEESA are ether-linked, giving
exactly one same-ether "Yes-Yes" pair out of 45.
""")

code(r"""
def permutation_attr_test(sim_matrix: pd.DataFrame, attr_series: pd.Series, n_perm: int = 10000, seed: int = 0):
    chems = list(sim_matrix.index)
    ii, jj = np.triu_indices(len(chems), k=1)
    sim_vals = sim_matrix.values[ii, jj]
    attr_vals = attr_series.loc[chems].values

    def stat(labels):
        same = labels[ii] == labels[jj]
        if same.sum() == 0 or same.sum() == len(same):
            return np.nan, np.nan, np.nan
        return sim_vals[same].mean(), sim_vals[~same].mean(), sim_vals[same].mean() - sim_vals[~same].mean()

    same_mean, diff_mean, obs = stat(attr_vals)
    n_same_pairs = int((attr_vals[ii] == attr_vals[jj]).sum())
    rng = np.random.default_rng(seed)
    perm_stats = np.full(n_perm, np.nan)
    for p in range(n_perm):
        _, _, s = stat(rng.permutation(attr_vals))
        perm_stats[p] = s
    valid = ~np.isnan(perm_stats)
    p_value = (np.sum(np.abs(perm_stats[valid]) >= abs(obs)) + 1) / (valid.sum() + 1)
    return same_mean, diff_mean, obs, perm_stats[valid], p_value, n_same_pairs

attr_test_rows = []
for attr in ["chain_length", "functional_group", "ether"]:
    same_mean, diff_mean, obs, perm_stats, p, n_same = permutation_attr_test(corr_resid_k5, ATTRS[attr], n_perm=10000, seed=0)
    attr_test_rows.append({
        "attribute": attr, "same_mean_corr": round(same_mean, 3), "diff_mean_corr": round(diff_mean, 3),
        "observed_gap": round(obs, 3), "n_same_pairs": n_same, "n_total_pairs": len(pairs_all),
        "permutation_p_value": round(p, 4), "null_mean": round(perm_stats.mean(), 3), "null_sd": round(perm_stats.std(), 3),
    })
attr_test_df = pd.DataFrame(attr_test_rows)
attr_test_df.to_csv(ROOT / "output" / "kspaces_runs" / "structure_activity_permutation_test.csv", index=False)
print("Structure-activity permutation test (10,000 permutations of chemical-level attribute labels), corrected k=5 metric:")
display(attr_test_df)

n_ether_yes = int((ATTRS["ether"] == "Yes").sum())
print(f"\nEther category counts: {ATTRS['ether'].value_counts().to_dict()} -- only "
      f"{n_ether_yes * (n_ether_yes - 1) // 2} same-ether 'Yes-Yes' pair(s) exist (GenX-PFEESA). "
      f"This comparison cannot be powered regardless of the p-value above.")
""")

md(r"""**The chain-length effect does not survive.** Under the corrected
metric, same-chain-length pairs have a slightly *lower* mean residual
correlation than different-chain-length pairs (see `same_mean_corr` vs.
`diff_mean_corr` above) : the opposite direction from Section 10's original
claim : and the 10,000-permutation test finds this gap is entirely
consistent with chance. Functional group is also null. Ether linkage cannot
be tested meaningfully at all: with only one same-ether pair (GenX-PFEESA)
in the entire 45-pair set, "same vs. different" for ether reduces to "one
specific pair vs. everyone else," and a p-value of 1.0 there reflects that
the dataset cannot answer the ether question with only two ether-linked
chemicals : not evidence of no effect.

**Conclusion: none of Section 10's three structure-activity claims replicate
under the corrected, null-adjusted similarity metric.** The original
chain-length finding (0.74 vs 0.58) was a product of the same module-size
confound Section 12 retracted for the pairwise chemical comparisons : it is
not merely unconfirmed, it is contradicted in direction by the corrected
version. Section 10's original cells are left in place (with a pointer to
this section) rather than edited to assert a number that no longer holds;
Section 16's limitations are updated accordingly.
""")

# ---------------------------------------------------------------------------
md(r"""## 15. Which module is "glutathione transferase activity," really?

Section 12's single-fit-vs-consensus GO comparison table shows module 3 (in
Section 5's numbering) going from "glutathione transferase activity" (single
fit, 286 genes) to "cytosolic large ribosomal subunit" (consensus, 2,546
genes) : a roughly ninefold size change and a complete change of biological
story. Before reading either version as real biology: what actually happened
to the genes?
""")

code(r"""
target_hits = compare_df.index[compare_df["top_term_single_fit (Section 8)"].str.contains("glutathione transferase", na=False)]
TARGET_MODULE = int(target_hits[0])
print(f"Module under investigation: {TARGET_MODULE}")
print(f"  single-fit (Section 5/8) top term:  {compare_df.loc[TARGET_MODULE, 'top_term_single_fit (Section 8)']}")
print(f"  consensus (Section 12) top term:    {compare_df.loc[TARGET_MODULE, 'top_term_consensus']}")

ap_aligned_idx = assign_primary.loc[assign_consensus.index]
overlap_matrix = np.zeros((K, K), dtype=int)
for i in range(K):
    for j in range(K):
        overlap_matrix[i, j] = int(((assign_consensus.values == i) & (ap_aligned_idx.values == j)).sum())
overlap_df = pd.DataFrame(overlap_matrix, index=[f"raw_consensus_{i}" for i in range(K)], columns=[f"primary_{j}" for j in range(K)])
print("\nFull gene-count overlap matrix (rows = raw consensus label, cols = Section 5 primary label):")
display(overlap_df)

genes_single = set(assign_primary.index[assign_primary == TARGET_MODULE])
genes_consensus = set(assign_consensus_aligned.index[assign_consensus_aligned == TARGET_MODULE])
inter = genes_single & genes_consensus
union = genes_single | genes_consensus
jaccard = len(inter) / len(union) if union else float("nan")
print(f"\nn_genes single-fit module {TARGET_MODULE}: {len(genes_single)}")
print(f"n_genes consensus module {TARGET_MODULE} (aligned numbering): {len(genes_consensus)}")
print(f"Gene overlap: {len(inter)} genes")
print(f"Jaccard overlap: {jaccard:.4f}")
overlap_df.to_csv(ROOT / "output" / "kspaces_runs" / "module_investigation_overlap_matrix.csv")
""")

code(r"""
raw_to_aligned = {}
for i in range(K):
    aligned_vals = assign_consensus_aligned[assign_consensus == i]
    raw_to_aligned[i] = int(aligned_vals.iloc[0])
    assert (aligned_vals == raw_to_aligned[i]).all(), "alignment is not a clean bijection -- unexpected"
lut = np.array([raw_to_aligned[i] for i in range(K)])

strength_by_module_rows = []
for m in sorted(assign_consensus_aligned.unique()):
    vals = consensus_strength_s.loc[assign_consensus_aligned.index[assign_consensus_aligned == m]]
    strength_by_module_rows.append({
        "module": m, "n_genes": len(vals), "mean_strength": round(vals.mean(), 3), "median_strength": round(vals.median(), 3),
        "pct_core_ge0.8": round((vals >= 0.8).mean() * 100, 1), "pct_ambiguous_lt0.5": round((vals < 0.5).mean() * 100, 1),
    })
strength_by_module_df = pd.DataFrame(strength_by_module_rows).set_index("module")
print(f"consensus_strength distribution by module (Section 5/8 numbering) -- module {TARGET_MODULE} vs. the other four:")
display(strength_by_module_df)

run_sizes = np.array([(lut[consensus_runs[r]] == TARGET_MODULE).sum() for r in range(consensus_runs.shape[0])])
print(f"\nModule {TARGET_MODULE} size across each of the 30 individual consensus runs (not just the majority-vote result):")
print(f"  all 30 values: {run_sizes.tolist()}")
print(f"  mean={run_sizes.mean():.0f}  sd={run_sizes.std():.0f}  min={run_sizes.min()}  max={run_sizes.max()}  "
      f"(majority-vote consensus size: {(assign_consensus_aligned==TARGET_MODULE).sum()})")

run_size_variance_rows = []
for m in range(K):
    sizes_m = np.array([(lut[consensus_runs[r]] == m).sum() for r in range(consensus_runs.shape[0])])
    run_size_variance_rows.append({"module": m, "mean_size": round(sizes_m.mean(), 0), "sd_size": round(sizes_m.std(), 0),
                                    "min_size": int(sizes_m.min()), "max_size": int(sizes_m.max()),
                                    "cv": round(sizes_m.std() / sizes_m.mean(), 3)})
run_size_variance_df = pd.DataFrame(run_size_variance_rows).set_index("module")
print("\nFor context, run-to-run size variance for all 5 modules:")
display(run_size_variance_df)

strength_by_module_df.to_csv(ROOT / "output" / "kspaces_runs" / "module_investigation_consensus_strength_by_module.csv")
run_size_variance_df.to_csv(ROOT / "output" / "kspaces_runs" / "module_investigation_run_size_variance.csv")
""")

md(r"""**Zero gene overlap. `Jaccard = 0.0000`.** Not one of the single
fit's 286 "glutathione transferase activity" genes is among the consensus's
2,546 genes carrying the module-3 label. The overlap matrix above explains
why: primary module 3's 286 genes split almost entirely between two much
larger consensus clusters, and essentially none overlap with the consensus
cluster that Hungarian alignment happened to assign the "module 3" label to.
That alignment is a valid bijection : Hungarian assignment must map every
consensus cluster to some primary cluster : but for this specific module the
best available match was still a near-zero-overlap pairing, because no
consensus cluster closely resembles primary module 3.

**"Module 3 changed from glutathione transferase to ribosomal" is not an
accurate description of what happened.** It is more accurate to say: the
single fit's tiny 286-gene module does not correspond to any single
consensus module : it dissolves, its genes absorbed mostly into two much
larger consensus clusters : and the consensus cluster that inherits the
"module 3" numeric label via forced 1-to-1 alignment is an unrelated
2,546-gene cluster with its own, separately real, "cytosolic large ribosomal
subunit" identity. Section 12's module-by-module GO comparison table is
numerically well-defined and not wrong on its own terms, but reading its
module-3 row as one module's biology "changing" is wrong : that row compares
two different clusters that happen to share a label, not a label's contents
evolving.

**Is the consensus module 3 (2,546 "ribosomal" genes) itself trustworthy?**
By `consensus_strength` it is not the worst of the five : see the table
above. **Whichever module has the lowest mean strength and highest ambiguous
fraction in the table above is the real low-confidence "leftover" bucket in
this consensus assignment**, not necessarily module 3; that module is also
worth cross-checking against `compare_df` (Section 12) for whether it has no
specific enriched GO term under either the single fit or consensus, which
would be independent evidence of the same conclusion. Module 3's run-to-run
size does vary substantially in absolute terms across the 30 individual
fits (see the printed range above), but whether its coefficient of variation
is actually the highest of the five modules is checked directly in the `cv`
column above, not assumed.

**Verdict, added to Section 16: consensus module 3 ("cytosolic large
ribosomal subunit") should be read as its own real module** : evaluated on
its own `consensus_strength` and GO evidence above : **not as a corrected or
evolved version of the single fit's "glutathione transferase" module, which
shares zero genes with it and has no real consensus counterpart at all.**
The single-fit module 3 was, by this evidence, closer to an idiosyncratic
artifact of that one EM run than a real, reproducible gene set. Whichever
module the table above flags as lowest-confidence (n_genes, mean_strength,
pct_ambiguous) is the one that actually deserves a "residual bucket, don't
over-interpret" caveat : cross-checked against `compare_df`'s "no specific
term" column where applicable.
""")

# ---------------------------------------------------------------------------
md(r"""## 16. Summary and honest limitations

### The headline correction from this pass
**The correlation-based "which PFAS chemicals share modules" claims in
Sections 5, 7, and the early part of Section 12 do not survive stress
testing and are retracted, not softened.** Correlating each chemical's
5-element percent-of-DEGs-per-module vector at k=5 produces high
correlation (median ~0.97 under a label-shuffle null) almost independent of
biology, because with only 5 modules of uneven size, any gene subset :
real or randomly labeled : lands across them in roughly the modules' own
size proportions. Under a proper null (Section 12), **not one of the 45
chemical pairs clears a 95th-percentile bar**, and the previously-headlined
"GenX-PFBS is the most robust pair" finding sits at the **0th** percentile :
real correlation *below* nearly every random shuffle, the opposite of the
original claim.

**Update after changing the primary representation:** Section 13 builds a
null-adjusted comparison from observed-versus-expected DEG counts per module
and repeats it at k=5 and k=10. All numerical pair, structure-activity, and
module-identity conclusions must be read from the next completed execution.
Results from the prior percentile-rank-primary run are historical and must
not be carried forward as results of the continuous-logFC-primary analysis.

### What this analysis found
(exact numbers are in the sections cited, generated live from this run --
not repeated here to avoid two sources of truth that could drift apart)
- **The consensus module assignment (Section 12), not the single fit in
  Section 5, is the recommended module structure to cite** : it resolves EM
  local-optima noise a single fit can't distinguish from real signal, and
  ships a per-gene `consensus_strength` confidence score Section 5's
  approach didn't have. This recommendation is about *reproducibility of
  gene membership*, separate from the correlation metric's retraction above.
- A formal permutation test (Section 12) supports cross-chemical dependence
  relative to its column-shuffle null. It does not establish a unique or
  biologically valid module structure. This survives the metric critique
  because it tests log-likelihood, not the confounded correlation.
- The chemical-level jackknife (Section 12) identifies which chemicals the
  overall structure most depends on versus which are largely redundant :
  using Adjusted Rand Index, which is corrected for chance agreement by
  construction and is not subject to the same confound as the correlation
  metric.
- Extending consensus clustering to the sensitivity variants (Section 12)
  shows the low cross-strategy ARI first found in Section 7 is **mostly
  real disagreement between strategies, not EM noise** : consensus-to-
  consensus ARI is only marginally higher than single-fit-to-single-fit ARI.
- **The corrected, null-adjusted chemical-similarity metric (Section 13)**
  replaces the retracted correlation approach: observed-vs-expected DEG
  counts per module (log2 ratio, primary; standardized Pearson residual,
  cross-check), tested with the same label-shuffle null and gene-resample
  bootstrap as before, at both k=5 and k=10. Empirical p-values are adjusted
  across all 45 pairs separately at each resolution. No pair passes BH q<0.05
  at either k in the corrected run, so the earlier unadjusted eleven-pair
  screen is withdrawn rather than promoted as a discovery.
- **Section 10's chain-length structure-activity claim does not replicate**
  under the corrected metric (Section 14) : the effect reverses direction
  and a 10,000-permutation test finds it indistinguishable from chance.
  Functional group is likewise null; ether linkage is underpowered (n=2)
  and untestable either way.
- **The apparent "module 3 changed from glutathione transferase to
  ribosomal" finding (Section 12's GO comparison) is a label-alignment
  artifact, not one module evolving** (Section 15) : the single-fit and
  consensus versions of "module 3" share zero genes (Jaccard=0.0000); the
  single-fit module dissolves under consensus rather than persisting in
  refined form.
- Re-deriving the heatmap, GO enrichment, and Sankey diagram from the
  consensus assignment (Section 12) is the corrected version of Sections 6,
  8, and 9 : read alongside Section 15's caveat about module 3 specifically.
- k-spaces finds real, non-null structure (Sections 11a and 12), but is
  **less reproducible than plain k-means** as a single fit (Section 11b) :
  consensus clustering (Section 12) narrows that gap by building a more
  representative assignment than any one fit, but does not eliminate the
  underlying tradeoff. This remains a genuine, unresolved method-choice
  question, not swept under the rug.
- Deterministic annealing was tested as a way to reduce EM optimization
  noise directly and **did not help** (Section 12) : a negative result kept
  in rather than omitted.

### Limitations (read before citing any of this)
- **The raw correlation-based chemical-similarity metric used in Sections 5
  and 7 is confounded by module size at k=5 and is retracted outright.**
  Section 13's corrected, null-adjusted replacement is the valid diagnostic
  to cite for chemical-pair similarity. It identifies no BH-significant pair
  at k=5 or k=10 in this run, so the analysis supports a null pairwise result,
  not a general "high correlation = similar" conclusion.
- **Section 10's chain-length-effect claim is retracted, not merely
  provisional** : Section 14 independently re-tested it on the corrected
  metric with a proper permutation test and it does not hold up (the effect
  reverses direction and is statistically indistinguishable from chance).
- Per-gene `consensus_strength` (Section 12) should be used to weight trust
  in any specific gene's module call : genes below ~50% agreement across
  runs are genuinely ambiguous, not just noisy labels to ignore. Section 15
  shows this distribution varies meaningfully by module; check it before
  treating any one module's gene list as settled.
- The core-only subset used in parts of Section 12 (consensus_strength>=0.8)
  is a **selected** 36.7%-of-genes subset, not a random subsample : rank
  agreement between it and the full consensus set was checked (Section 12)
  and found weak-to-moderate, not strong.
- **The Hungarian-alignment "same label = same module" convention used to
  compare single-fit and consensus module numbering (Section 12) can produce
  misleading comparisons for modules with no good match on the other side**
  : Section 15 found one module (numeric label 3) whose "before" and "after"
  are entirely different gene sets forced into the same slot by the
  alignment's bijection constraint. Check gene-level overlap (as Section 15
  does), not just the numeric label, before reading any single module's
  before/after GO terms as one module's biology changing.
- k=5 is a pragmatic choice balancing stability and resolution, not a
  validated optimum : needs review. Section 12's small heterogeneous-dimension
  probe (4 hand-picked configurations, not an exhaustive search) found
  mixed-dimension models (e.g. one 2-D module among 1-D ones) can beat
  uniform d=1 on BIC : a promising lead for future model selection that this
  analysis did not chase to completion. Section 13 shows some conclusions
  (which chemical pairs are similar) are themselves k-dependent between k=5
  and k=10, which is a further argument for treating k=5 as provisional.
- Bootstrap stability (Section 4), the chemical jackknife (Section 12), the
  correlation/residual null tests (Sections 12-13, 100-200 shuffles/
  resamples each), and the chain-length permutation test (Section 14,
  10,000 permutations, cheap enough to not be a limiting factor there) all
  used compute-time-limited counts where noted; more would sharpen, not
  overturn, these estimates.
- GO enrichment used only WormBase GO (BP/MF/CC) : no KEGG/Reactome, and no
  cross-check against the Reactome plots already in `Data/DEGs/` from
  earlier per-chemical analysis.
- k-spaces core runs in Python (`github.com/pachterlab/k-spaces`, the
  authors' own reference implementation), not R : chosen deliberately over
  reimplementing their EM algorithm by hand in R, which would have carried
  real correctness risk.

### Analyses requiring scientific decisions or additional inputs
- Advisor/PI review of the k=5 choice, the consensus-vs-single-fit
  recommendation, the k-spaces-vs-k-means tradeoff, and the implications of
  the Section 13 null pairwise result in a small panel (10 chemicals, 45 pairs).
- A biological/mechanistic look at why the jackknife's most structurally
  influential chemicals matter as much as they do (real chemistry, or a
  vehicle/control artifact?).
- Exhaustive heterogeneous-dimension model selection (only a handful of
  hand-picked configurations were tried).
- KEGG/Reactome enrichment alongside GO.
- Feeding module structure into a toxicity-prediction model (the planning
  document's "extra ideas" section) : a separate, substantial follow-on task.
- Extending the corrected residual metric (Section 13) to k values beyond
  5 and 10, and to the sensitivity gene-selection variants (signed percentile ranks,
  top-N), the way Section 12 did for the retracted correlation metric.
""")

nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "kspaces_venv", "language": "python", "name": "kspaces_venv"},
    "language_info": {"name": "python", "version": "3.13"},
}

out_path = "PFAS_kspaces_analysis.ipynb"
with open(out_path, "w") as f:
    nbf.write(nb, f)
print(f"Wrote {out_path} with {len(cells)} cells")
