# PFAS k-spaces analysis: corrected findings

This summary is generated from the corrected workflow represented by `build_notebook.py` and the executed `PFAS_master_analysis.ipynb`. The executed master notebook is authoritative when this summary and a stored artifact disagree.

## Data and primary representation

The ten edgeR contrast tables share an identical 13,852-gene measured universe after identifier normalization. The primary analysis restricts to the 10,255 genes detected at FDR below 0.05 in at least one contrast, retains continuous log2 fold changes, and row-standardizes each gene across the ten chemical conditions. Signed within-chemical percentile ranks and top-N selections are sensitivity analyses, not the primary representation.

The contrast tables are not raw count matrices. They support cross-condition pattern analysis but cannot reconstruct replicate-level variance, dose response, time response, or causal effects.

## Module discovery and uncertainty

The recommended k=5 assignment is the majority vote from 30 independently seeded k-spaces fits after label alignment. Mean per-gene consensus strength is 0.951; 91.4% of genes have at least 80% agreement. These values measure optimizer agreement under the declared model, not external biological reproducibility.

Dropping one chemical changes the partition to different degrees. PFOSA is most influential in the current jackknife, with ARI 0.547 against the full ten-condition consensus. The remaining leave-one-condition ARIs range from 0.766 to 0.886. This sensitivity must accompany any module interpretation.

## Retracted raw similarity claim

The original correlation of five-element percent-of-DEGs-per-module profiles is not a valid chemical-similarity statistic. Its label-shuffle null is intrinsically high because uneven module sizes dominate the profile. Raw pair rankings and their legacy figures are retained only for audit and must not be cited as findings.

The corrected exploratory metric subtracts the module-size expectation before comparing residual profiles. Pairwise empirical p-values are adjusted across all 45 chemical pairs separately at k=5 and k=10. No pair meets BH q below 0.05 at either k in this run, so `output/kspaces_runs/residual_metric_cross_k_robust_pairs.csv` is an empty, header-only result. The earlier eleven-pair screen used unadjusted pair-specific thresholds and is withdrawn. The analysis does not identify a statistically supported pairwise chemical-similarity claim.

## Functional interpretation

Consensus-module GO enrichment uses, for each namespace, the intersection of the 10,255-gene clustered DEG union and genes mapped in the pinned WormBase annotation. Raw p-values across all module, namespace, and term tests are corrected together with one Benjamini-Hochberg family. The consensus result contains 216 displayed significant module-term records after this correction. Chemical-exclusive and chemical-shared studies instead use the measured mapped universe, while chemical-by-module studies use the mapped genes within their corresponding consensus module; each of those three declared families receives its own family-wide correction.

Examples of leading Biological Process terms are autophagy for module 0, nucleic-acid and RNA metabolism for module 1, translation for module 2, xenobiotic metabolism and cilium organization for module 3, and innate immune or defense response for module 4. These are enrichment associations. Annotation incompleteness and module-assignment uncertainty are not fully propagated into their p-values.

The historical claim that a corrected top module was labeled “glutathione transferase activity” does not reproduce under the continuous-effect primary matrix and family-wide GO procedure. The master notebook records that absence instead of forcing the old module number onto the new result.

## Modeling and analysis planning

PyTorch and TensorFlow implement the same compact synthetic workflow-classification contract. The executed notebook regenerated both backends with 2,092 parameters each. On the 16-record held-out smoke set, PyTorch achieved accuracy 0.875 and macro F1 0.861; TensorFlow achieved 1.0 for both. This narrow, synthetic check detects implementation-specific problems and does not support a PFAS biological claim or framework-superiority claim.

Natural-language requests pass through typed interpretation, dataset compatibility checks, verified evidence retrieval, allowlisted planning, parameter and scientific validation, unsupported-claim detection, and a human-review boundary. The notebook demonstrates a review-ready plan and an abstention but does not fabricate approval or execute a live scientific request. The adapted language-model artifact is disclosed as loss-only support and is not the execution authority.

## Reproduction

Never edit either generated notebook directly. Regenerate and execute the master notebook from this directory:

```bash
.venv/bin/python build_notebook.py
.venv/bin/python build_master_notebook.py
.venv/bin/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=3600 \
  --ExecutePreprocessor.kernel_name=kspaces_venv \
  PFAS_master_analysis.ipynb
```

The scientific environment is pinned in `requirements.txt`. `GLOBAL_SEED=0` controls the NumPy global state used by k-spaces. Exact replay evidence and inferential stability are reported separately.

## Analyses requiring scientific decisions or additional inputs

- Advisor or PI review of k=5, the consensus recommendation, and the null pairwise result.
- Raw counts and complete sample metadata for replicate-aware reanalysis.
- External data for biological replication and causal or regulatory claims.
- A future pairwise testing plan that propagates clustering uncertainty and is powered for ten conditions.
- Pathway resources beyond the pinned WormBase GO annotation, with their own provenance and testing family.
