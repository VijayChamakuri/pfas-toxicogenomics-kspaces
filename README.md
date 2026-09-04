# PFAS toxicogenomics with consensus k-spaces

[![CI](https://github.com/VijayChamakuri/pfas-toxicogenomics-kspaces/actions/workflows/integrated-pfas-ci.yml/badge.svg)](https://github.com/VijayChamakuri/pfas-toxicogenomics-kspaces/actions/workflows/integrated-pfas-ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-3776AB.svg)](https://www.python.org/)
[![Research status](https://img.shields.io/badge/status-exploratory-orange.svg)](#scientific-scope-and-limitations)

A reproducible analysis of shared and chemical-associated transcriptional responses to ten PFAS
exposures in *Caenorhabditis elegans*. The project combines edgeR contrast outputs, consensus
k-spaces module discovery, uncertainty analysis, Gene Ontology enrichment, and an evidence-grounded
workflow system with PyTorch and TensorFlow engineering checks.

> **Primary scientific artifact:**
> [`K-spaces Analysis/PFAS_master_analysis.ipynb`](K-spaces%20Analysis/PFAS_master_analysis.ipynb)
> is executed end to end and contains the complete analytical narrative, code, tables, and figures.

![Consensus k=5 gene-module heatmap](K-spaces%20Analysis/output/figures/consensus_k5_heatmap.png)

## Research question

Can transcriptional responses across structurally diverse PFAS exposures be separated into robust
shared modules and chemical-associated patterns that support group-based toxicological reasoning?

Ten PFAS chemicals were analyzed at their experiment-specific EC50 values: PFOSA, PFBSA, PFOS,
PFNA, PFOA, GenX, PFEESA, PFBS, PFPeA, and PFBA. The starting point is ten precomputed
PFAS-versus-control edgeR tables sharing the same filtered 13,852-gene universe.

## Study design and analytical contract

| Component | Scope |
|---|---|
| Organism | *C. elegans* |
| Conditions | 10 PFAS exposures, each compared with its control |
| Measured universe | 13,852 genes in every contrast |
| Clustered universe | 10,255 genes with FDR below 0.05 in at least one contrast |
| Primary representation | Continuous log2 fold change, row standardized across chemicals |
| Module recommendation | k=5 consensus from 30 independently seeded, label-aligned fits |
| Functional annotation | Gene Ontology BP, MF, and CC with declared eligible backgrounds |
| Reproducibility | Fixed global seed, pinned requirements, input hashes, executed notebook, tests, CI |

Signed percentile ranks and top-N selections are sensitivity analyses. Detection thresholds are not
used to erase continuous effects in the primary matrix.

## Verified findings

- PFOSA produced the largest detected response with 9,161 DEGs, followed by PFPeA with 5,997 and
  PFBSA with 2,601. PFOA had the smallest detected response with 154 DEGs.
- The recommended k=5 consensus modules contain 1,915, 1,496, 2,386, 2,944, and 1,514 genes.
- Mean per-gene consensus strength is 0.951, and 91.4% of clustered genes have at least 0.8 agreement
  across the 30 fits. This measures optimizer agreement, not external biological replication.
- Chemical jackknifing identifies PFOSA as the most influential condition, with ARI 0.547 against
  the full consensus. The remaining leave-one-chemical ARIs range from 0.766 to 0.886.
- **No chemical pair passes BH q below 0.05 at either k=5 or k=10 across the 45 pairwise tests.**
  The earlier unadjusted eleven-pair screen is withdrawn.
- Consensus GO enrichment contains 216 displayed family-wide significant module-term records.
  Leading Biological Process themes include autophagy, RNA and nucleic-acid metabolism, translation,
  xenobiotic metabolism and cilium organization, and innate immune or defense responses.

These are descriptive and exploratory associations. They do not establish causal mechanisms,
human-health effects, regulatory groups, or framework superiority.

## Results gallery

### Exact cross-chemical DEG intersections

The UpSet-style view shows the 30 largest exact detection patterns. Complete unsigned, upregulated,
and downregulated intersection tables are available under `output/deep_analysis/`.

![Exact DEG intersection patterns](K-spaces%20Analysis/output/figures/deep_analysis/upset_exact_DEG_intersections.png)

### Response scale, exclusive genes, and shared genes

![Per-chemical DEG composition](K-spaces%20Analysis/output/figures/deep_analysis/per_chemical_DEG_composition.png)

### Consensus module activity across ten PFAS conditions

The heatmap and PCA summarize condition-level activity using the consensus assignment. The ten
conditions are observations in the PCA; genes are not treated as independent chemical replicates.

![Module activity and chemical PCA](K-spaces%20Analysis/output/figures/deep_analysis/module_activity_and_PCA.png)

### Chemical-by-module functional depth

All 50 chemical-by-module gene sets are exported. GO testing records every eligible or skipped
study across three namespaces and applies one BH correction across each declared family.

![Significant GO terms per chemical and module](K-spaces%20Analysis/output/figures/deep_analysis/GO_per_chemical_per_module_heatmap.png)

### Ontology-backed module interpretation

These panels use actual GO parent relationships rather than manually assigned pathway trees.

![Consensus-module Gene Ontology relationships](K-spaces%20Analysis/output/figures/deep_analysis/ontology_trees_consensus_modules.png)

Additional figures include signed UpSet plots, chemical-module preference, per-module contribution,
exclusive/shared ontology panels, and Sankey summaries in
[`K-spaces Analysis/output/figures/deep_analysis/`](K-spaces%20Analysis/output/figures/deep_analysis/).

## Important scientific correction

The original correlation between five-element module-percentage profiles is retracted. A label-null
stress test showed that high correlations arise largely from uneven module sizes. The replacement
analysis subtracts module-size expectation and applies empirical pairwise tests with BH correction
across all 45 pairs separately at k=5 and k=10. No pair survives at either resolution.

This negative result is retained as a finding. The repository does not promote unadjusted rankings,
the earlier GenX-PFBS emphasis, or the chain-length structure-activity claim. Reactome results are
not reconstructed because a pinned gene-set release and reproducible historical result table were
not available.

## Repository structure

```text
.
├── Data/DEGs/                         # 10 approved edgeR contrast tables
├── K-spaces Analysis/
│   ├── PFAS_master_analysis.ipynb     # executed, single-entry scientific report
│   ├── build_master_notebook.py       # notebook generator
│   ├── build_notebook.py              # maintained canonical cell source
│   ├── src/pfas_master/               # tested scientific helpers
│   ├── tests/                         # scientific invariants and QA
│   ├── data_external/                 # frozen GO and compressed WormBase snapshots
│   └── output/                        # selected tables, gene lists, PNGs, and PDFs
├── integrated_system/
│   ├── src/pfas_workflow/             # typed planning and execution boundary
│   ├── llm_system/                    # RAG, evaluation, and bounded adaptation code
│   ├── ml_frameworks/                 # PyTorch/TensorFlow parity contract
│   ├── tests/                         # unit and integration tests
│   ├── docker/                        # API and notebook images
│   └── infra/                         # approval-gated AWS design and Terraform
└── .github/workflows/                 # automated verification
```

## Quick start

Python 3.13 is recommended.

```bash
git clone https://github.com/VijayChamakuri/pfas-toxicogenomics-kspaces.git
cd pfas-toxicogenomics-kspaces

python3.13 -m venv "K-spaces Analysis/.venv"
"K-spaces Analysis/.venv/bin/pip" install -r "K-spaces Analysis/requirements.txt"
"K-spaces Analysis/.venv/bin/pip" install nbformat nbconvert jupyterlab ipykernel pytest
"K-spaces Analysis/.venv/bin/python" -m ipykernel install --user \
  --name kspaces_venv --display-name "Python (PFAS k-spaces)"

python3.13 -m venv integrated_system/.venv
integrated_system/.venv/bin/pip install -e './integrated_system[dev,ml]'

"K-spaces Analysis/.venv/bin/python" "K-spaces Analysis/scripts/prepare_annotations.py"
```

Open the finished report:

```bash
"K-spaces Analysis/.venv/bin/jupyter" lab "K-spaces Analysis/PFAS_master_analysis.ipynb"
```

## Rebuild the master notebook

```bash
cd "K-spaces Analysis"
PFAS_SYSTEM_PYTHON="../integrated_system/.venv/bin/python" \
  .venv/bin/python build_master_notebook.py
PFAS_SYSTEM_PYTHON="../integrated_system/.venv/bin/python" \
  .venv/bin/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=5400 \
  --ExecutePreprocessor.kernel_name=kspaces_venv \
  PFAS_master_analysis.ipynb
```

The full run typically takes 15 to 20 minutes on the development machine. Runtime varies by CPU,
memory, BLAS implementation, and package platform.

## Verification

```bash
"K-spaces Analysis/.venv/bin/python" -m pytest -q "K-spaces Analysis/tests"
(cd integrated_system && .venv/bin/ruff check .)
(cd integrated_system && .venv/bin/mypy src/pfas_workflow llm_system)
(cd integrated_system && .venv/bin/pytest -q)
```

The published notebook was accepted only after all 56 code cells executed with zero error outputs,
10 scientific tests passed, and 35 integrated-system tests passed. CI additionally checks linting,
typing, workflow evaluation, both ML backends, Terraform syntax, and container health.

## Engineering extensions

- **PyTorch and TensorFlow:** equivalent compact classifiers are run on a separate 80-record synthetic
  workflow-request dataset. Their 16-record held-out scores are smoke tests for implementation parity,
  not PFAS prediction accuracy.
- **RAG:** reviewed local findings and methods are chunked with source hashes, retrieved for supported
  requests, and passed through typed validation and an abstention boundary.
- **Fine-tuning:** a bounded synthetic router adaptation is documented through a model card and metrics.
  It is not the scientific execution authority.
- **Docker and CI/CD:** multi-stage API/notebook images and GitHub Actions provide repeatable checks.
- **AWS:** Terraform is an approval-gated deployment starter. No live cloud deployment is claimed.

## Scientific scope and limitations

- Raw RNA-seq counts, full sample metadata, and sequencing QC are not included. The upstream edgeR
  analysis cannot be independently reconstructed from this repository.
- GenX used a different vehicle or control, confounding direct cross-chemical interpretation.
- Each chemical was studied at its own EC50, so results do not represent a common external dose.
- There is no dose-response or time-course analysis.
- k=5 is a pragmatic resolution/stability recommendation, not a uniquely proven model.
- Consensus strength captures computational agreement, not biological reproducibility.
- GO p-values do not propagate annotation incompleteness or module-assignment uncertainty.
- No external dataset validates the modules or functional themes.
- Findings must not be interpreted as causal mechanisms or direct evidence for human regulation.

See [`FINDINGS.md`](K-spaces%20Analysis/FINDINGS.md), the
[`scientific audit`](integrated_system/docs/SCIENTIFIC_AUDIT.md), and
[`data-use statement`](integrated_system/docs/DATA_USE.md) before citing results.

## Data and software provenance

The ten approved DGE tables, GO snapshot, WormBase annotation snapshot, versions, checksums, and
third-party terms are documented in [`DATA_SOURCES.md`](DATA_SOURCES.md). GO data are redistributed
under CC BY 4.0. The project code is MIT licensed. The k-spaces package is BSD 2-Clause licensed.

## Citation

Use [`CITATION.cff`](CITATION.cff) to cite this software release. The clustering method is based on:

> Markarian N, Engelhardt BE, Pierce NA, Sternberg PW, Pachter L. *k-spaces: Mixtures of Gaussian
> latent variable models*. bioRxiv. 2025. doi:10.1101/2025.11.24.690254.

When using the ontology or annotations, cite the Gene Ontology Consortium and state the data release
recorded in [`DATA_SOURCES.md`](DATA_SOURCES.md).

## License

Original code is available under the [MIT License](LICENSE). Vendored ontology and annotation data
retain their upstream terms. See [`DATA_SOURCES.md`](DATA_SOURCES.md).
