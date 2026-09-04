# Reproducibility record

## Reproducibility contract

Every reported result should be traceable to immutable inputs, a source revision, an environment, declared parameters, random seeds, and a machine-readable output. Generated notebooks are built from source scripts and executed end to end. Python modules hold reusable logic and tests exercise the same interfaces used by the notebook and service.

## Verified input state

All ten canonical differential-expression files contain 13,852 rows, no missing values, and no duplicate WormBase identifiers after identifier normalization. Their filtered gene universe is identical across contrasts. GenX labels its identifier column `WB_id`; the other nine label it `Unnamed: 0`. The loader normalizes this schema difference without modifying the source files.

The original raw data and result artifacts remain outside this package and are treated as read-only inputs. The legacy ZIP inventory is preserved as a coverage manifest so each archived output has an explicit retain, replace, or diagnostic-only disposition.

## Environment separation

The scientific k-spaces analysis uses `K-spaces Analysis/.venv` and the `kspaces_venv` Jupyter kernel. The analysis-planning package uses `integrated_system/.venv`. This separation prevents framework dependencies from silently changing the validated scientific environment.

Record exact packages before a release:

```bash
"../K-spaces Analysis/.venv/bin/python" -m pip freeze \
  > scientific-environment.txt
.venv/bin/python -m pip freeze > system-environment.txt
```

Environment files are evidence snapshots and should be reviewed before publication because they may include local editable paths.

## Deterministic controls

- The k-spaces workflow sets a global NumPy seed because the package draws from NumPy's global random state.
- Consensus assignment uses 30 independent fits and reports per-gene assignment strength.
- Framework comparison uses shared group-disjoint splits, preprocessing, capacity, seeds, stopping rules, and prediction schemas.
- Workflow evaluation uses versioned, hand-authored cases and deterministic scoring where possible.
- Artifact execution records SHA-256 hashes instead of relying only on filenames.

Determinism does not remove statistical uncertainty. Optimizer variation, observation resampling, clustering uncertainty, and multiple-testing scope are reported separately.

## Verification commands

```bash
cd pfas-toxicogenomics-kspaces/integrated_system
.venv/bin/ruff check .
.venv/bin/python -m pytest -q
PYTHONPATH=src .venv/bin/python scripts/run_evals.py

cd "../K-spaces Analysis"
.venv/bin/python -m pytest -q tests/test_master_science.py
.venv/bin/python build_master_notebook.py
.venv/bin/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=kspaces_venv PFAS_master_analysis.ipynb
```

Full notebook execution is intentionally separate from fast fixture-based CI. A release record should include test counts, execution duration, notebook error count, artifact hashes, environment fingerprints, and any skipped checks.

## Reproducibility matrix

| Component | Inputs | Deterministic controls | Verification | Current boundary |
|---|---|---|---|---|
| DGE table validation | Ten canonical contrast CSV files | Schema, row-count, identifier, and hash checks | `pytest -q tests/test_master_science.py` in the scientific environment | Begins from precomputed edgeR results |
| Module analysis | Validated continuous logFC matrix | Fixed global seed, 30-fit consensus, confidence scores | Clean execution of `PFAS_master_analysis.ipynb` | Modules and pair results remain exploratory where stated |
| Workflow planning | Typed request plus artifact catalog | Allowlisted operations, deterministic rules, versioned evidence | `pytest -q tests/test_workflow.py tests/test_llm_rag.py` | Human review remains mandatory |
| Backend verification | Versioned synthetic request dataset | Shared splits, preprocessing, capacity, seed, and metrics | `python -m ml_frameworks --smoke --framework both` | Software portability evidence only |
| Model adaptation | Synthetic grouped workflow records | Duplicate checks, group-disjoint splits, fixed seed | Training metrics plus schema-level held-out evaluation | Not approved as an execution gate |
| Containers | Source tree and pinned direct requirements | Multi-stage build and nonroot runtime | API health check and notebook image smoke | Requires a running Docker daemon |
| AWS design | Terraform variables and approved infrastructure identifiers | Versioned plan, provider lock, least-privilege review | `terraform fmt`, `init -backend=false`, and `validate` | No deployment or cloud result is claimed |

## Irreducible limits

The RNA-seq differential-expression analysis cannot be independently reproduced from this folder because raw counts, complete sample metadata, QC artifacts, and the original R environment are absent. The exploratory R script contains local paths and is not a clean pipeline. Therefore reproducibility begins at the precomputed edgeR contrast tables and must be described that way.

No AWS run is part of the local scientific evidence. Cloud execution becomes reproducible only after an approved deployment records the container digest, Terraform plan, region, data object versions, IAM boundary, and execution logs.
