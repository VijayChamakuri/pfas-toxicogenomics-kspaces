# PFAS evidence-grounded analysis system

This package turns research questions about ten PFAS exposures in *C. elegans* into validated, traceable operations over verified toxicogenomics artifacts. It connects natural-language request interpretation, evidence retrieval, scientific guardrails, reproducible execution, cross-backend model verification, and deployment controls without expanding the biological claims supported by the available data.

## Research and system objective

The scientific analysis asks which transcriptional responses are shared across PFAS chemicals and which are chemical-specific. The software layer makes those analyses easier to request and safer to execute. A request is translated into a typed plan, checked against the available datasets and corrected methods, grounded in versioned evidence, reviewed when scientific judgment is required, and executed only through supported operations.

The package begins from ten precomputed edgeR contrast tables and consensus k-spaces outputs. It does not process raw reads or independently reconstruct differential expression because raw counts and complete sample metadata are not present.

## System flow

```text
research question
    -> request interpretation and compatibility checks
    -> retrieval from verified methods and dataset evidence
    -> typed workflow plan with citations and parameters
    -> scientific-rule validation and abstention when unsupported
    -> human-review checkpoint
    -> allowlisted artifact operation
    -> result with hashes, provenance, limitations, and review status
```

The modeling layer uses equivalent compact classifiers in PyTorch and TensorFlow to check backend portability and prediction consistency on a separate synthetic workflow-request dataset. These tests validate software behavior only. They do not establish biological prediction performance or neural-model superiority.

## Supported operations

| Operation | Evidence source | Scientific boundary |
|---|---|---|
| DEG summary | Ten PFAS-versus-control edgeR tables | FDR is computed within each contrast |
| Consensus module lookup | Majority-voted k=5 assignment and confidence | Modules remain provisional |
| Functional enrichment lookup | Existing consensus enrichment artifacts | Background definition and correction family require review |
| Exploratory chemical-pair lookup | Residual cross-k comparison table | No pair passes BH q<0.05 at k=5 or k=10 |

Requests for raw-read QC, unsupported chemicals, arbitrary code, causal conclusions, missing fields, or the retracted raw module-usage correlation are rejected with a reason.

## Repository structure

| Path | Role |
|---|---|
| `src/pfas_workflow/` | Typed request planning, artifact catalog, execution, and FastAPI boundary |
| `llm_system/` | Evidence corpus, retrieval, planning evaluation, adaptation pipeline, and review rubric |
| `ml_frameworks/` | Shared PyTorch and TensorFlow verification contract |
| `evals/` | Versioned adversarial workflow cases |
| `tests/` | Unit and integration tests for planning, retrieval, modeling, and execution |
| `docker/` | Reproducible container targets and run instructions |
| `infra/` | Approval-gated AWS deployment design and Terraform starter |
| `docs/` | Methods, data-use, evaluation, security, reproducibility, and troubleshooting records |

The merged scientific narrative is generated separately at `../K-spaces Analysis/PFAS_master_analysis.ipynb`. Reusable implementation remains in Python modules so the notebook is an executable report rather than the only source of logic.

## Quick start

Python 3.11 or later is required for this package. The k-spaces analysis uses its own Python 3.13 environment.

```bash
cd pfas-toxicogenomics-kspaces/integrated_system
python3.13 -m venv .venv
.venv/bin/pip install -e '.[dev,ml,finetune]'
.venv/bin/pytest -q
PYTHONPATH=src .venv/bin/python scripts/run_evals.py
PFAS_PROJECT_ROOT="$PWD/.." .venv/bin/uvicorn pfas_workflow.api:app --reload
```

For a lighter installation, omit `ml` and `finetune`. Those extras are needed only for backend-parity and model-adaptation runs. The service exposes `GET /health`, `POST /plan`, and `POST /execute`.

## Verification

```bash
.venv/bin/ruff check .
.venv/bin/python -m pytest -q
PYTHONPATH=src .venv/bin/python scripts/run_evals.py
PYTHONPATH=. .venv/bin/python -m ml_frameworks --smoke --framework both \
  --output artifacts/framework_parity_smoke
```

The evaluation set covers supported, ambiguous, incomplete, incompatible, causal, hallucinated, and retracted-method requests. Automated checks measure plan validity, retrieval grounding, citation integrity, abstention, and deterministic reproducibility. The human-review rubric governs questions that require scientific judgment.

## Documentation

- [Data use and governance](docs/DATA_USE.md)
- [Methods and system architecture](docs/METHODOLOGY.md)
- [Evaluation and results](docs/EVALUATION.md)
- [Scientific results narrative](docs/RESULTS.md)
- [Scientific workflow diagram](docs/scientific_workflow.mmd)
- [Security](docs/SECURITY.md)
- [Reproducibility](docs/BASELINE_AND_REPRODUCIBILITY.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [AWS deployment design](infra/README.md)
- [System and AWS architecture diagram](infra/architecture.mmd)
- [Scientific audit](docs/SCIENTIFIC_AUDIT.md)

## Scientific scope and validation status

The local analysis, service, evaluation, and infrastructure definitions are reproducible project artifacts. No AWS deployment is claimed, and no cloud resources should be created without explicit approval and cost review. The absence of raw RNA-seq counts and complete metadata prevents independent reconstruction of the original edgeR model. No chemical pair passes BH q below 0.05 at k=5 or k=10, and pathway associations remain exploratory because annotation and clustering uncertainty are not fully propagated.

Original project software is released under the repository MIT license. Vendored ontology and annotation data retain their upstream terms documented in `../DATA_SOURCES.md`.
