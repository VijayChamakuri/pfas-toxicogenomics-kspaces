# Troubleshooting

## Python version or installation failure

Use Python 3.11 or later for this package. The adjacent k-spaces environment requires Python 3.10 or later and is verified with Python 3.13 on the development machine.

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/pip install -e '.[dev]'
```

Install `ml` or `finetune` extras only for the corresponding runs. TensorFlow and PyTorch downloads are large and may require architecture-specific wheels.

## Project root or missing artifacts

If planning succeeds but execution cannot find data, set `PFAS_PROJECT_ROOT` to the repository root:

```bash
export PFAS_PROJECT_ROOT="$PWD"
```

Do not point it at `integrated_system` or directly at `Data`. Verify that all ten `Data/DEGs/*vsControl_DGE_results.csv` files exist and remain unchanged.

## Import errors

For source-checkout commands, either install the package in editable mode or set the documented path explicitly:

```bash
.venv/bin/pip install -e '.[dev]'
PYTHONPATH=src .venv/bin/python scripts/run_evals.py
```

The framework and adaptation packages live at the repository root and are included by the editable install.

## TensorFlow or PyTorch backend failure

Run one backend at a time to isolate installation or memory problems:

```bash
PYTHONPATH=. .venv/bin/python -m ml_frameworks --smoke --framework pytorch \
  --output artifacts/pytorch_smoke
PYTHONPATH=. .venv/bin/python -m ml_frameworks --smoke --framework tensorflow \
  --output artifacts/tensorflow_smoke
```

Compare package versions, device selection, seed, split manifest, and preprocessing before comparing predictions. A backend-specific pass does not replace the cross-backend agreement check.

## Notebook execution failure

Regenerate notebooks from their builder scripts rather than editing generated `.ipynb` files. Use the registered `kspaces_venv` kernel and run from `K-spaces Analysis`. Full k-spaces execution can take substantially longer than fixture-based CI.

## Docker unavailable

Confirm both the Docker client and daemon are available:

```bash
docker version
cd pfas-toxicogenomics-kspaces
docker build --target api -f integrated_system/docker/Dockerfile -t pfas-workflow-api:local .
```

If `docker version` cannot connect to the daemon, start Docker Desktop or the approved container runtime. Do not report a successful image build until the build and health check both complete.

## Terraform validation failure

Run initialization without a backend before local validation:

```bash
cd infra/terraform
terraform init -backend=false
terraform fmt -check -recursive
terraform validate
```

Planning also requires reviewed networking and load-balancer identifiers in an ignored `terraform.tfvars`. Do not run `terraform apply` without explicit approval and cost review.
