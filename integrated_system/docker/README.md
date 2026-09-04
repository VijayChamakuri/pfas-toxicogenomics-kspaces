# Container workflows

The API image is the production artifact. It runs as UID 10001, has a health check, drops Linux capabilities in Compose, and sees the research project through a read-only mount.

```bash
cd pfas-toxicogenomics-kspaces
docker build --target api -f integrated_system/docker/Dockerfile -t pfas-workflow-api:local .
docker run --rm -p 8000:8000 --read-only --tmpfs /tmp \
  --cap-drop ALL --security-opt no-new-privileges \
  -v "$PWD/..:/data:ro" pfas-workflow-api:local
curl --fail http://localhost:8000/health
```

No credentials are baked into either image. Pass local secrets with `--env-file` from an ignored file, or mount a short-lived credential file read-only. In AWS, use the ECS task role and Secrets Manager injection instead of environment files.

## Scientific notebook and modeling image

Direct notebook dependencies are exactly pinned. Before publishing the image, generate and review a transitive hash lock on the target architecture. Local builds use the complete repository root so the scientific analysis, verified data inputs, and planning system share one filesystem contract:

```bash
docker build --target notebook -f integrated_system/docker/Dockerfile -t pfas-workflow-notebook:local .
docker run --rm pfas-workflow-notebook:local \
  python -c "import torch, tensorflow as tf; print(torch.__version__, tf.__version__)"
docker compose -f integrated_system/docker/compose.yaml --profile learning up notebook
```

For a publication image, run `pip-compile --generate-hashes` from `docker/requirements-notebook.in` and change the install to `--require-hashes`. Framework imports are checked by the manual CI job because this image is large.
