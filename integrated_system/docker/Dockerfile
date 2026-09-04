# syntax=docker/dockerfile:1.7
ARG PYTHON_IMAGE=python:3.12.11-slim-bookworm

FROM ${PYTHON_IMAGE} AS base
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1
WORKDIR /app
RUN groupadd --system --gid 10001 app && \
    useradd --system --uid 10001 --gid app --create-home app

FROM base AS api
COPY integrated_system/pyproject.toml integrated_system/README.md ./
COPY integrated_system/src ./src
RUN python -m pip install --upgrade "pip==25.2" && \
    python -m pip install .
USER 10001:10001
ENV PFAS_PROJECT_ROOT=/data
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/health', timeout=3)" || exit 1
CMD ["uvicorn", "pfas_workflow.api:app", "--host", "0.0.0.0", "--port", "8000", "--no-access-log"]

FROM base AS notebook
COPY integrated_system/docker/requirements-notebook.txt /tmp/requirements-notebook.txt
COPY ["K-spaces Analysis/requirements.txt", "/tmp/requirements-science.txt"]
RUN python -m pip install --upgrade "pip==25.2" && \
    python -m pip install -r /tmp/requirements-science.txt -r /tmp/requirements-notebook.txt
COPY --chown=10001:10001 integrated_system /workspace/integrated_system
COPY --chown=10001:10001 Data/DEGs /workspace/Data/DEGs
COPY --chown=10001:10001 ["K-spaces Analysis", "/workspace/K-spaces Analysis"]
RUN python -m pip install /workspace/integrated_system
WORKDIR "/workspace/K-spaces Analysis"
USER 10001:10001
ENV JUPYTER_CONFIG_DIR=/tmp/jupyter \
    MPLCONFIGDIR=/tmp/matplotlib \
    PFAS_KSPACES_DIR="/workspace/K-spaces Analysis"
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--ServerApp.root_dir=/workspace"]
