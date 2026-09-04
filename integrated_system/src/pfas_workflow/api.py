import os
from pathlib import Path

from fastapi import FastAPI, HTTPException

from .catalog import ArtifactCatalog
from .executor import WorkflowExecutor
from .modeling import ModelingContract
from .models import ApprovedExecution, BiologicalRequest, ExecutionResult, WorkflowPlan
from .planner import RuleBasedPlanner

app = FastAPI(title="PFAS Evidence-Grounded Workflow Service", version="0.2.0")
project_root = Path(os.environ.get("PFAS_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
catalog = ArtifactCatalog(project_root)
planner = RuleBasedPlanner(catalog)
executor = WorkflowExecutor(catalog)


@app.get("/health")
def health():
    return {"status": "ok", "dataset": catalog.validate()}


@app.get("/modeling-contract", response_model=ModelingContract)
def modeling_contract():
    return ModelingContract()


@app.post("/plan", response_model=WorkflowPlan)
def plan(request: BiologicalRequest):
    return planner.plan(request)


@app.post("/execute", response_model=ExecutionResult)
def execute(approved: ApprovedExecution):
    try:
        return executor.execute(approved)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
