from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class ModelingContract(BaseModel):
    contract_version: str = "1.0"
    objective: str = "classify workflow requests for independent implementation verification"
    input_features: str = "L2-normalized bag-of-token vectors fit on training records only"
    targets: tuple[str, ...] = ("accepted", "abstain", "incompatible", "needs_clarification")
    partitions: str = "group-disjoint train, validation, and test partitions"
    capacity: str = "single hidden layer with matched dimensions"
    metrics: tuple[str, ...] = ("accuracy", "macro_f1", "per_class", "confusion_matrix")
    seed: int = 2026
    checkpoint_behavior: str = "retain weights from the lowest validation loss"
    prediction_format: str = "class label, class probabilities, and framework identifier"
    uncertainty_outputs: str = "class probabilities and maximum-probability confidence"
    supported_backends: tuple[str, ...] = ("pytorch", "tensorflow")
    scientific_scope: str = "implementation verification only; outputs are not biological evidence"


class AdaptedRouterMetadata(BaseModel):
    component: str = "structured workflow routing candidate"
    model_name: str
    seed: int
    training_steps: int
    train_loss: float
    held_out_loss: float
    evaluation_status: str = Field(default="loss-only smoke evaluation; not approved for workflow decisions")


def load_router_metadata(metrics_path: Path) -> AdaptedRouterMetadata:
    metrics: dict[str, Any] = json.loads(metrics_path.read_text(encoding="utf-8"))
    return AdaptedRouterMetadata(
        model_name=str(metrics["model_name"]), seed=int(metrics["seed"]),
        training_steps=int(metrics["max_steps"]), train_loss=float(metrics["train_loss"]),
        held_out_loss=float(metrics["test_loss"]),
    )


def router_disclosure(metrics_path: Path | None = None) -> dict[str, Any]:
    disclosure: dict[str, Any] = {
        "selected_router": "deterministic evidence rules",
        "reason": "validated rules remain the execution gate",
    }
    if metrics_path and metrics_path.is_file():
        disclosure["adapted_model"] = load_router_metadata(metrics_path).model_dump()
    else:
        disclosure["adapted_model"] = {
            "evaluation_status": "training metrics unavailable; model is not used for workflow decisions"
        }
    return disclosure
