"""PyTorch implementation of the parity benchmark."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np

from .common import (
    BenchmarkTimer,
    build_vocabulary,
    classification_metrics,
    error_analysis,
    vectorize,
)
from .dataset import LABELS, DatasetBundle


def run_pytorch(bundle: DatasetBundle, output_dir: Path, seed: int = 2026, epochs: int = 40, patience: int = 5) -> dict[str, object]:
    try:
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise RuntimeError("PyTorch dependency blocker: install torch to run the PyTorch benchmark") from exc

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    vocabulary = build_vocabulary(bundle.train)
    x_train, y_train = vectorize(bundle.train, vocabulary)
    x_val, y_val = vectorize(bundle.validation, vocabulary)
    x_test, y_test = vectorize(bundle.test, vocabulary)
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(TensorDataset(torch.from_numpy(x_train), torch.from_numpy(y_train)), batch_size=16, shuffle=True, generator=generator)

    class Classifier(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.layers = nn.Sequential(nn.Linear(len(vocabulary), 24), nn.ReLU(), nn.Linear(24, len(LABELS)))

        def forward(self, values: torch.Tensor) -> torch.Tensor:
            return self.layers(values)

    model = Classifier()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.02)
    loss_fn = nn.CrossEntropyLoss()
    best_loss = float("inf")
    best_state = None
    stale = 0
    history = []
    with BenchmarkTimer() as timer:
        for epoch in range(epochs):
            model.train()
            for batch_x, batch_y in loader:
                optimizer.zero_grad()
                loss = loss_fn(model(batch_x), batch_y)
                loss.backward()
                optimizer.step()
            model.eval()
            with torch.no_grad():
                val_loss = float(loss_fn(model(torch.from_numpy(x_val)), torch.from_numpy(y_val)))
            history.append({"epoch": epoch + 1, "validation_loss": val_loss})
            if val_loss < best_loss - 1e-6:
                best_loss, stale = val_loss, 0
                best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
            else:
                stale += 1
                if stale >= patience:
                    break
    if best_state is None:
        raise RuntimeError("PyTorch training did not produce a checkpoint")
    model.load_state_dict(best_state)
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = output_dir / "pytorch_checkpoint.pt"
    torch.save({"state_dict": best_state, "vocabulary": vocabulary, "labels": LABELS}, checkpoint)
    model.eval()
    with torch.no_grad():
        predicted = model(torch.from_numpy(x_test)).argmax(dim=1).numpy()
    metrics = classification_metrics(y_test, predicted)
    return {
        "framework": "pytorch", "framework_version": torch.__version__, "seed": seed,
        "parameter_count": sum(p.numel() for p in model.parameters()), "epochs_completed": len(history),
        "best_validation_loss": best_loss, "runtime_seconds": timer.seconds,
        "python_peak_memory_bytes": timer.peak_bytes, "checkpoint": str(checkpoint),
        "metrics": metrics, "errors": error_analysis(bundle.test, y_test, predicted), "history": history,
    }
