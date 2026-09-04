import importlib.util
from pathlib import Path

import pytest

from ml_frameworks.dataset import generate_dataset, split_dataset
from ml_frameworks.pytorch_benchmark import run_pytorch
from ml_frameworks.tensorflow_benchmark import run_tensorflow


@pytest.mark.parametrize(
    ("module", "runner", "message"),
    [
        ("torch", run_pytorch, "PyTorch dependency blocker: install torch to run the PyTorch benchmark"),
        ("tensorflow", run_tensorflow, "TensorFlow dependency blocker: install tensorflow to run the TensorFlow benchmark"),
    ],
)
def test_framework_runs_or_reports_exact_dependency_blocker(tmp_path, module, runner, message):
    bundle = split_dataset(generate_dataset(examples_per_group=2))
    if importlib.util.find_spec(module) is None:
        with pytest.raises(RuntimeError, match=message):
            runner(bundle, tmp_path, epochs=1, patience=1)
    else:
        report = runner(bundle, tmp_path, epochs=2, patience=1)
        assert report["framework"] in {"pytorch", "tensorflow"}
        assert report["metrics"]["accuracy"] >= 0.0
        assert report["parameter_count"] > 0
        assert Path(report["checkpoint"]).exists()
