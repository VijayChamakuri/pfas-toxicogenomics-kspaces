"""Command line entry point for dataset generation and parity runs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .common import write_report
from .dataset import generate_dataset, save_bundle, split_dataset
from .pytorch_benchmark import run_pytorch
from .tensorflow_benchmark import run_tensorflow


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the deterministic ML framework parity benchmark")
    parser.add_argument("--framework", choices=("pytorch", "tensorflow", "both"), default="both")
    parser.add_argument("--output", type=Path, default=Path("artifacts/framework_parity"))
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--smoke", action="store_true", help="Use two examples per group and at most three epochs")
    args = parser.parse_args(argv)
    bundle = split_dataset(generate_dataset(args.seed, 2 if args.smoke else 4), args.seed)
    save_bundle(bundle, args.output / "dataset")
    runners = {"pytorch": run_pytorch, "tensorflow": run_tensorflow}
    selected = tuple(runners) if args.framework == "both" else (args.framework,)
    report: dict[str, object] = {"dataset": bundle.manifest, "smoke": args.smoke, "results": {}, "blockers": {}}
    for name in selected:
        try:
            result = runners[name](bundle, args.output, args.seed, epochs=3 if args.smoke else 40, patience=2 if args.smoke else 5)
            report["results"][name] = result  # type: ignore[index]
        except RuntimeError as exc:
            report["blockers"][name] = str(exc)  # type: ignore[index]
    write_report(report, args.output / "parity_report.json")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 2 if report["blockers"] else 0


if __name__ == "__main__":
    sys.exit(main())
