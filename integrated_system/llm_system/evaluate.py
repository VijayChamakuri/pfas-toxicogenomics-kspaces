from __future__ import annotations

import argparse
import json
from pathlib import Path

from .rag import GroundedWorkflowSystem, PromptOnlySystem


def evaluate(cases_path: Path) -> dict:
    cases = json.loads(cases_path.read_text(encoding="utf-8"))
    systems: dict[str, dict] = {}
    report = {"schema_version": "1.0", "case_count": len(cases), "systems": systems}
    for name, system in (("prompt_only", PromptOnlySystem()), ("rag", GroundedWorkflowSystem())):
        rows = []
        totals = {"status": 0, "intent": 0, "unsupported_claims": 0, "workflow_validity": 0}
        accepted_count = 0
        citation_coverage = 0
        citation_correctness = 0
        for case in cases:
            result = system.propose(case["request"])
            if result.status == "accepted":
                accepted_count += 1
                citation_coverage += int(bool(result.citations))
                citation_correctness += int(bool(result.citations) and all(
                    len(citation.sha256) == 64 and not citation.source.startswith("/")
                    for citation in result.citations
                ))
            checks = {
                "status": result.status == case["expected_status"],
                "intent": result.intent == case.get("expected_intent", "unknown"),
                "unsupported_claims": not any(term in result.explanation.lower() for term in ("proves that", "causes toxicity", "confirmed mechanism")),
                "workflow_validity": bool(result.workflow_steps) == (result.status == "accepted"),
            }
            for key, passed in checks.items():
                totals[key] += int(passed)
            rows.append({"id": case["id"], "checks": checks, "result": result.as_dict()})
        metrics = {key: round(value / len(cases), 4) for key, value in totals.items()}
        metrics["groundedness"] = round(citation_coverage / max(accepted_count, 1), 4)
        metrics["citation_correctness"] = round(citation_correctness / max(accepted_count, 1), 4)
        systems[name] = {
            "metrics": metrics,
            "cases": rows,
        }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    default = Path(__file__).with_name("eval_cases.json")
    parser.add_argument("--cases", type=Path, default=default)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = evaluate(args.cases)
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
