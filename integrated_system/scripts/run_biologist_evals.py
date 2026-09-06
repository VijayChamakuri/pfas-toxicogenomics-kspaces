"""Grade the rule-based planner against the biologist eval set.

Reports overall and per-category accuracy, safety-critical refusal behaviour,
and every divergence between expected and observed planner outcomes.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from pfas_workflow.models import BiologicalRequest
from pfas_workflow.planner import RuleBasedPlanner

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "evals" / "biologist_cases.json"
HOLDOUT = ROOT / "evals" / "holdout_cases.json"

# Categories where a wrong answer is a scientific-integrity failure, not a UX miss.
SAFETY_CATEGORIES = {"abstain_causal", "abstain_retracted"}


def grade(case: dict, plan) -> dict:
    observed_status = plan.status.value
    observed_analysis = plan.analysis.value if plan.analysis else None
    observed_codes = {finding.code for finding in plan.findings}

    status_ok = observed_status == case["expected_status"]
    analysis_ok = observed_analysis == case["expected_analysis"]
    codes_ok = set(case["expected_codes"]).issubset(observed_codes) if case["expected_codes"] else True
    warns_ok = set(case["expected_warning_codes"]).issubset(observed_codes) if case["expected_warning_codes"] else True
    expected_params = case.get("expected_parameters") or {}
    params_ok = all(plan.parameters.get(k) == v for k, v in expected_params.items())
    refusal_ok = (not plan.accepted) == (not case["accepted"])

    return {
        "id": case["id"],
        "category": case["category"],
        "question": case["question"],
        "expected_status": case["expected_status"],
        "observed_status": observed_status,
        "expected_analysis": case["expected_analysis"],
        "observed_analysis": observed_analysis,
        "expected_codes": case["expected_codes"],
        "expected_warning_codes": case["expected_warning_codes"],
        "observed_codes": sorted(observed_codes),
        "checks": {
            "status": status_ok,
            "analysis": analysis_ok,
            "codes": codes_ok,
            "warnings": warns_ok,
            "refusal": refusal_ok,
            "parameters": params_ok,
        },
        "expected_parameters": expected_params,
        "observed_parameters": plan.parameters,
        "passed": all([status_ok, analysis_ok, codes_ok, warns_ok, refusal_ok, params_ok]),
        "rationale": case["rationale"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--holdout", action="store_true", help="Grade the held-out set instead.")
    parser.add_argument("--fail-under", type=float, default=None,
                        help="Exit non-zero if overall accuracy falls below this fraction.")
    parser.add_argument("--fail-on-unsafe", action="store_true",
                        help="Exit non-zero if any safety-critical case is not refused.")
    args = parser.parse_args()

    cases = json.loads((HOLDOUT if args.holdout else CASES).read_text())
    planner = RuleBasedPlanner()
    results = [grade(case, planner.plan(BiologicalRequest(question=case["question"]))) for case in cases]

    by_category: dict[str, list[dict]] = defaultdict(list)
    for result in results:
        by_category[result["category"]].append(result)

    passed = sum(r["passed"] for r in results)
    check_totals = Counter()
    for result in results:
        for name, ok in result["checks"].items():
            check_totals[name] += int(ok)

    # Safety-critical: a causal or retracted request that was NOT refused.
    unsafe = [r for r in results
              if r["category"] in SAFETY_CATEGORIES and r["observed_status"] == "ready_for_review"]
    # False refusals: a valid scientific question the planner blocked.
    false_refusals = [r for r in results
                      if r["category"].startswith("valid_") and r["observed_status"] != "ready_for_review"]

    print(f"cases: {len(results)}   passed: {passed}   accuracy: {passed / len(results):.1%}\n")

    print("per-check accuracy")
    for name, count in check_totals.most_common():
        print(f"  {name:12} {count / len(results):6.1%}  ({count}/{len(results)})")

    print("\nper-category accuracy")
    for category in sorted(by_category):
        group = by_category[category]
        hits = sum(r["passed"] for r in group)
        flag = "  <-- SAFETY" if category in SAFETY_CATEGORIES and hits < len(group) else ""
        print(f"  {category:32} {hits:>3}/{len(group):<3} {hits / len(group):6.1%}{flag}")

    if unsafe:
        print(f"\nSAFETY FAILURES ({len(unsafe)}): accepted a request that must be refused")
        for result in unsafe:
            print(f"  [{result['id']}] {result['question']}")
            print(f"      -> {result['observed_status']} / {result['observed_analysis']}"
                  f"  (expected {result['expected_status']})")

    if false_refusals:
        print(f"\nFALSE REFUSALS ({len(false_refusals)}): blocked a valid scientific question")
        for result in false_refusals:
            print(f"  [{result['id']}] {result['question']}")
            print(f"      -> {result['observed_status']} {result['observed_codes']}")

    other = [r for r in results if not r["passed"] and r not in unsafe and r not in false_refusals]
    if other:
        print(f"\nOTHER DIVERGENCES ({len(other)})")
        for result in other:
            failed = [name for name, ok in result["checks"].items() if not ok]
            print(f"  [{result['id']}] {','.join(failed)}: "
                  f"expected {result['expected_status']}/{result['expected_analysis']} "
                  f"got {result['observed_status']}/{result['observed_analysis']}")
            missing = sorted(set(result["expected_codes"] + result["expected_warning_codes"])
                             - set(result["observed_codes"]))
            if missing:
                print(f"      missing codes: {missing}")
            if not result["checks"]["parameters"]:
                print(f"      parameters: expected {result['expected_parameters']} "
                      f"got {result['observed_parameters']}")

    summary = {
        "planner_version": "evidence-rules-2.0",
        "case_count": len(results),
        "accuracy": passed / len(results),
        "per_check": {name: count / len(results) for name, count in check_totals.items()},
        "per_category": {c: sum(r["passed"] for r in g) / len(g) for c, g in by_category.items()},
        "safety_failures": [r["id"] for r in unsafe],
        "false_refusals": [r["id"] for r in false_refusals],
        "results": results,
    }
    if args.json_out:
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.json_out}")

    if args.fail_on_unsafe and unsafe:
        return 1
    if args.fail_under is not None and passed / len(results) < args.fail_under:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
