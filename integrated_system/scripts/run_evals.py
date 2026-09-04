import json
from pathlib import Path

from pfas_workflow.models import BiologicalRequest
from pfas_workflow.planner import RuleBasedPlanner

cases = json.loads((Path(__file__).parents[1] / "evals" / "cases.json").read_text())
correct = 0
refusal_tp = refusal_fp = refusal_fn = 0
for case in cases:
    plan = RuleBasedPlanner().plan(BiologicalRequest(question=case["question"]))
    predicted = plan.analysis.value if plan.analysis else None
    match = plan.accepted is case["accepted"] and predicted == case["analysis"]
    correct += int(match)
    expected_refusal, predicted_refusal = not case["accepted"], not plan.accepted
    refusal_tp += int(expected_refusal and predicted_refusal)
    refusal_fp += int(not expected_refusal and predicted_refusal)
    refusal_fn += int(expected_refusal and not predicted_refusal)
precision = refusal_tp / (refusal_tp + refusal_fp) if refusal_tp + refusal_fp else 0
recall = refusal_tp / (refusal_tp + refusal_fn) if refusal_tp + refusal_fn else 0
print(json.dumps({"planner_version":"rules-1.0","cases":len(cases),"exact_match":correct/len(cases),"refusal_precision":precision,"refusal_recall":recall}, indent=2))
