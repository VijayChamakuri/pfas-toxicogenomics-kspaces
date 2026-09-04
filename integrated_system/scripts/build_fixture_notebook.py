"""Build a deterministic notebook smoke fixture for CI execution."""

from __future__ import annotations

import hashlib
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "tests" / "fixtures" / "workflow_smoke.ipynb"

cells = [
    nbf.v4.new_markdown_cell(
        "# Workflow lifecycle smoke check\n\n"
        "This fixture validates planning, review, execution, and provenance without full-data fitting."
    ),
    nbf.v4.new_code_cell(
        "from pathlib import Path\n"
        "import sys\n"
        "project_root = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'integrated_system').is_dir())\n"
        "sys.path.insert(0, str(project_root / 'integrated_system' / 'src'))\n"
        "from pfas_workflow.catalog import ArtifactCatalog\n"
        "from pfas_workflow.executor import WorkflowExecutor\n"
        "from pfas_workflow.models import ApprovedExecution, BiologicalRequest, HumanReview, ReviewDecision\n"
        "from pfas_workflow.planner import RuleBasedPlanner\n"
        "catalog = ArtifactCatalog(project_root)\n"
        "plan = RuleBasedPlanner(catalog).plan(BiologicalRequest(question='Summarize PFOS differential expression.'))\n"
        "assert plan.accepted and plan.citations\n"
        "review = HumanReview(decision=ReviewDecision.APPROVED, reviewer='ci scientific review', plan_digest=plan.review_digest())\n"
        "result = WorkflowExecutor(catalog).execute(ApprovedExecution(plan=plan, review=review))\n"
        "assert result.records[0]['chemical'] == 'PFOS'\n"
        "assert any('result_input' in item.roles for item in result.provenance.values())\n"
        "result.model_dump(mode='json')"
    ),
]
for index, cell in enumerate(cells):
    cell["id"] = hashlib.sha1(f"{index}:{cell.source}".encode()).hexdigest()[:12]

notebook = nbf.v4.new_notebook(cells=cells)
notebook.metadata.kernelspec = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3",
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(notebook, OUTPUT)
print(f"Wrote {OUTPUT}")
