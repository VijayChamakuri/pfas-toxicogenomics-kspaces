# Evaluation and results

## Evaluation layers

The project treats evaluation as a scientific quality-control system rather than a single model score.

| Layer | Primary checks | Interpretation |
|---|---|---|
| Data integrity | Required files, schemas, identifiers, hashes, shared universe | Confirms inputs match the declared contract |
| Scientific reconstruction | Legacy output coverage, corrected metrics, null checks, consensus stability | Distinguishes reproduced behavior from supported conclusions |
| Workflow planning | Exact plan fields, compatible operation, required parameters | Confirms that requests become executable plans |
| Retrieval grounding | Evidence relevance, citation presence and validity | Confirms that important decisions can be traced |
| Safety | Unsupported claims, hallucinated identifiers, missing metadata, abstention | Confirms that invalid requests do not reach execution |
| Modeling | Held-out metrics, backend agreement, error slices, reproducibility | Detects implementation-specific behavior |
| Execution | Allowlist enforcement, artifact lineage, deterministic output | Confirms that approved plans run through bounded code paths |

## Adversarial cases

Versioned cases cover correct, ambiguous, incomplete, dataset-incompatible, statistically invalid, causal, hallucinated-column, hallucinated-gene, unsupported-tool, incorrect-citation, overconfident, and retracted-method requests. Each case declares the expected operation or refusal reason. Scientific judgment cases are evaluated with `llm_system/HUMAN_REVIEW_RUBRIC.md`, not an automated model judge.

## Reported evidence

Machine-readable reports live beside the system that produces them:

- `llm_system/evaluation_report.json` records workflow-planning and grounding results.
- `llm_system/training_output/metrics.json` records the model-adaptation run.
- `artifacts/framework_parity_smoke/parity_report.json` records PyTorch and TensorFlow smoke results.
- The master notebook records scientific tables, figures, execution outputs, and limitations.

These artifacts are snapshots, not timeless claims. Compare their dataset, corpus, model, environment, and code identifiers before comparing results.

## Interpretation rules

- Passing synthetic workflow classification does not validate PFAS biology.
- Lower adaptation loss does not establish valid structured generation. Held-out schema validity, routing accuracy, abstention, and error categories are also required.
- Endpoint health does not establish scientific correctness.
- Retrieval relevance does not establish that a citation supports the generated claim. Citation correctness is checked separately.
- Exploratory biological signals remain exploratory when multiplicity or clustering uncertainty is unresolved.

## Reproducing evaluation

```bash
cd pfas-toxicogenomics-kspaces/integrated_system
PYTHONPATH=src .venv/bin/python scripts/run_evals.py
.venv/bin/python -m pytest -q tests/test_workflow.py tests/test_llm_rag.py
PYTHONPATH=. .venv/bin/python -m ml_frameworks --smoke --framework both \
  --output artifacts/framework_parity_smoke
```

Review generated JSON artifacts and the human-review rubric together. A release candidate should not proceed when any required case regresses, citations cannot be resolved, or an unsupported request reaches execution.
