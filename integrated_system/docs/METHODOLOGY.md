# Methods and system architecture

## Scientific analysis boundary

The project starts with precomputed edgeR contrasts. It validates their common gene universe, reconstructs the original k-spaces workflow for audit, and applies corrected preprocessing, consensus fitting, stability checks, null comparisons, and qualified interpretation. The canonical scientific implementation and executable narrative live in the adjacent `K-spaces Analysis` directory.

The recommended module assignment is the 30-fit majority-voted consensus with per-gene `consensus_strength`. A single expectation-maximization fit is retained only as a reconstruction artifact because local optima produce meaningful run-to-run variation. The original correlation of each chemical's percent-per-module vector is retracted because a null stress test showed that high correlations arise largely from the five-part composition.

## Analysis-planning contract

A research request passes through five explicit boundaries:

1. Pydantic request models normalize the question and requested scope.
2. Compatibility checks confirm that chemicals, files, fields, and requested conclusions exist and are supportable.
3. Retrieval selects versioned evidence describing the dataset, corrected methods, assumptions, known failures, and interpretation rules.
4. The planner emits an allowlisted structured plan with evidence identifiers, parameters, limitations, and a human-review flag.
5. The executor revalidates the plan and input hashes before calling a bounded artifact operation.

The service abstains when metadata is missing, a requested operation is incompatible, sample size is insufficient, an identifier is unknown, a causal claim is requested, or the plan cannot be validated safely.

## Evidence grounding

The retrieval corpus is file-backed because the verified collection is small and versioned. Exact and deterministic retrieval makes evidence selection inspectable and reproducible. Important planning decisions identify the retrieved evidence. A vector database is not required at the current corpus size and would not replace corpus curation, citation validation, or scientific review.

## Model adaptation

The adaptation path targets structured workflow generation, not biological prediction. Training examples encode supported operations, compatibility failures, required parameters, and abstention behavior. Data generation records provenance, removes duplicates, separates groups across train, validation, and held-out partitions, and reports training configuration and errors. The adapted model remains behind the same typed validation and execution boundary as every other planner.

## Cross-backend verification

PyTorch and TensorFlow implement the same compact classification contract over a synthetic workflow-request dataset. Shared definitions cover features, targets, group-disjoint partitions, preprocessing, capacity, seed, early stopping, metrics, checkpoints, and prediction output. Agreement and error differences help detect framework-specific implementation defects. This synthetic benchmark is an engineering control and is never used as PFAS biological evidence.

## Human review

Automated validation cannot resolve biological plausibility, interpretation scope, or publication readiness. Plans that cross those boundaries remain pending until a domain reviewer confirms the input provenance, statistical assumptions, correction family, uncertainty language, and requested claim.
