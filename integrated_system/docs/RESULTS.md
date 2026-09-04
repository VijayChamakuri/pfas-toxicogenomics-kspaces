# Results narrative

## Data integrity and analytical basis

The ten canonical PFAS-versus-control contrast tables share a validated 13,852-gene filtered
universe. Biological interpretation begins from these precomputed edgeR results because raw
counts, complete sample metadata, and sequencing quality-control records are unavailable. The
module analysis preserves continuous log-fold changes rather than replacing effects with zero
when a contrast does not meet a detection threshold.

## Module analysis and uncertainty

The recommended module assignment is the 30-fit, majority-voted k=5 consensus with a
per-gene consensus-strength score. The single-fit reconstruction remains useful for explaining
the initial workflow, but its expectation-maximization local optimum is not the recommended
scientific result. Optimizer variation, observation resampling, chemical jackknife checks, and
null comparisons are treated as distinct uncertainty questions.

The original chemical-similarity claim based on correlations among five-part module-usage
profiles is retracted. A null stress test showed that high correlations arise readily from the
composition itself. Residual comparisons now adjust empirical p-values across all 45 pairs
separately at k=5 and k=10. No pair meets BH q below 0.05 at either k, and the earlier
eleven-pair unadjusted screen is withdrawn. Full clustering uncertainty is still not propagated.

## Functional interpretation

Functional summaries use locally pinned ontology and annotation inputs and report identifier
mapping loss. Consensus-module enrichment uses clustered DEG-union genes mapped in each namespace.
Chemical-exclusive and chemical-shared studies use the measured mapped universe, while each
chemical-by-module study uses mapped genes within its module. Each declared analysis family has
one Benjamini-Hochberg correction across all terms it tests. Terms remain exploratory associations
because annotation and module-assignment uncertainty are not fully propagated and they do not
establish causality.

## Analysis-planning quality control

The production planner accepts supported requests only after artifact compatibility checks,
evidence retrieval, typed-plan validation, and scientific-rule checks. It returns clarification,
incompatibility, or abstention for missing metadata, unknown entities, unavailable columns,
unsupported tools, causal claims, and retracted methods. The evaluation adapter now invokes
that same production planner, so citation, grounding, hallucination, and abstention results
exercise the API decision path rather than a parallel implementation.

PyTorch and TensorFlow implement the same compact classifier over a separate synthetic
workflow-request dataset. Their results test backend portability and implementation consistency
only. They do not contribute to PFAS biological findings. The adapted language model remains a
documented routing candidate behind deterministic validation and is not an approved execution
gate.

The executed master notebook is the detailed source for numerical tables and figures. Readers
should interpret every output together with the status and uncertainty text adjacent to it.
