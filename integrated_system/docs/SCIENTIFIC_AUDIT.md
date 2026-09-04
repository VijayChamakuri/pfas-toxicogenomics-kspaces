# Scientific audit

## Confirmed high-priority issues

1. RNA-seq results are internally consistent but not independently regenerable from raw counts. GenX also uses a different vehicle/control, so its contrast is confounded with vehicle when chemicals are compared.
2. Residual pair tests now control BH FDR across all 45 pairs separately at k=5 and k=10. No pair meets q below 0.05 at either k, so the earlier eleven-pair unadjusted screen is withdrawn. Conditional nulls still hold the consensus clustering fixed and do not propagate model uncertainty.
3. Consensus-module GO enrichment uses a namespace-specific eligible population formed by intersecting the clustered DEG union with mapped genes. Chemical-exclusive and chemical-shared studies use the measured mapped universe, while chemical-by-module studies use the mapped genes within each module. Each declared family has one Benjamini-Hochberg correction. Annotation and module uncertainty remain outside those calculations.
4. The 30-permutation likelihood test yields a plus-one Monte Carlo p-value of 0.0323. It supports cross-chemical dependence under that shuffle null, not a unique or biologically valid module structure.

## Corrections made

`FINDINGS.md` and the notebook builder withdraw the unadjusted eleven-pair screen, report the multiplicity-controlled null result, and disclose missing clustering-uncertainty propagation. The executed master notebook records the corrected continuous-effect primary matrix and declared family-wide GO procedures.

## Publication and sharing gate

Only the curated files in this public repository are approved for redistribution. Do not copy material from the broader private research workspace without a new ownership, privacy, and data-use review.
