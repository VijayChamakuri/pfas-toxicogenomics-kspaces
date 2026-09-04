# Data use and governance

## Data inventory

The analysis consumes ten precomputed PFAS-versus-control differential-expression tables from `../Data/DEGs/`. Each table contains the same filtered universe of 13,852 genes. The k-spaces analysis uses the 10,255-gene union detected as differentially expressed in at least one contrast under the declared threshold.

The repository also contains derived module assignments, enrichment tables, evaluation cases, synthetic workflow requests, model checkpoints, and execution reports. Synthetic workflow requests test software behavior and are kept separate from biological inputs and conclusions.

## Permitted claims

The available artifacts support descriptive differential-expression summaries, provisional consensus-module exploration, qualified functional enrichment, and exploratory chemical-pair comparisons. They do not support independent reconstruction of sequencing QC or the edgeR design, causal toxicological conclusions, or population-level risk estimates.

## Missing provenance

Raw RNA-seq counts, complete sample metadata, sequencing QC records, and the original reproducible R environment are not present. GenX uses a different vehicle or control, which confounds direct cross-chemical interpretation. These limits must travel with exported results.

## Sharing boundary

This curated public repository was assembled after project-owner confirmation of principal-investigator and collaborator approval. The broader private research directory is not part of the release because it contains individual-level screening data, local session artifacts, collaborator documents, and materials with unclear redistribution rights. Future releases must:

1. identify the owner and license for every data and reference artifact;
2. remove local paths, session files, caches, credentials, and personal material;
3. document each input's source, version, access date, and permitted use;
4. confirm whether derived tables and figures may be redistributed;
5. attach the scientific limitations and provenance manifest to released outputs.

## Integrity controls

Inputs are treated as immutable. The catalog validates expected chemicals, schemas, row counts, identifiers, and hashes before execution. Generated results should record source paths, SHA-256 hashes, code revision, environment, parameters, seed, and human-review status. Never replace a source artifact in place while retaining its prior identifier or hash.

## Model-data separation

The cross-framework validation and workflow-planning adaptation data do not contain experimental observations and must not be merged into scientific tables. Production requests must not be reused for training without consent, governance approval, redaction, and a documented retention policy.
