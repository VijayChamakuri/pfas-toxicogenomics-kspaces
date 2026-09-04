# Reviewed scientific and tool sources

This catalog identifies the local sources allowed to support workflow planning. It records what each source can support and what it cannot establish. Copyrighted papers remain in their original files and are not copied into the retrieval corpus verbatim.

| Source | Local record | SHA-256 or version | Permitted support | Boundary |
| --- | --- | --- | --- | --- |
| Markarian et al. k-spaces paper | `Background Reading/Markarian2025.pdf` | `2a813b8236271103d0bd61bceaa196d5ab194494c77288ceadbe636139b0ac78` | Definition and intended use of k-spaces; rationale for affine-subspace clustering | Does not validate this dataset, selected k, or biological interpretation |
| kspaces Python implementation | `K-spaces Analysis/requirements.txt` | `kspaces==0.1.5` | Executable implementation version used by the notebook | Package output still requires stability and null checks |
| Gene Ontology basic ontology | `K-spaces Analysis/data_external/go-basic.obo` | `c72fc198a86983d55e43aac585d1ffdbeb6e3601475b3f18b6045acdc0a0734c` | GO identifiers, namespaces, names, and ontology relationships | Ontology terms are annotations, not causal findings |
| WormBase GO annotation file | `K-spaces Analysis/data_external/wb.gaf.gz` | `bda92f6642cc68a16e485ba32eb5890320855895185b9cbcd448f52cf698fed2` | Mapping eligible WormBase identifiers to GO annotations | Unmapped genes and annotation date must be reported |
| GOATOOLS implementation | `K-spaces Analysis/requirements.txt` | Version pinned in the environment file | Enrichment computation and ontology parsing | Multiple-testing family and population remain project decisions |
| Corrected project findings | `K-spaces Analysis/FINDINGS.md` | Hash recorded when the corpus is built | Project decisions, corrected claims, and retractions | Summary defers to the executed notebook when they disagree |

## Retrieval policy

Planning retrieval uses short, project-authored review records plus their source hashes. It does not ingest entire copyrighted PDFs or treat retrieved text as new biological evidence. A hash change invalidates the existing corpus manifest until a reviewer accepts the new source version.
