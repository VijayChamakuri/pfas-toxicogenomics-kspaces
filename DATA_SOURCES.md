# Data sources and third-party terms

## Differential-expression inputs

The ten files in `Data/DEGs/` are PFAS-versus-control edgeR result tables approved for this
public release by the project owner and collaborators. They share a 13,852-gene filtered universe.
They are derived statistical results, not raw sequencing reads or count matrices. This repository
therefore reproduces the cross-chemical analysis but cannot independently reconstruct sequencing
quality control, sample-level variance, or the original edgeR design.

Each file contains WormBase gene identifiers and the columns `logFC`, `logCPM`, `F`, `PValue`, and
`FDR`. Input SHA-256 hashes are validated inside the executed notebook.

## Gene Ontology snapshot

`K-spaces Analysis/data_external/go-basic.obo` is a vendored Gene Ontology basic ontology snapshot:

- data version: `releases/2026-06-15`
- SHA-256: `c72fc198a86983d55e43aac585d1ffdbeb6e3601475b3f18b6045acdc0a0734c`
- upstream: <https://geneontology.org/>
- license: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

Gene Ontology data are copyright the Gene Ontology Consortium and are redistributed under CC BY
4.0. The repository MIT license applies to original project software, not this third-party data.

## WormBase annotations distributed through GO

`K-spaces Analysis/data_external/wb.gaf.gz` is the compressed *C. elegans* Gene Association File:

- GAF generation date: `2026-05-21T09:09`
- SHA-256: `bda92f6642cc68a16e485ba32eb5890320855895185b9cbcd448f52cf698fed2`
- upstream annotation download documentation: <https://geneontology.org/docs/download-go-annotations/>
- format: GAF 2.2

Run `python "K-spaces Analysis/scripts/prepare_annotations.py"` to validate and decompress this
snapshot before rebuilding the notebook. Annotation provenance and limitations remain visible in
the notebook and `K-spaces Analysis/FINDINGS.md`.

## Method software

The analysis pins `kspaces==0.1.5`, the Python implementation from the
[Pachter Lab k-spaces repository](https://github.com/pachterlab/k-spaces), which is distributed
under the BSD 2-Clause license. Cite Markarian et al., “k-spaces: Mixtures of Gaussian latent
variable models,” bioRxiv 2025.11.24.690254.

