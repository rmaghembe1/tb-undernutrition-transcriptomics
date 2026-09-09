# Source data and external resources

The transcriptomic source datasets analyzed by this project are public:

- GSE152218 — discovery and internal active-TB convergence analyses
- GSE79362 — prospective progression projection
- GSE94438 / GC6 — assay-portability-gated progression analyses

Large source matrices and upstream annotation/pathway resources are intentionally
not redistributed in this repository.

Exact locally used upstream resources are fingerprinted in:

`provenance/UPSTREAM_RESOURCE_MANIFEST.tsv`

The analysis used GENCODE v49 for identifier harmonization and Reactome gene
sets for ranked pathway analysis.

External projection matrices were obtained using curatedTBData. Export scripts
and corresponding environment/session provenance are retained in this release.

Expected source-input locations used by the analysis scripts are preserved under
the `data/` directory structure.
