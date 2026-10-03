# TB undernutrition transcriptomics

Historical v1.0.0 Zenodo DOI:
[10.5281/zenodo.22685970](https://doi.org/10.5281/zenodo.22685970)

Reproducibility materials for the study:

**Distinct transcriptional remodeling of immune and metabolic pathways in severe undernutrition during latent tuberculosis infection does not consistently predict tuberculosis progression across independent cohorts**

This repository contains the canonical analysis scripts, frozen seven-module programme, exact R-confirmed derived results, five main figures, Supplementary Methods and Tables S1-S6, Supplementary Data S1, and the v1.1.0 GC6 overlap-robustness addendum.

Public source datasets: GSE152218, GSE79362, and GSE94438/GC6.

Large upstream source matrices and GENCODE/Reactome resources are not redistributed. Their exact locally used copies are fingerprinted in `provenance/UPSTREAM_RESOURCE_MANIFEST.tsv`.

Canonical confirmation:

`./run_phaseA_confirmation.sh`

`./run_phaseB_R_refresh.sh`

Supplementary Data S1 internal checksum verification: 47/47 PASS.

Release lineage:

- v1.0.0: original public reproducibility release, archived at DOI 10.5281/zenodo.22685970.
- v1.1.0: publicly released on 2026-10-03, adding GC6 gene-identity robustness analyses and reconciling release metadata to the 14-author manuscript-v28 author list. Version DOI: 10.5281/zenodo.23119331.
- All versions: Zenodo concept DOI 10.5281/zenodo.22685969.
