#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$SCRIPT_DIR}"

echo "Phase B project root: $ROOT"

REQUIRED=(
  "$ROOT/data/raw/GSE152218_features_combined.txt.gz"
  "$ROOT/data/mapping/gencode.v49.ENSG_to_symbol.tsv.gz"
  "$ROOT/data/pathways/ReactomePathways.gmt.zip"
  "$ROOT/data/metadata/gse152218_severe_undernutrition_TB_axis_metadata.tsv"
  "$ROOT/data/modules/exact_confirmed_seven_module_programme.tsv"
  "$ROOT/results/R_confirmation/limma_voom/limma_voom_sensitivity_excluding_102-00459-B.tsv"
  "$ROOT/results/R_confirmation/reactome_fgsea/reactome_fgsea_sensitivity_excluding_102-00459-B.tsv"
  "$ROOT/results/R_confirmation/reactome_fgsea/reactome_fgsea_pathway_robustness.tsv"
  "$ROOT/data/external_validation/GSE79362/curatedTBData/GSE79362_assay_curated.tsv.gz"
  "$ROOT/data/external_validation/GSE79362/curatedTBData/GSE79362_colData.tsv"
  "$ROOT/data/external_validation/GSE94438/curatedTBData/GSE94438_assay_curated.tsv.gz"
  "$ROOT/data/external_validation/GSE94438/curatedTBData/GSE94438_colData.tsv"
)

for f in "${REQUIRED[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: Missing or empty required input: $f" >&2
    exit 1
  fi
done

mkdir -p "$ROOT/results/R_confirmation_phaseB"

echo "1/3 Refreshing exact active-TB convergence gate..."
Rscript "$ROOT/scripts/R/03_refresh_activeTB_gate_exact_fgsea.R" "$ROOT"

echo "2/3 Confirming GSE79362 projection in R..."
Rscript "$ROOT/scripts/R/04_confirm_GSE79362_exact_modules.R" "$ROOT"

echo "3/3 Confirming GC6 portability-gated analysis in R..."
Rscript "$ROOT/scripts/R/05_confirm_GC6_portability_gated_analysis.R" "$ROOT"

cd "$ROOT"
rm -f TB_Undernutrition_PhaseB_R_Refresh_Results_v1.zip
zip -r TB_Undernutrition_PhaseB_R_Refresh_Results_v1.zip \
  results/R_confirmation_phaseB \
  data/modules/exact_confirmed_seven_module_programme.tsv \
  data/metadata/gse152218_severe_undernutrition_TB_axis_metadata.tsv \
  scripts/R/phaseB_helpers.R \
  scripts/R/03_refresh_activeTB_gate_exact_fgsea.R \
  scripts/R/04_confirm_GSE79362_exact_modules.R \
  scripts/R/05_confirm_GC6_portability_gated_analysis.R

echo
echo "DONE. Upload:"
echo "$ROOT/TB_Undernutrition_PhaseB_R_Refresh_Results_v1.zip"
