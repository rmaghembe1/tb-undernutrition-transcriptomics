#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${1:-$SCRIPT_DIR}"
SCRIPT_DIR="$PROJECT_ROOT/scripts/R"

echo "Project root: $PROJECT_ROOT"

for required in \
  "$PROJECT_ROOT/data/raw/GSE152218_features_combined.txt.gz" \
  "$PROJECT_ROOT/data/metadata/metadata_gse152218_adult_ltbi_primary.tsv" \
  "$PROJECT_ROOT/data/pathways/ReactomePathways.gmt.zip" \
  "$PROJECT_ROOT/data/mapping/gencode.v49.ENSG_to_symbol.tsv.gz"; do
  if [[ ! -s "$required" ]]; then
    echo "ERROR: Missing or empty required input: $required" >&2
    exit 1
  fi
done

mkdir -p "$PROJECT_ROOT/results/R_confirmation"

echo "Running exact limma-voom discovery confirmation..."
Rscript "$SCRIPT_DIR/01_confirm_GSE152218_limma_voom.R" "$PROJECT_ROOT"

echo "Running exact Reactome fgsea confirmation..."
Rscript "$SCRIPT_DIR/02_confirm_GSE152218_Reactome_fgsea.R" "$PROJECT_ROOT"

echo "Confirmation outputs:"
find "$PROJECT_ROOT/results/R_confirmation" -maxdepth 2 -type f -printf "%p\n" | sort

echo "DONE: Upload the result folders or zip them for review."
