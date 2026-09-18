#!/usr/bin/env bash
# Run the complete _t population test pipeline except Fbranch.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$HOME/hare_work}"
BASE_PREFIX="${BASE_PREFIX:-MyHare}"
DATA_PREFIX="${BASE_PREFIX}_t"
ANALYSIS_DIR="${ANALYSIS_DIR:-./population_analysis_t}"
LOG_FILE="${LOG_FILE:-$WORK_DIR/test_population_pipeline_t.log}"
FINAL_VCF="$WORK_DIR/${DATA_PREFIX}_with_all_samples.vcf.gz"

mkdir -p "$WORK_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1
cd "$WORK_DIR"

export BASE_PREFIX DATA_PREFIX ANALYSIS_DIR

echo "Test population pipeline started: $(date)"
echo "Working directory: $WORK_DIR"
echo "Outputs: ${DATA_PREFIX}_* and $ANALYSIS_DIR"

echo "[1/4] Calling per-sample VCFs and creating ${DATA_PREFIX}_with_all_samples.vcf.gz"
bash "$SCRIPT_DIR/scr_11_bam_to_vcf_test.sh"

if [ ! -s "$FINAL_VCF" ] || [ ! -s "${FINAL_VCF}.tbi" ]; then
    echo "Error: test VCF stage did not produce a complete all-sample VCF: $FINAL_VCF" >&2
    exit 1
fi
if command -v bcftools >/dev/null 2>&1; then
    echo "Completed test VCF sites: $(bcftools index -n "$FINAL_VCF")"
fi

echo "[2/4] Running PCA and ADMIXTURE"
bash "$SCRIPT_DIR/scr_12_pca_admixture.sh"

echo "[3/4] Building PCA and ADMIXTURE plots"
Rscript "$SCRIPT_DIR/scr_13_population_plots.R"

echo "[4/4] Running Dsuite Dtrios"
bash "$SCRIPT_DIR/scr_14_dsuite.sh" "$FINAL_VCF"

echo "Test population pipeline finished: $(date)"
echo "Fbranch was not run. Use scr_15_fbranch_new_trees.sh separately with TREE_SOURCE."
