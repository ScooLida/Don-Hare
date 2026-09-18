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
WAIT_FOR_VCF="${WAIT_FOR_VCF:-0}"
VCF_WAIT_INTERVAL="${VCF_WAIT_INTERVAL:-30}"
VCF_WAIT_TIMEOUT="${VCF_WAIT_TIMEOUT:-0}"

mkdir -p "$WORK_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1
cd "$WORK_DIR"

export BASE_PREFIX DATA_PREFIX ANALYSIS_DIR

echo "Test population pipeline started: $(date)"
echo "Working directory: $WORK_DIR"
echo "Outputs: ${DATA_PREFIX}_* and $ANALYSIS_DIR"

vcf_is_complete() {
    [ -s "$FINAL_VCF" ] || return 1
    [ -s "${FINAL_VCF}.tbi" ] || return 1
    [ "$(stat -c %Y "$FINAL_VCF")" -ge "$WAIT_STARTED" ] || return 1
    [ "$(stat -c %Y "${FINAL_VCF}.tbi")" -ge "$WAIT_STARTED" ] || return 1
    bcftools index -n "$FINAL_VCF" >/dev/null
}

wait_for_vcf() {
    local elapsed=0
    WAIT_STARTED=$(date +%s)
    echo "Waiting for an externally running scr_11_bam_to_vcf_test.sh"
    echo "Expected VCF: $FINAL_VCF"
    while true; do
        if vcf_is_complete; then
            echo "Test VCF stage completed: $(bcftools index -n "$FINAL_VCF") sites"
            return 0
        fi
        if [ "$VCF_WAIT_TIMEOUT" -gt 0 ] && [ "$elapsed" -ge "$VCF_WAIT_TIMEOUT" ]; then
            echo "Error: timed out waiting for test VCF: $FINAL_VCF" >&2
            return 1
        fi
        sleep "$VCF_WAIT_INTERVAL"
        elapsed=$((elapsed + VCF_WAIT_INTERVAL))
    done
}

if [ "$WAIT_FOR_VCF" = "1" ]; then
    command -v bcftools >/dev/null 2>&1 || {
        echo "Error: bcftools is required in WAIT_FOR_VCF mode." >&2
        exit 1
    }
    wait_for_vcf
else
    echo "[1/4] Calling per-sample VCFs and creating ${DATA_PREFIX}_with_all_samples.vcf.gz"
    bash "$SCRIPT_DIR/scr_11_bam_to_vcf_test.sh"
fi

if [ "$WAIT_FOR_VCF" != "1" ] && { [ ! -s "$FINAL_VCF" ] || [ ! -s "${FINAL_VCF}.tbi" ]; }; then
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
