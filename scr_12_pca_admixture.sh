#!/bin/bash
# Run PCA on modern and all-sample datasets. Train ADMIXTURE on modern samples,
# then project all samples using the fixed modern allele-frequency model.

set -euo pipefail

# Configuration
DATA_PREFIX="MyHare"
MODERN_VCF="${DATA_PREFIX}_modern.vcf.gz"
ALL_VCF="${DATA_PREFIX}_with_all_samples.vcf.gz"
ANALYSIS_DIR="./population_analysis"
PLINK="$HOME/plink"
ADMIXTURE="$HOME/admixture/dist/admixture_linux-1.3.0/admixture"
PLINK_THREADS=8
ADMIXTURE_THREADS=4
MISSINGNESS=0.2
PCA_COMPONENTS=4
K_MIN=2
K_MAX=15
ANCIENT_SAMPLES="1k,3k,4k,5kS8"

prepare_dataset() {
    local label=$1
    local input_vcf=$2
    local prefix="$ANALYSIS_DIR/$label"

    if [ ! -f "$input_vcf" ]; then
        echo "Error: VCF not found for $label: $input_vcf" >&2
        return 1
    fi

    if [ "$label" = "modern" ]; then
        "$PLINK" --vcf "$input_vcf" --geno "$MISSINGNESS" --make-bed \
            --set-missing-var-ids '@:#_$1_$2' \
            --threads "$PLINK_THREADS" --out "$prefix" --allow-extra-chr
    else
        # The VCF already contains modern-selected sites. Do not filter on
        # missing ancient genotypes here.
        "$PLINK" --vcf "$input_vcf" --make-bed \
            --set-missing-var-ids '@:#_$1_$2' \
            --threads "$PLINK_THREADS" --out "$prefix" --allow-extra-chr
    fi

    "$PLINK" --bfile "$prefix" --pca "$PCA_COMPONENTS" \
        --threads "$PLINK_THREADS" --allow-extra-chr --out "$prefix"
}

run_modern_admixture() {
    local prefix="$ANALYSIS_DIR/modern"
    local cv_table="$ANALYSIS_DIR/modern_admixture_cv.tsv"

    printf "K\tCV_error\n" > "$cv_table"
    for K in $(seq "$K_MIN" "$K_MAX"); do
        log_file="$ANALYSIS_DIR/modern_admixture_K${K}.log"
        "$ADMIXTURE" -j"$ADMIXTURE_THREADS" --cv "$prefix.bed" "$K" \
            > "$log_file" 2>&1
        cv_error=$(awk '/CV error/ { value=$NF } END { print value }' "$log_file")
        if [ -n "$cv_error" ]; then
            printf "%s\t%s\n" "$K" "$cv_error" >> "$cv_table"
        fi
    done

    best_k=$(awk 'NR > 1 && $2 != "" { print }' "$cv_table" |
        sort -k2,2n | head -n 1 | cut -f1)
    if [ -z "$best_k" ]; then
        echo "Error: could not determine optimal K for modern samples." >&2
        return 1
    fi

    echo "$best_k" > "$ANALYSIS_DIR/modern_optimal_K.txt"
    echo "Optimal K for modern samples: $best_k"
}

project_all_samples() {
    local modern_prefix="$ANALYSIS_DIR/modern"
    local all_prefix="$ANALYSIS_DIR/with_all_samples"
    local best_k
    best_k=$(tr -d '[:space:]' < "$ANALYSIS_DIR/modern_optimal_K.txt")

    # Projection requires identical SNP IDs and order in both datasets.
    if ! cmp -s "${modern_prefix}.bim" "${all_prefix}.bim"; then
        echo "Error: modern and all-sample PLINK SNP sets differ." >&2
        return 1
    fi

    cp "${modern_prefix}.${best_k}.P" "${all_prefix}.${best_k}.P.in"
    cp "$ANALYSIS_DIR/modern_optimal_K.txt" \
        "$ANALYSIS_DIR/with_all_samples_optimal_K.txt"
    projection_log="$ANALYSIS_DIR/with_all_samples_projection_K${best_k}.log"
    "$ADMIXTURE" -j"$ADMIXTURE_THREADS" -P "$all_prefix.bed" "$best_k" \
        > "$projection_log" 2>&1

    projection_q="${all_prefix}.${best_k}.Q"
    if [ ! -s "$projection_q" ]; then
        echo "Error: projection Q file was not generated: $projection_q" >&2
        return 1
    fi

    ancient_table="$ANALYSIS_DIR/ancient_projection_K${best_k}.tsv"
    awk -v ancient="$ANCIENT_SAMPLES" -v k="$best_k" '
        BEGIN {
            split(ancient, names, ",")
            for (i in names) ancient_sample[names[i]] = 1
            printf "sample"
            for (i = 1; i <= k; i++) printf "\tQ%d", i
            print ""
        }
        NR == FNR {
            sample[FNR] = $2
            keep[FNR] = ($2 in ancient_sample)
            next
        }
        keep[FNR] {
            printf "%s", sample[FNR]
            for (i = 1; i <= NF; i++) printf "\t%s", $i
            print ""
        }
    ' "${all_prefix}.fam" "$projection_q" > "$ancient_table"

    echo "Ancient projection: $ancient_table"
    echo "Fixed modern model: ${modern_prefix}.${best_k}.P"
}

mkdir -p "$ANALYSIS_DIR"
prepare_dataset "modern" "$MODERN_VCF"
run_modern_admixture

prepare_dataset "with_all_samples" "$ALL_VCF"
project_all_samples

echo "PCA and fixed-model ADMIXTURE analyses complete."
