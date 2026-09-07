#!/bin/bash
# Continue scr_12 without recalculating modern PCA or ADMIXTURE.
# Run from hare_work after MyHare_with_all_samples.vcf.gz is ready.

set -euo pipefail

ANALYSIS_DIR="./population_analysis"
ALL_VCF="MyHare_with_all_samples.vcf.gz"
MODERN_PREFIX="$ANALYSIS_DIR/modern"
ALL_PREFIX="$ANALYSIS_DIR/with_all_samples"
ADMIXTURE="$HOME/admixture/dist/admixture_linux-1.3.0/admixture"
PLINK="$HOME/plink"
PLINK_THREADS=8
ADMIXTURE_THREADS=4
PCA_COMPONENTS=4
K_MIN=2
K_MAX=15
ANCIENT_SAMPLES="1k,3k,4k,5kS8"

if [ ! -s "$ALL_VCF" ] || [ ! -s "${ALL_VCF}.tbi" ]; then
    echo "Error: all-sample VCF is not ready: $ALL_VCF" >&2
    exit 1
fi
if [ ! -s "$MODERN_PREFIX.bed" ] || [ ! -s "$MODERN_PREFIX.bim" ] || [ ! -s "$MODERN_PREFIX.fam" ]; then
    echo "Error: modern PLINK files are missing. Run modern conversion first." >&2
    exit 1
fi

CV_TABLE="$ANALYSIS_DIR/modern_admixture_cv.tsv"
printf "K\tCV_error\n" > "$CV_TABLE"
for K in $(seq "$K_MIN" "$K_MAX"); do
    log_file="$ANALYSIS_DIR/modern_admixture_K${K}.log"
    if [ ! -s "$log_file" ]; then
        echo "Error: missing modern ADMIXTURE log: $log_file" >&2
        exit 1
    fi
    cv_error=$(awk '/CV error/ { value=$NF } END { print value }' "$log_file")
    if [ -z "$cv_error" ]; then
        echo "Error: no CV error found in $log_file" >&2
        echo "Inspect this log before continuing." >&2
        exit 1
    fi
    printf "%s\t%s\n" "$K" "$cv_error" >> "$CV_TABLE"
done

BEST_K=$(awk 'NR > 1 { print }' "$CV_TABLE" | sort -k2,2n | awk 'NR == 1 { print $1 }')
if [ -z "$BEST_K" ]; then
    echo "Error: could not determine the best K." >&2
    exit 1
fi
echo "Best modern K: $BEST_K"
printf '%s\n' "$BEST_K" > "$ANALYSIS_DIR/modern_optimal_K.txt"

if [ ! -s "$ALL_PREFIX.bed" ] || [ ! -s "$ALL_PREFIX.bim" ] || [ ! -s "$ALL_PREFIX.fam" ]; then
    "$PLINK" --vcf "$ALL_VCF" --make-bed \
        --set-missing-var-ids '@:#_$1_$2' \
        --threads "$PLINK_THREADS" --out "$ALL_PREFIX" --allow-extra-chr
fi

if [ ! -s "$ALL_PREFIX.eigenvec" ] || [ ! -s "$ALL_PREFIX.eigenval" ]; then
    "$PLINK" --bfile "$ALL_PREFIX" --pca "$PCA_COMPONENTS" \
        --threads "$PLINK_THREADS" --allow-extra-chr --out "$ALL_PREFIX"
fi

if ! cmp -s "$MODERN_PREFIX.bim" "$ALL_PREFIX.bim"; then
    echo "Error: modern and all-sample PLINK variant sets differ." >&2
    exit 1
fi

MODERN_P="$MODERN_PREFIX.${BEST_K}.P"
if [ ! -s "$MODERN_P" ]; then
    echo "Error: modern ADMIXTURE P file is missing: $MODERN_P" >&2
    exit 1
fi

cp "$MODERN_P" "$ALL_PREFIX.${BEST_K}.P.in"
printf '%s\n' "$BEST_K" > "$ANALYSIS_DIR/with_all_samples_optimal_K.txt"
PROJECTION_LOG="$ANALYSIS_DIR/with_all_samples_projection_K${BEST_K}.log"
PROJECTION_Q="$ALL_PREFIX.${BEST_K}.Q"
if [ ! -s "$PROJECTION_Q" ]; then
    "$ADMIXTURE" -j"$ADMIXTURE_THREADS" -P "$ALL_PREFIX.bed" "$BEST_K" \
        > "$PROJECTION_LOG" 2>&1
fi
if [ ! -s "$PROJECTION_Q" ]; then
    echo "Error: projection Q file was not generated: $PROJECTION_Q" >&2
    exit 1
fi

ANCIENT_TABLE="$ANALYSIS_DIR/ancient_projection_K${BEST_K}.tsv"
awk -v ancient="$ANCIENT_SAMPLES" -v k="$BEST_K" '
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
' "$ALL_PREFIX.fam" "$PROJECTION_Q" > "$ANCIENT_TABLE"

echo "Continuation complete."
echo "Best K: $BEST_K"
echo "All-sample PCA: $ALL_PREFIX.eigenvec"
echo "Ancient projection: $ANCIENT_TABLE"
