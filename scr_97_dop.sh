#!/bin/bash
# Continue scr_12 without recalculating modern PCA. Reuse valid modern
# ADMIXTURE results, or calculate them once with numeric chromosome codes.
# Run from hare_work after MyHare_with_all_samples.vcf.gz is ready.

set -euo pipefail

ANALYSIS_DIR="./population_analysis"
ALL_VCF="MyHare_with_all_samples.vcf.gz"
MODERN_PREFIX="$ANALYSIS_DIR/modern"
ALL_PREFIX="$ANALYSIS_DIR/with_all_samples"
MODERN_ADMIX_PREFIX="$ANALYSIS_DIR/modern_admix"
ALL_ADMIX_PREFIX="$ANALYSIS_DIR/with_all_samples_admix"
CHROM_MAP="$ANALYSIS_DIR/admix_chromosomes.tsv"
UPDATE_CHR="$ANALYSIS_DIR/admix_update_chr.tsv"
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

# ADMIXTURE requires numeric chromosome codes. Keep the original PLINK files
# unchanged and create separate numeric-chromosome copies for ADMIXTURE.
if [ ! -s "$CHROM_MAP" ]; then
    awk '!seen[$1]++ { print $1, ++n }' "$MODERN_PREFIX.bim" > "$CHROM_MAP"
fi
if [ ! -s "$UPDATE_CHR" ]; then
    awk 'NR == FNR { new[$1] = $2; next } ($1 in new) { print $2, new[$1] }' \
        "$CHROM_MAP" "$MODERN_PREFIX.bim" > "$UPDATE_CHR"
fi

if [ ! -s "$MODERN_ADMIX_PREFIX.bed" ] || [ ! -s "$MODERN_ADMIX_PREFIX.bim" ] || [ ! -s "$MODERN_ADMIX_PREFIX.fam" ]; then
    "$PLINK" --bfile "$MODERN_PREFIX" --update-chr "$UPDATE_CHR" --make-bed \
        --allow-extra-chr --out "$MODERN_ADMIX_PREFIX"
fi

CV_TABLE="$ANALYSIS_DIR/modern_admixture_cv.tsv"
printf "K\tCV_error\n" > "$CV_TABLE"
for K in $(seq "$K_MIN" "$K_MAX"); do
    log_file="${MODERN_ADMIX_PREFIX}_K${K}.log"
    if [ ! -s "$MODERN_ADMIX_PREFIX.$K.P" ] || [ ! -s "$MODERN_ADMIX_PREFIX.$K.Q" ] || ! grep -q 'CV error' "$log_file" 2>/dev/null; then
        (
            cd "$ANALYSIS_DIR"
            "$ADMIXTURE" -j"$ADMIXTURE_THREADS" --cv "$(basename "$MODERN_ADMIX_PREFIX").bed" "$K" \
                > "$(basename "$log_file")" 2>&1
        )
    fi
    cv_error=$(awk '/CV error/ { value=$NF } END { print value }' "$log_file")
    if [ -z "$cv_error" ]; then
        echo "Error: no CV error found in $log_file" >&2
        echo "Inspect this log before continuing." >&2
        exit 1
    fi
    printf "%s\t%s\n" "$K" "$cv_error" >> "$CV_TABLE"
    cp "$MODERN_ADMIX_PREFIX.$K.P" "$MODERN_PREFIX.$K.P"
    cp "$MODERN_ADMIX_PREFIX.$K.Q" "$MODERN_PREFIX.$K.Q"
done

echo "ADMIXTURE completed for K=${K_MIN}..${K_MAX}; K will be selected manually from the CV plot."

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

if [ ! -s "$ALL_ADMIX_PREFIX.bed" ] || [ ! -s "$ALL_ADMIX_PREFIX.bim" ] || [ ! -s "$ALL_ADMIX_PREFIX.fam" ]; then
    "$PLINK" --bfile "$ALL_PREFIX" --update-chr "$UPDATE_CHR" --make-bed \
        --allow-extra-chr --out "$ALL_ADMIX_PREFIX"
fi
if ! cmp -s "$MODERN_ADMIX_PREFIX.bim" "$ALL_ADMIX_PREFIX.bim"; then
    echo "Error: numeric modern and all-sample PLINK variant sets differ." >&2
    exit 1
fi

for K in $(seq "$K_MIN" "$K_MAX"); do
    MODERN_P="$MODERN_ADMIX_PREFIX.${K}.P"
    if [ ! -s "$MODERN_P" ]; then
        echo "Error: modern ADMIXTURE P file is missing: $MODERN_P" >&2
        exit 1
    fi

    cp "$MODERN_P" "$ALL_ADMIX_PREFIX.${K}.P.in"
    PROJECTION_LOG="$ANALYSIS_DIR/with_all_samples_projection_K${K}.log"
    PROJECTION_Q="$ALL_ADMIX_PREFIX.${K}.Q"
    if [ ! -s "$PROJECTION_Q" ]; then
        (
            cd "$ANALYSIS_DIR"
            "$ADMIXTURE" -j"$ADMIXTURE_THREADS" -P "$(basename "$ALL_ADMIX_PREFIX").bed" "$K" \
                > "$(basename "$PROJECTION_LOG")" 2>&1
        )
    fi
    if [ ! -s "$PROJECTION_Q" ]; then
        echo "Error: projection Q file was not generated: $PROJECTION_Q" >&2
        exit 1
    fi
    cp "$PROJECTION_Q" "$ALL_PREFIX.${K}.Q"

    ANCIENT_TABLE="$ANALYSIS_DIR/ancient_projection_K${K}.tsv"
    awk -v ancient="$ANCIENT_SAMPLES" -v k="$K" '
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
    echo "Ancient projection K=$K: $ANCIENT_TABLE"
done

echo "Continuation complete."
echo "All-sample PCA: $ALL_PREFIX.eigenvec"
echo "All K projections: $ANALYSIS_DIR/ancient_projection_K*.tsv"
