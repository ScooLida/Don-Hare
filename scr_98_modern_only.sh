#!/bin/bash
# Run modern-only PCA and ADMIXTURE without requiring the all-sample VCF.

set -euo pipefail

ANALYSIS_DIR="./population_analysis"
MODERN_PREFIX="$ANALYSIS_DIR/modern"
MODERN_ADMIX_PREFIX="$ANALYSIS_DIR/modern_admix"
MODERN_VCF="${MODERN_VCF:-MyHare_modern.vcf.gz}"
CHROM_MAP="$ANALYSIS_DIR/admix_chromosomes.tsv"
UPDATE_CHR="$ANALYSIS_DIR/admix_update_chr.tsv"
PLINK="$HOME/plink"
ADMIXTURE="$HOME/admixture/dist/admixture_linux-1.3.0/admixture"
PLINK_THREADS=8
ADMIXTURE_THREADS=4
MISSINGNESS=0.2
PCA_COMPONENTS=4
K_MIN=2
K_MAX=15
LOG_FILE="$ANALYSIS_DIR/scr_98_modern_only.log"

mkdir -p "$ANALYSIS_DIR"
exec > >(tee "$LOG_FILE") 2>&1

echo "Modern-only analysis started: $(date)"
echo "Main log: $LOG_FILE"

REBUILD_MODERN="${REBUILD_MODERN:-0}"
if [ "$REBUILD_MODERN" = "1" ]; then
    if [ ! -s "$MODERN_VCF" ]; then
        echo "Error: modern VCF is missing: $MODERN_VCF" >&2
        exit 1
    fi
    echo "Rebuilding modern PLINK files from $MODERN_VCF"
    "$PLINK" --vcf "$MODERN_VCF" --geno "$MISSINGNESS" --make-bed \
        --set-missing-var-ids '@:#_$1_$2' --threads "$PLINK_THREADS" \
        --allow-extra-chr --out "$MODERN_PREFIX"
fi

if [ ! -s "$MODERN_PREFIX.bed" ] || [ ! -s "$MODERN_PREFIX.bim" ] || [ ! -s "$MODERN_PREFIX.fam" ]; then
    echo "Error: modern PLINK files are missing: $MODERN_PREFIX.*" >&2
    exit 1
fi

if [ "$REBUILD_MODERN" = "1" ] || [ ! -s "$MODERN_PREFIX.eigenvec" ] || [ ! -s "$MODERN_PREFIX.eigenval" ]; then
    echo "Modern PCA files are missing; running PCA."
    "$PLINK" --bfile "$MODERN_PREFIX" --pca "$PCA_COMPONENTS" \
        --threads "$PLINK_THREADS" --allow-extra-chr --out "$MODERN_PREFIX"
else
    echo "Using existing modern PCA files."
fi

# ADMIXTURE requires numeric chromosome codes. Keep the original PLINK files
# unchanged and create a separate numeric-chromosome copy for ADMIXTURE.
if [ "$REBUILD_MODERN" = "1" ] || [ ! -s "$CHROM_MAP" ]; then
    awk '!seen[$1]++ { print $1, ++n }' "$MODERN_PREFIX.bim" > "$CHROM_MAP"
fi
if [ "$REBUILD_MODERN" = "1" ] || [ ! -s "$UPDATE_CHR" ]; then
    awk 'NR == FNR { new[$1] = $2; next } ($1 in new) { print $2, new[$1] }' \
        "$CHROM_MAP" "$MODERN_PREFIX.bim" > "$UPDATE_CHR"
fi

if [ "$REBUILD_MODERN" = "1" ] || [ ! -s "$MODERN_ADMIX_PREFIX.bed" ] || [ ! -s "$MODERN_ADMIX_PREFIX.bim" ] || [ ! -s "$MODERN_ADMIX_PREFIX.fam" ]; then
    echo "Creating numeric-chromosome PLINK files for ADMIXTURE."
    "$PLINK" --bfile "$MODERN_PREFIX" --update-chr "$UPDATE_CHR" --make-bed \
        --allow-extra-chr --out "$MODERN_ADMIX_PREFIX"
else
    echo "Using existing numeric-chromosome PLINK files."
fi

CV_TABLE="$ANALYSIS_DIR/modern_admixture_cv.tsv"
printf "K\tCV_error\n" > "$CV_TABLE"
for K in $(seq "$K_MIN" "$K_MAX"); do
    log_file="${MODERN_ADMIX_PREFIX}_K${K}.log"
    if [ ! -s "$MODERN_ADMIX_PREFIX.$K.P" ] || [ ! -s "$MODERN_ADMIX_PREFIX.$K.Q" ] || ! grep -q 'CV error' "$log_file" 2>/dev/null; then
        echo "Starting modern ADMIXTURE K=$K"
        (
            cd "$ANALYSIS_DIR"
            "$ADMIXTURE" -j"$ADMIXTURE_THREADS" --cv "$(basename "$MODERN_ADMIX_PREFIX").bed" "$K" \
                > "$(basename "$log_file")" 2>&1
        )
    else
        echo "Using existing modern ADMIXTURE K=$K result."
    fi

    cv_error=$(awk '/CV error/ { value=$NF } END { print value }' "$log_file")
    if [ -z "$cv_error" ]; then
        echo "Error: no CV error found in $log_file" >&2
        exit 1
    fi
    printf "%s\t%s\n" "$K" "$cv_error" >> "$CV_TABLE"
    cp "$MODERN_ADMIX_PREFIX.$K.P" "$MODERN_PREFIX.$K.P"
    cp "$MODERN_ADMIX_PREFIX.$K.Q" "$MODERN_PREFIX.$K.Q"
    cp "$log_file" "$ANALYSIS_DIR/modern_admixture_K${K}.log"
done

echo "ADMIXTURE completed for K=${K_MIN}..${K_MAX}; K will be selected manually from the CV plot."

Rscript --vanilla - "$ANALYSIS_DIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- args[[1]]
plot_dir <- file.path(analysis_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}

pca_file <- file.path(analysis_dir, "modern.eigenvec")
pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(pca) < 4) stop("Unexpected PCA format: ", pca_file)
names(pca)[1:4] <- c("Family", "Sample", "PC1", "PC2")

ggplot2::ggsave(
  file.path(plot_dir, "modern_PCA.png"),
  ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::labs(title = "PCA: modern", x = "PC1", y = "PC2") +
    ggplot2::theme_bw(),
  width = 8,
  height = 6,
  dpi = 300
)

message("Modern PCA plot saved to: ", plot_dir)
RSCRIPT

echo "Modern-only analysis complete: $(date)"
echo "PCA plot: $ANALYSIS_DIR/plots/modern_PCA.png"
echo "ADMIXTURE results: $ANALYSIS_DIR/modern_admix_K{${K_MIN}..${K_MAX}}.Q"
