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
PLINK2="${PLINK2:-$HOME/plink2}"
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

    echo "${label} VCF variants: $(bcftools index -n "$input_vcf")"

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

    echo "${label} PLINK variants after conversion: $(wc -l < "${prefix}.bim")"
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

    echo "ADMIXTURE completed for K=${K_MIN}..${K_MAX}; K will be selected manually from the CV plot."
}

run_pca_projection() {
    local modern_prefix="$ANALYSIS_DIR/modern"
    local all_prefix="$ANALYSIS_DIR/with_all_samples"
    local pca_prefix="$ANALYSIS_DIR/modern_pca"
    local projection_prefix="$ANALYSIS_DIR/with_all_samples_pca_projection"
    local score_end=$((5 + PCA_COMPONENTS))

    if [ ! -x "$PLINK2" ]; then
        echo "Error: PLINK 2 is required for PCA projection: $PLINK2" >&2
        echo "Set PLINK2=/path/to/plink2 or install it at $HOME/plink2." >&2
        return 1
    fi
    if ! cmp -s "${modern_prefix}.bim" "${all_prefix}.bim"; then
        echo "Error: modern and all-sample PLINK SNP sets differ." >&2
        return 1
    fi

    "$PLINK2" --bfile "$modern_prefix" --freq counts \
        --pca allele-wts "$PCA_COMPONENTS" vcols=chrom,ref,alt \
        --out "$pca_prefix"
    "$PLINK2" --bfile "$all_prefix" --read-freq "${pca_prefix}.acount" \
        --score "${pca_prefix}.eigenvec.allele" 2 5 header-read \
        no-mean-imputation variance-standardize \
        --score-col-nums "6-${score_end}" --out "$projection_prefix"

    python3 - "$projection_prefix.sscore" "$ANCIENT_SAMPLES" \
        "$ANALYSIS_DIR/ancient_pca_projection.tsv" <<'PY'
import sys

score_file, ancient_text, output_file = sys.argv[1:]
ancient = set(ancient_text.split(","))
with open(score_file) as source:
    header = source.readline().rstrip("\n").split("\t")
    clean_header = [value.lstrip("#") for value in header]
    iid_index = clean_header.index("IID")
    pc_indices = [i for i, value in enumerate(clean_header) if value.endswith("_AVG")]
    with open(output_file, "w") as output:
        output.write("sample\t" + "\t".join(clean_header[i] for i in pc_indices) + "\n")
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if fields[iid_index] in ancient:
                output.write(fields[iid_index] + "\t" + "\t".join(fields[i] for i in pc_indices) + "\n")
PY

    echo "Modern PCA reference: ${pca_prefix}.eigenvec"
    echo "All-sample PCA projection: ${projection_prefix}.sscore"
    echo "Ancient PCA projection: $ANALYSIS_DIR/ancient_pca_projection.tsv"
}

project_all_samples() {
    local modern_prefix="$ANALYSIS_DIR/modern"
    local all_prefix="$ANALYSIS_DIR/with_all_samples"

    # Projection requires identical SNP IDs and order in both datasets.
    if ! cmp -s "${modern_prefix}.bim" "${all_prefix}.bim"; then
        echo "Error: modern and all-sample PLINK SNP sets differ." >&2
        return 1
    fi

    for K in $(seq "$K_MIN" "$K_MAX"); do
        cp "${modern_prefix}.${K}.P" "${all_prefix}.${K}.P.in"
        projection_log="$ANALYSIS_DIR/with_all_samples_projection_K${K}.log"
        "$ADMIXTURE" -j"$ADMIXTURE_THREADS" -P "$all_prefix.bed" "$K" \
            > "$projection_log" 2>&1

        projection_q="${all_prefix}.${K}.Q"
        if [ ! -s "$projection_q" ]; then
            echo "Error: projection Q file was not generated: $projection_q" >&2
            return 1
        fi

        ancient_table="$ANALYSIS_DIR/ancient_projection_K${K}.tsv"
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
        ' "${all_prefix}.fam" "$projection_q" > "$ancient_table"
        echo "Ancient ADMIXTURE projection K=$K: $ancient_table"
    done
}

mkdir -p "$ANALYSIS_DIR"
prepare_dataset "modern" "$MODERN_VCF"
run_modern_admixture

prepare_dataset "with_all_samples" "$ALL_VCF"
run_pca_projection
project_all_samples

echo "PCA and fixed-model ADMIXTURE analyses complete."
