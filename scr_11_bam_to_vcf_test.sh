#!/usr/bin/env bash
# Test pipeline: call each sample in parallel, build a union of variant sites,
# then re-genotype every sample at the union sites before modern filtering.
# All outputs use the _t suffix and are isolated from the production pipeline.

set -euo pipefail

PARALLEL_JOBS="${PARALLEL_JOBS:-8}"
BCFTOOLS_THREADS="${BCFTOOLS_THREADS:-1}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-2}"
REFERENCE="${REFERENCE:-$HOME/hare_work/krol_g.fasta}"
BAM_DIR="${BAM_DIR:-$HOME/hare_work/My_new_krol/my_genome}"
SAMPLE_LIST="${SAMPLE_LIST:-sample_list.txt}"
CHROM_LIST="${CHROM_LIST:-chr.txt}"
BASE_PREFIX="${BASE_PREFIX:-MyHare}"
DATA_PREFIX="${BASE_PREFIX}_t"
VCF_DIR="${VCF_DIR:-./per_sample_vcfs_t}"
SELECTED_SITES="./${DATA_PREFIX}_modern_selected_sites.tsv"
QUAL=20
DEPTH=3
ANCIENT_SAMPLES=("1k" "3k" "4k" "5kS8")
SCRIPT_PATH="$(readlink -f "$0")"

export PARALLEL_JOBS BCFTOOLS_THREADS SAMTOOLS_THREADS

call_sample() {
    local sample=$1
    local mode=$2
    local union_sites=${3:-}
    local bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    local output

    if [ "$mode" = "initial" ]; then
        output="$VCF_DIR/${sample}_t_variants.vcf.gz"
        bcftools mpileup --threads "$BCFTOOLS_THREADS" -a FORMAT/DP -d 250 \
            -R "$CHROM_LIST" -f "$REFERENCE" "$bam" |
            bcftools call --threads "$BCFTOOLS_THREADS" -mv -Oz -o "$output"
    else
        output="$VCF_DIR/${sample}_t_on_union.vcf.gz"
        bcftools mpileup --threads "$BCFTOOLS_THREADS" -a FORMAT/DP -d 250 \
            -R "$union_sites" -f "$REFERENCE" "$bam" |
            bcftools call --threads "$BCFTOOLS_THREADS" -m -Oz -o "$output"
    fi
    bcftools index --force --tbi "$output"
}

if [ "${1:-}" = "__call_sample" ]; then
    call_sample "${2:-}" "${3:-initial}" "${4:-}"
    exit 0
fi

if ! [[ "$PARALLEL_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: PARALLEL_JOBS must be a positive integer: $PARALLEL_JOBS" >&2
    exit 1
fi
if ! [[ "$BCFTOOLS_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: BCFTOOLS_THREADS must be a positive integer: $BCFTOOLS_THREADS" >&2
    exit 1
fi

for command in bcftools samtools parallel; do
    command -v "$command" >/dev/null 2>&1 || {
        echo "Error: command not found: $command" >&2
        exit 1
    }
done
for file in "$SAMPLE_LIST" "$CHROM_LIST" "$REFERENCE"; do
    [ -f "$file" ] || {
        echo "Error: required file not found: $file" >&2
        exit 1
    }
done

mapfile -t SAMPLES < <(tr -d '\r' < "$SAMPLE_LIST" | awk 'NF {print $1}')
if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "Error: sample list is empty: $SAMPLE_LIST" >&2
    exit 1
fi

is_ancient() {
    local candidate=$1
    local ancient
    for ancient in "${ANCIENT_SAMPLES[@]}"; do
        [ "$candidate" = "$ancient" ] && return 0
    done
    return 1
}

mkdir -p "$VCF_DIR"
for sample in "${SAMPLES[@]}"; do
    bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    [ -f "$bam" ] || {
        echo "Error: BAM not found for $sample: $bam" >&2
        exit 1
    }
    if [ ! -f "$bam.bai" ] && [ ! -f "${bam%.bam}.bai" ] && [ ! -f "$bam.csi" ]; then
        samtools index -@ "$SAMTOOLS_THREADS" "$bam"
    fi
done

declare -a INITIAL_VCFS=()
for sample in "${SAMPLES[@]}"; do
    INITIAL_VCFS+=("$VCF_DIR/${sample}_t_variants.vcf.gz")
done
printf '%s\n' "${INITIAL_VCFS[@]}" > "${DATA_PREFIX}_initial_vcfs.txt"

printf '%s\n' "${SAMPLES[@]}" |
    parallel --tmpdir . --bar --halt soon,fail=1 -j "$PARALLEL_JOBS" \
    bash "$SCRIPT_PATH" __call_sample {} initial ""

UNION_VCF="${DATA_PREFIX}_union.vcf.gz"
bcftools merge --threads "$BCFTOOLS_THREADS" -Oz -o "$UNION_VCF" \
    "${INITIAL_VCFS[@]}"
bcftools index --force --tbi "$UNION_VCF"
bcftools query -f '%CHROM\t%POS\t%POS\n' "$UNION_VCF" |
    sort -k1,1 -k2,2n -u > "${DATA_PREFIX}_union_sites.tsv"
UNION_SITES="${DATA_PREFIX}_union_sites.tsv"
if [ ! -s "$UNION_SITES" ]; then
    echo "Error: no union variant sites were produced." >&2
    exit 1
fi

declare -a REGENOTYPED_VCFS=()
for sample in "${SAMPLES[@]}"; do
    REGENOTYPED_VCFS+=("$VCF_DIR/${sample}_t_on_union.vcf.gz")
done
printf '%s\n' "${REGENOTYPED_VCFS[@]}" > "${DATA_PREFIX}_regenotyped_vcfs.txt"

printf '%s\n' "${SAMPLES[@]}" |
    parallel --tmpdir . --bar --halt soon,fail=1 -j "$PARALLEL_JOBS" \
    bash "$SCRIPT_PATH" __call_sample {} union "$UNION_SITES"

REGENOTYPED_VCF="${DATA_PREFIX}_regenotyped_all.vcf.gz"
bcftools merge --threads "$BCFTOOLS_THREADS" -Oz -o "$REGENOTYPED_VCF" \
    "${REGENOTYPED_VCFS[@]}"
bcftools index --force --tbi "$REGENOTYPED_VCF"

ALL_SAMPLE_FILE="${DATA_PREFIX}_samples_all.txt"
MODERN_SAMPLE_FILE="${DATA_PREFIX}_samples_modern.txt"
ANCIENT_SAMPLE_FILE="${DATA_PREFIX}_samples_ancient.txt"
printf '%s\n' "${SAMPLES[@]}" > "$ALL_SAMPLE_FILE"
{
    for sample in "${SAMPLES[@]}"; do
        if ! is_ancient "$sample"; then printf '%s\n' "$sample"; fi
    done
} > "$MODERN_SAMPLE_FILE"
{
    for sample in "${SAMPLES[@]}"; do
        if is_ancient "$sample"; then printf '%s\n' "$sample"; fi
    done
} > "$ANCIENT_SAMPLE_FILE"

MERGED_MODERN_VCF="${DATA_PREFIX}_merged_modern.vcf.gz"
bcftools view -S "$MODERN_SAMPLE_FILE" -Oz -o "$MERGED_MODERN_VCF" "$REGENOTYPED_VCF"
bcftools index --force --tbi "$MERGED_MODERN_VCF"
bcftools view -e "QUAL < ${QUAL} || MIN(FMT/DP) < ${DEPTH}" -Oz \
    -o "${DATA_PREFIX}_modern.vcf.gz" "$MERGED_MODERN_VCF"
bcftools index --force --tbi "${DATA_PREFIX}_modern.vcf.gz"

bcftools query -f '%CHROM\t%POS\t%POS\n' "${DATA_PREFIX}_modern.vcf.gz" > "$SELECTED_SITES"
if [ ! -s "$SELECTED_SITES" ]; then
    echo "Error: no modern variants passed filtering." >&2
    exit 1
fi

ANCIENT_SELECTED_RAW="${DATA_PREFIX}_ancient_selected_raw.vcf.gz"
ANCIENT_SELECTED_MASKED="${DATA_PREFIX}_ancient_selected_masked.vcf.gz"
bcftools view -S "$ANCIENT_SAMPLE_FILE" -R "$SELECTED_SITES" -Oz \
    -o "$ANCIENT_SELECTED_RAW" "$REGENOTYPED_VCF"
bcftools index --force --tbi "$ANCIENT_SELECTED_RAW"
bcftools +setGT "$ANCIENT_SELECTED_RAW" -Oz \
    -o "$ANCIENT_SELECTED_MASKED" -- -t q -n . -i 'FMT/DP<2'
bcftools index --force --tbi "$ANCIENT_SELECTED_MASKED"

bcftools merge --threads "$BCFTOOLS_THREADS" -Oz \
    -o "${DATA_PREFIX}_with_all_samples.vcf.gz" \
    "${DATA_PREFIX}_modern.vcf.gz" "$ANCIENT_SELECTED_MASKED"
bcftools index --force --tbi "${DATA_PREFIX}_with_all_samples.vcf.gz"

echo "Test pipeline complete. Outputs use the _t suffix."
echo "Modern VCF: ${DATA_PREFIX}_modern.vcf.gz"
echo "All-sample VCF: ${DATA_PREFIX}_with_all_samples.vcf.gz"
