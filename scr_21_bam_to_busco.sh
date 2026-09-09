#!/bin/bash
# Extract BUSCO gene loci from BAM files and prepare two FASTA datasets:
# one containing modern samples only, and one with modern plus ancient samples.

set -euo pipefail

# Configuration
# Put one sample ID per line in SAMPLES. Empty lines and lines starting with
# '#' are ignored.
SAMPLES="${SAMPLES:-samples.txt}"
if [ ! -f "$SAMPLES" ] && [ -f "sample_list.txt" ]; then
    SAMPLES="sample_list.txt"
fi
# Previous inline list:
# SAMPLES="SRR14535670,SRR32541919"
ANCIENT_SAMPLES="1k,3k,4k,5kS8"
DATASET_MODERN="modern"
DATASET_WITH_ALL="with_all_samples"
REF_GENOME="${REF_GENOME:-$HOME/hare_work/krol_g.fasta}"
BED_FILE="${BED_FILE:-$HOME/hare_work/krol_genes_fixed.bed}"
BAM_DIR="${BAM_DIR:-$HOME/hare_work/My_new_krol/my_genome}"
OLD_GENES_DIR="./subset_parallel/genes_by_locus"

MIN_COVERAGE=2
NUM_THREADS=4
SLEEP_INTERVAL=0.1
MIN_SEQ_LENGTH=10

WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
EXTRACTED_DIR="$WORK_DIR/extracted_fasta"
EXTRACTED_GENES_LIST="$WORK_DIR/extracted_genes_list.txt"
COMBINED_MODERN_DIR="$WORK_DIR/combined_unaligned_${DATASET_MODERN}"
COMBINED_WITH_ALL_DIR="$WORK_DIR/combined_unaligned_${DATASET_WITH_ALL}"

mkdir -p "$EXTRACTED_DIR" "$COMBINED_MODERN_DIR" "$COMBINED_WITH_ALL_DIR"

display_time() {
    local seconds=$1
    printf "%02d:%02d:%02d" $((seconds / 3600)) $(((seconds % 3600) / 60)) $((seconds % 60))
}

is_ancient_sample() {
    local sample=$1
    [[ ",${ANCIENT_SAMPLES}," == *",${sample},"* ]]
}

filter_fasta_samples() {
    local input_file=$1
    local output_file=$2
    awk -v ancient="$ANCIENT_SAMPLES" '
        BEGIN {
            split(ancient, names, ",")
            for (i in names) ancient_sample[names[i]] = 1
        }
        /^>/ {
            sample = $1
            sub(/^>/, "", sample)
            keep = !(sample in ancient_sample)
        }
        keep { print }
    ' "$input_file" >> "$output_file"
}

> "$EXTRACTED_GENES_LIST"
if [ ! -f "${REF_GENOME}.fai" ]; then samtools faidx "$REF_GENOME"; fi

if [ ! -f "$SAMPLES" ]; then
    echo "Error: sample list not found: $SAMPLES" >&2
    exit 1
fi

mapfile -t SAMPLE_ARRAY < <(
    awk '{
        sub(/\r$/, "", $1)
        if ($1 != "" && $1 !~ /^#/) print $1
    }' "$SAMPLES"
)
if [ "${#SAMPLE_ARRAY[@]}" -eq 0 ]; then
    echo "Error: sample list is empty: $SAMPLES" >&2
    exit 1
fi
TOTAL_GENES=$(wc -l < "$BED_FILE")

for sample in "${SAMPLE_ARRAY[@]}"; do
    bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    if [ ! -f "$bam" ]; then
        echo "Error: BAM file not found: $bam"
        exit 1
    fi
    if [ ! -f "$bam.bai" ] && [ ! -f "${bam%.bam}.bai" ] && [ ! -f "$bam.csi" ]; then
        samtools index "$bam"
    fi
    mkdir -p "$EXTRACTED_DIR/$sample"

    job_count=0
    start_time=$(date +%s)
    while read -r chrom start end gene; do
        [ -z "$gene" ] && continue
        chrom=${chrom//$'\r'/}
        start=${start//$'\r'/}
        end=${end//$'\r'/}
        region="${chrom}:${start}-${end}"

        (
            temp_vcf="$EXTRACTED_DIR/$sample/${gene}_temp.vcf.gz"
            mask_bed="$EXTRACTED_DIR/$sample/${gene}_mask.bed"
            new_fasta="$EXTRACTED_DIR/$sample/${gene}.fasta"
            temp_ref="$EXTRACTED_DIR/$sample/${gene}_ref.fasta"

            if [ -s "$new_fasta" ] &&
                [ "$(grep -v '^>' "$new_fasta" | tr -d '\n' | wc -c)" -gt "$MIN_SEQ_LENGTH" ]; then
                exit 0
            fi

            samtools faidx "$REF_GENOME" "$region" > "$temp_ref"
            bcftools mpileup -Ou -f "$REF_GENOME" -r "$region" "$bam" |
                bcftools call -c -Oz -o "$temp_vcf"
            bcftools index -t "$temp_vcf"
            samtools depth -a -r "$region" "$bam" |
                awk -v min_cov="$MIN_COVERAGE" '$3 < min_cov { print $1 "\t" $2 - 1 "\t" $2 }' > "$mask_bed"

            printf ">%s\n" "$sample" > "$new_fasta"
            bcftools consensus -m "$mask_bed" -f "$temp_ref" "$temp_vcf" | grep -v '^>' >> "$new_fasta"
            rm -f "$temp_vcf" "${temp_vcf}.tbi" "$mask_bed" "$temp_ref"
        ) &

        sleep "$SLEEP_INTERVAL"
        job_count=$((job_count + 1))
        if (( job_count % NUM_THREADS == 0 )); then
            wait
            elapsed=$(( $(date +%s) - start_time ))
            [ "$elapsed" -gt 0 ] || elapsed=1
            remaining=$((TOTAL_GENES - job_count))
            echo "${sample}: $job_count/$TOTAL_GENES; elapsed $(display_time "$elapsed"); remaining approximately $remaining"
        fi
    done < "$BED_FILE"
    wait
done

for sample in "${SAMPLE_ARRAY[@]}"; do
    for fasta in "$EXTRACTED_DIR/$sample"/*.fasta; do
        [ -e "$fasta" ] || continue
        length=$(grep -v '^>' "$fasta" | tr -d '\n' | wc -c)
        if [ "$length" -gt "$MIN_SEQ_LENGTH" ]; then
            basename "$fasta" .fasta >> "$EXTRACTED_GENES_LIST"
        else
            rm -f "$fasta"
        fi
    done
done
sort -u "$EXTRACTED_GENES_LIST" -o "$EXTRACTED_GENES_LIST"

# Include loci already present in the legacy FASTA collection, even when a
# newly processed sample has no sequence for that locus.
for old_fasta in "$OLD_GENES_DIR"/*.fasta; do
    [ -e "$old_fasta" ] || continue
    basename "$old_fasta" .fasta >> "$EXTRACTED_GENES_LIST"
done
sort -u "$EXTRACTED_GENES_LIST" -o "$EXTRACTED_GENES_LIST"

while read -r gene; do
    [ -z "$gene" ] && continue
    combined_modern="$COMBINED_MODERN_DIR/${gene}.fasta"
    combined_all="$COMBINED_WITH_ALL_DIR/${gene}.fasta"
    : > "$combined_modern"
    : > "$combined_all"

    old_fasta="$OLD_GENES_DIR/${gene}.fasta"
    if [ -f "$old_fasta" ]; then
        cat "$old_fasta" >> "$combined_all"
        filter_fasta_samples "$old_fasta" "$combined_modern"
    fi

    for sample in "${SAMPLE_ARRAY[@]}"; do
        fasta="$EXTRACTED_DIR/$sample/${gene}.fasta"
        [ -f "$fasta" ] || continue
        cat "$fasta" >> "$combined_all"
        if ! is_ancient_sample "$sample"; then cat "$fasta" >> "$combined_modern"; fi
    done
done < "$EXTRACTED_GENES_LIST"

echo "BUSCO extraction and two FASTA datasets are ready."
echo "Modern samples only: $COMBINED_MODERN_DIR"
echo "With all samples: $COMBINED_WITH_ALL_DIR"
