#!/usr/bin/env bash

set -euo pipefail

BAM_ROOT="${BAM_ROOT:-$HOME/hare_work/My_new_krol/my_genome}"
OUTPUT="${OUTPUT:-$HOME/hare_work/ancient_rescaled_summary.tsv}"
SAMPLES=(1k 3k 4k 5kS8)

printf 'sample\tbam\ttotal_reads\tprimary_reads\tmapped_reads\tmapped_percent\tduplicates\tproperly_paired\treference_bases\tcovered_bases\tbreadth_percent\tmean_depth\tmean_baseq\tmean_mapq\n' > "$OUTPUT"

for sample in "${SAMPLES[@]}"; do
    bam="$BAM_ROOT/$sample/$sample.rescaled.bam"
    if [[ ! -s "$bam" ]]; then
        echo "Error: BAM not found: $bam" >&2
        exit 1
    fi
    if ! samtools quickcheck "$bam"; then
        echo "Error: BAM failed quickcheck: $bam" >&2
        exit 1
    fi

    flagstat=$(samtools flagstat "$bam")
    total_reads=$(awk '/^[0-9]+ \+ [0-9]+ in total/ {print $1; exit}' <<< "$flagstat")
    primary_reads=$(awk '/^[0-9]+ \+ [0-9]+ primary$/ {print $1; exit}' <<< "$flagstat")
    mapped_reads=$(awk '/^[0-9]+ \+ [0-9]+ mapped \(/ {print $1; exit}' <<< "$flagstat")
    mapped_percent=$(awk '/^[0-9]+ \+ [0-9]+ mapped \(/ {value=$5; gsub(/[()%]/, "", value); print value; exit}' <<< "$flagstat")
    duplicates=$(awk '/^[0-9]+ \+ [0-9]+ duplicates$/ {print $1; exit}' <<< "$flagstat")
    properly_paired=$(awk '/^[0-9]+ \+ [0-9]+ properly paired/ {print $1; exit}' <<< "$flagstat")

    coverage=$(samtools coverage "$bam" | awk '
        /^#/ { next }
        NF >= 9 {
            length_bp = $3 - $2 + 1
            reference_bases += length_bp
            covered_bases += $5
            depth_sum += $7 * length_bp
            baseq_sum += $8 * length_bp
            mapq_sum += $9 * length_bp
        }
        END {
            if (reference_bases == 0) {
                print "0\t0\t0\t0\t0\t0"
            } else {
                printf "%.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\n", reference_bases, covered_bases, 100 * covered_bases / reference_bases, depth_sum / reference_bases, baseq_sum / reference_bases, mapq_sum / reference_bases
            }
        }')

    IFS=$'\t' read -r reference_bases covered_bases breadth_percent mean_depth mean_baseq mean_mapq <<< "$coverage"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$bam" "$total_reads" "$primary_reads" "$mapped_reads" "$mapped_percent" \
        "$duplicates" "$properly_paired" "$reference_bases" "$covered_bases" "$breadth_percent" \
        "$mean_depth" "$mean_baseq" "$mean_mapq" >> "$OUTPUT"

    echo "Processed: $sample"
done

echo "Saved: $OUTPUT"
