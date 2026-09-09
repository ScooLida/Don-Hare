#!/usr/bin/env bash
set -euo pipefail

BAM_ROOT="${BAM_ROOT:-$HOME/hare_work/My_new_krol/my_genome}"
SAMPLE_LIST="${SAMPLE_LIST:-}"
FILTERED_DIR="${FILTERED_DIR:-./sex_filtered_bams}"
OUTPUT="${OUTPUT:-sex_results.txt}"
MAPQ=20
SAMPLES="${SAMPLES:-1k,4k,3k}"

# X-хромосома
X="NC_091453.1"

# Митохондриальная хромосома,
MT=""

# Общий файл с результатами
result_file="$OUTPUT"

# Заголовок
echo -e "Sample\tX_depth\tAutosomal_mean_depth\tX_autosome_ratio\tSex" > "$result_file"

# Разделяем список по запятым
if [[ -n "$SAMPLE_LIST" && -s "$SAMPLE_LIST" ]]; then
    mapfile -t sample_array < <(awk '!/^[[:space:]]*#/ && NF {print $1}' "$SAMPLE_LIST")
else
    IFS=',' read -ra sample_array <<< "$SAMPLES"
fi
mkdir -p "$FILTERED_DIR"

for sample in "${sample_array[@]}"
do
    echo "Processing: $sample"

    bam="$BAM_ROOT/$sample/$sample.rescaled.bam"
    filtered_bam="$FILTERED_DIR/${sample}.filtered.bam"
    coverage_file="$FILTERED_DIR/${sample}.coverage.txt"
    if [[ ! -s "$bam" ]]; then
        echo "Error: BAM not found: $bam" >&2
        exit 1
    fi

    # 1. Фильтрация BAM
    # MAPQ >= 20
    # убрать unmapped, secondary, QC-fail, duplicates, supplementary
    samtools view -b -q "$MAPQ" -F 3844 "$bam" > "$filtered_bam"

    # 2. Индексация
    samtools index "$filtered_bam"

    # 3. Coverage по хромосомам
    samtools coverage "$filtered_bam" \
        > "$coverage_file"
    # 4. X/autosome ratio и запись результата
    awk -v SAMPLE="$sample" \
        -v X="$X" \
        -v MT="$MT" '
    $1 == X {
        x=$7
        x_seen=1
    }

    $1 ~ /^NC_/ && $1 != X && (MT=="" || $1 != MT) {
        sum += $7
        n++
    }

    END {
        if (!x_seen || n == 0 || sum == 0) exit 1
        auto=sum/n
        ratio=x/auto

        if (ratio < 0.65)
            sex="MALE"
        else if (ratio > 0.80)
            sex="FEMALE"
        else
            sex="AMBIGUOUS"

        print SAMPLE "\t" x "\t" auto "\t" ratio "\t" sex
    }' "$coverage_file" >> "$result_file"

done

echo "Done. Results:"
column -t "$result_file"
