#!/usr/bin/env bash
# Run the complete mitochondrial PCG stage in order.
# scr_31 creates and retains all intermediate data; scr_32 then builds the
# primary partitioned concatenated ML tree from those alignments.
set -Eeuo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT="${OUT:-/NatureUsers/ltursunova/hare_work/mito_pcg_analysis}"
PCG_BED="${PCG_BED:-$SCRIPT_DIR/for_data/scr_31_hare_l_europaeus_NC_004028.1_PCGs.bed}"
REBUILD="${REBUILD:-0}"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

mkdir -p "$OUT"
LOG_FILE="$OUT/mito_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

printf '[mito pipeline] Output: %s\n' "$OUT"
printf '[mito pipeline] Intermediate files are preserved; no cleanup is performed.\n'
[[ -s "$PCG_BED" ]] || die "PCG BED not found: $PCG_BED"

has_files() {
    local directory="$1"
    [[ -d "$directory" ]] && compgen -G "$directory/*" >/dev/null
}

has_complete_pcg_inputs() {
    local gene
    [[ -s "$OUT/astral/astral.tree" ]] || return 1
    while IFS= read -r gene; do
        [[ -s "$OUT/alignments/${gene}.aln.fa" ]] || return 1
    done < <(awk 'BEGIN {FS="\t"} !/^#/ && NF >= 5 {print $1}' "$PCG_BED")
}

if [[ "$REBUILD" != 1 ]] && has_complete_pcg_inputs; then
    printf '[mito pipeline] Existing PCG intermediates found; reusing them.\n'
else
    if [[ "$REBUILD" != 1 ]]; then
        for directory in vcf consensus gene_fastas alignments gene_trees astral; do
            if has_files "$OUT/$directory"; then
                die "Existing mitochondrial output detected in $OUT/$directory; refusing to overwrite. Set REBUILD=1 only for an intentional rerun."
            fi
        done
    fi
    export OUT PCG_BED
    bash "$SCRIPT_DIR/scr_31_mito_pcg_astral.sh"
fi

for directory in vcf consensus gene_fastas alignments gene_trees astral; do
    [[ -d "$OUT/$directory" ]] || die "Missing output directory after scr_31: $OUT/$directory"
done

[[ -s "$OUT/astral/astral.tree" ]] || die "ASTRAL tree was not produced: $OUT/astral/astral.tree"

PCG_BED="$PCG_BED" bash "$SCRIPT_DIR/scr_32_mito_pcg_concat.sh" "$OUT"

[[ -s "$OUT/concatenated/pcg_concat.treefile" ]] || \
    die "Concatenated ML tree was not produced: $OUT/concatenated/pcg_concat.treefile"

printf '\n[mito pipeline] Finished successfully.\n'
printf '[mito pipeline] Preserved intermediate data: %s/{vcf,consensus,gene_fastas,alignments,gene_trees,astral}\n' "$OUT"
printf '[mito pipeline] Concatenated tree: %s\n' "$OUT/concatenated/pcg_concat.treefile"
