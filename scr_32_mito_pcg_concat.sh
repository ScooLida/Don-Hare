#!/usr/bin/env bash
# Concatenate the 13 non-overlapping mitochondrial PCGs with gene partitions
# and build the primary partitioned ML tree.
set -Eeuo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT="${1:-${MITO_OUT:-/NatureUsers/ltursunova/hare_work/mito_pcg_analysis}}"
ALIGN_DIR="${ALIGN_DIR:-$OUT/alignments}"
CONCAT_DIR="${CONCAT_DIR:-$OUT/concatenated}"
UNIQUE_PCG_BED="${UNIQUE_PCG_BED:-$OUT/pcg_unique_sites.bed}"
THREADS="${THREADS:-8}"
BOOTSTRAPS="${BOOTSTRAPS:-1000}"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

[[ -d "$ALIGN_DIR" ]] || die "Alignment directory not found: $ALIGN_DIR"
[[ -s "$UNIQUE_PCG_BED" ]] || die "Unique-site PCG BED not found: $UNIQUE_PCG_BED"
command -v python3 >/dev/null 2>&1 || die "python3 is required"

if [[ -n "${IQTREE_CMD:-}" ]]; then
    command -v "$IQTREE_CMD" >/dev/null 2>&1 || die "IQTREE_CMD not found: $IQTREE_CMD"
elif command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_CMD=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_CMD=iqtree
else
    CONDA_BIN=""
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
        CONDA_BIN="$CONDA_EXE"
    elif command -v conda >/dev/null 2>&1; then
        CONDA_BIN=$(command -v conda)
    fi
    [[ -n "$CONDA_BIN" ]] || die "Neither iqtree2 nor iqtree is available"

    mapfile -t CONDA_PREFIXES < <(
        "$CONDA_BIN" env list --json 2>/dev/null \
          | awk -F'"' '{for (i=2; i<=NF; i+=2) if ($i ~ /^\//) print $i}'
    )
    IQTREE_CMD=""
    for prefix in "${CONDA_PREFIXES[@]}"; do
        if [[ -x "$prefix/bin/iqtree2" ]]; then
            IQTREE_CMD="$prefix/bin/iqtree2"
            break
        elif [[ -x "$prefix/bin/iqtree" ]]; then
            IQTREE_CMD="$prefix/bin/iqtree"
            break
        fi
    done
    [[ -n "$IQTREE_CMD" ]] || die "Neither iqtree2 nor iqtree was found in conda environments"
fi

mapfile -t GENES < <(awk 'BEGIN {FS="\t"} !/^#/ && NF >= 5 {print $1}' "$UNIQUE_PCG_BED")
(( ${#GENES[@]} == 13 )) || die "Expected 13 PCGs, found ${#GENES[@]}"

mkdir -p "$CONCAT_DIR"
python3 "$SCRIPT_DIR/concat_fasta_partitions.py" \
    --alignment-dir "$ALIGN_DIR" \
    --genes "${GENES[@]}" \
    --output "$CONCAT_DIR/pcg_concat.fa" \
    --partitions "$CONCAT_DIR/pcg_partitions.txt" \
    --manifest "$CONCAT_DIR/pcg_manifest.tsv"

"$IQTREE_CMD" \
    -s "$CONCAT_DIR/pcg_concat.fa" \
    -st DNA \
    -p "$CONCAT_DIR/pcg_partitions.txt" \
    -m MFP+MERGE \
    -B "$BOOTSTRAPS" \
    --alrt "$BOOTSTRAPS" \
    -T "$THREADS" \
    --prefix "$CONCAT_DIR/pcg_concat"

printf '\nFinished.\nConcatenated alignment: %s\nPartition manifest: %s\nML tree: %s.treefile\n' \
    "$CONCAT_DIR/pcg_concat.fa" \
    "$CONCAT_DIR/pcg_manifest.tsv" \
    "$CONCAT_DIR/pcg_concat"
