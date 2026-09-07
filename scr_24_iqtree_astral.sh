#!/bin/bash
# Build a modern-only backbone tree first, then build a second tree after
# adding ancient samples to the same modern-selected loci.

set -euo pipefail

# Configuration
# Usage: bash scr_24_iqtree_astral.sh [parallel_jobs]
THREADS="${1:-9}"
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
MODEL="GTR+G"
ASTRAL_JAR="${ASTRAL_JAR:-$HOME/Astral/astral.5.7.8.jar}"
GAP_THRESHOLD=30
DATASET_MODERN="modern"
DATASET_WITH_ALL="with_all_samples"
SCRIPT_PATH="$(readlink -f "$0")"

build_gene_tree() {
    local dataset=$1
    local gene=$2
    local align_dir="$WORK_DIR/final_alignments_${dataset}"
    local trees_dir="$WORK_DIR/gene_trees_${dataset}"
    local input="$align_dir/${gene}.fasta"
    local prefix="$trees_dir/$gene"
    local iqtree_input="$trees_dir/.${gene}.iqtree.fasta"

    if [ ! -s "$input" ]; then
        echo "Error: alignment not found or empty: $input" >&2
        return 1
    fi
    if [ -s "${prefix}.treefile" ]; then
        return 0
    fi

    # IQ-TREE rejects records made entirely of missing data. Omit only those
    # records; missing taxa are valid in the other gene trees.
    awk '
        function write_record(   informative) {
            if (header == "") return
            informative = sequence
            gsub(/[Nn?-.]/, "", informative)
            if (informative != "") print header "\n" sequence
        }
        /^>/ {
            write_record()
            header = $0
            sequence = ""
            next
        }
        { sequence = sequence $0 }
        END { write_record() }
    ' "$input" > "$iqtree_input"

    if ! iqtree -s "$iqtree_input" -st DNA -m "$MODEL" -nt 1 \
        --prefix "$prefix" -quiet; then
        rm -f "$iqtree_input"
        return 1
    fi
    rm -f "$iqtree_input"
}

if [ "${1:-}" = "__gene_tree" ]; then
    build_gene_tree "${2:-}" "${3:-}"
    exit $?
fi

build_dataset_trees() {
    local dataset=$1
    local align_dir="$WORK_DIR/final_alignments_${dataset}"
    local gene_list="$WORK_DIR/list_${GAP_THRESHOLD}_percent_gap_${dataset}.txt"
    local trees_dir="$WORK_DIR/gene_trees_${dataset}"
    local all_trees="$WORK_DIR/all_gene_trees_${dataset}.treefile"
    local astral_output="$WORK_DIR/astral_species_tree_${dataset}.tre"

    if [ ! -d "$align_dir" ] || [ ! -f "$gene_list" ]; then
        echo "Error: missing alignment directory or gene list for $dataset" >&2
        return 1
    fi
    mkdir -p "$trees_dir"

    grep -v '^$' "$gene_list" |
        xargs -P "$THREADS" -I {} bash "$SCRIPT_PATH" __gene_tree "$dataset" "{}"

    : > "$all_trees"
    while read -r gene; do
        [ -z "$gene" ] && continue
        treefile="$trees_dir/$gene.treefile"
        if [ ! -s "$treefile" ]; then
            echo "Error: gene tree not found or empty: $treefile" >&2
            return 1
        fi
        cat "$treefile" >> "$all_trees"
    done < "$gene_list"

    if [ ! -s "$all_trees" ]; then
        echo "Error: no gene trees generated for $dataset" >&2
        return 1
    fi
    if [ ! -f "$ASTRAL_JAR" ]; then
        echo "Error: ASTRAL jar not found: $ASTRAL_JAR" >&2
        return 1
    fi

    java -Xshare:off -Djava.io.tmpdir="${TMPDIR:-$HOME/tmp}" \
        -jar "$ASTRAL_JAR" -i "$all_trees" -o "$astral_output" \
        2> "$WORK_DIR/astral_${dataset}.log"
    echo "Species tree saved to: $astral_output"
}

build_dataset_trees "$DATASET_MODERN"
build_dataset_trees "$DATASET_WITH_ALL"
