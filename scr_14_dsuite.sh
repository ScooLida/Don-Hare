#!/bin/bash
# Run Dsuite Dtrios for both population/species assignments.

set -euo pipefail

DSUITE="$HOME/Dsuite/Build/Dsuite"
VCF="${1:-MyHare_t_with_all_samples.vcf.gz}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SET_DIR="$SCRIPT_DIR/for_data"
OUT_DIR="${OUT_DIR:-./dsuite_results_t}"

if [ ! -x "$DSUITE" ]; then
    echo "Error: Dsuite executable not found or not executable: $DSUITE" >&2
    exit 1
fi
if [ ! -f "$VCF" ]; then
    echo "Error: VCF not found: $VCF" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

for sets in "$SET_DIR/sets_hare.txt" "$SET_DIR/sets_krol.txt"; do
    if [ ! -s "$sets" ]; then
        echo "Error: SETS file not found or empty: $sets" >&2
        exit 1
    fi

    name=$(basename "$sets" .txt)
    prefix="$OUT_DIR/${name}_t"
    echo "Running Dsuite with $sets"
    "$DSUITE" Dtrios -o "$prefix" "$VCF" "$sets"
done

echo "Dsuite analyses complete: $OUT_DIR"
