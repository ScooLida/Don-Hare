#!/usr/bin/env bash
# Post-mapping mitochondrial pipeline:
# BAM -> variants -> low-depth-masked consensus -> 13 PCGs -> MAFFT ->
# IQ-TREE gene trees and ASTRAL species tree.
set -Eeuo pipefail
shopt -s nullglob globstar

# ---------- Edit these paths ----------
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REF="${REF:-/NatureUsers/ltursunova/hare_work/myto_hare.fasta}"
PCG_BED="${PCG_BED:-$SCRIPT_DIR/for_data/scr_31_hare_l_europaeus_NC_004028.1_PCGs.bed}"
BAM_ROOT="${BAM_ROOT:-/NatureUsers/ltursunova/hare_work/MyHare_myto/my_genome}"
SAMPLE_LIST="${SAMPLE_LIST:-/NatureUsers/ltursunova/hare_work/list_myto.txt}"

# Set this only when it is a BAM for one biological sample. Leave empty when
# all BAMs are already under BAM_ROOT. A merged multi-sample BAM must first be
# split by read-group/sample (SM tag), otherwise it will be called as one sample.
COMBINED_BAM="${COMBINED_BAM:-}"

OUT="${OUT:-/NatureUsers/ltursunova/hare_work/mito_pcg_analysis}"
THREADS="${THREADS:-8}"

# Ancient-DNA consensus filters. Adjust to the experiment if necessary.
MIN_MQ="${MIN_MQ:-25}"
MIN_BQ="${MIN_BQ:-20}"
MIN_DP="${MIN_DP:-2}"
MIN_QUAL="${MIN_QUAL:-30}"

# ASTRAL jar available on sy2. It can be overridden at runtime:
# ASTRAL_JAR=/actual/path/astral.jar bash scr_31_mito_pcg_astral.sh
ASTRAL_JAR="${ASTRAL_JAR:-$HOME/Astral/astral.5.7.8.jar}"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

# Run a program from the current environment when available. If it is absent,
# search all conda environments and run it with `conda run` from the first
# environment containing that executable. This avoids activating one global
# environment for the whole pipeline.
CONDA_BIN=""
if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    CONDA_BIN="$CONDA_EXE"
elif command -v conda >/dev/null 2>&1; then
    CONDA_BIN=$(command -v conda)
else
    for candidate in \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/mambaforge/bin/conda" \
        "$HOME/anaconda3/bin/conda"; do
        if [[ -x "$candidate" ]]; then
            CONDA_BIN="$candidate"
            break
        fi
    done
fi

declare -a CONDA_PREFIXES=()
if [[ -n "$CONDA_BIN" ]]; then
    mapfile -t CONDA_PREFIXES < <(
        "$CONDA_BIN" env list --json 2>/dev/null \
          | awk -F'"' '{for (i=2; i<=NF; i+=2) if ($i ~ /^\//) print $i}'
    )
fi
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    CONDA_PREFIXES=("$CONDA_PREFIX" "${CONDA_PREFIXES[@]}")
fi

declare -A TOOL_PREFIXES=()
resolve_tool() {
    local tool prefix
    tool="$1"
    [[ ${TOOL_PREFIXES[$tool]+present} ]] && return 0

    if command -v "$tool" >/dev/null 2>&1; then
        TOOL_PREFIXES[$tool]=""
        return 0
    fi

    for prefix in "${CONDA_PREFIXES[@]}"; do
        if [[ -x "$prefix/bin/$tool" ]]; then
            TOOL_PREFIXES[$tool]="$prefix"
            return 0
        fi
    done
    return 1
}

run_tool() {
    local tool prefix
    tool="$1"
    shift
    resolve_tool "$tool" || die "Command not found in current PATH or conda environments: $tool"
    prefix="${TOOL_PREFIXES[$tool]}"
    if [[ -n "$prefix" ]]; then
        printf '[env] %s <- %s\n' "$tool" "$prefix" >&2
        "$CONDA_BIN" run --no-capture-output -p "$prefix" "$tool" "$@"
    else
        command "$tool" "$@"
    fi
}

for cmd in samtools bcftools mafft java realpath awk tr rev sort; do
    resolve_tool "$cmd" || die "Command not found in current PATH or conda environments: $cmd"
done
if resolve_tool iqtree2; then
    IQTREE_CMD=iqtree2
elif resolve_tool iqtree; then
    IQTREE_CMD=iqtree
else
    die "Neither iqtree2 nor iqtree was found in current PATH or conda environments"
fi

[[ -s "$REF" ]] || die "Reference not found: $REF"
[[ -s "$PCG_BED" ]] || die "PCG BED not found: $PCG_BED"
[[ -d "$BAM_ROOT" ]] || die "BAM directory not found: $BAM_ROOT"
[[ -s "$SAMPLE_LIST" ]] || die "Sample list not found: $SAMPLE_LIST"
[[ -s "$ASTRAL_JAR" ]] || die "ASTRAL jar not found: $ASTRAL_JAR"

# Intermediate files are deliberately retained for concatenated analysis,
# troubleshooting, and reproducibility. This script never cleans OUT.
mkdir -p "$OUT" "$OUT/vcf" "$OUT/consensus" "$OUT/gene_fastas" \
    "$OUT/alignments" "$OUT/gene_trees" "$OUT/astral"

[[ -s "${REF}.fai" ]] || run_tool samtools faidx "$REF"

# BED format: gene<TAB>contig<TAB>start<TAB>end<TAB>strand
# Coordinates are 0-based and half-open. Keep overlapping genes as separate
# rows; each row is processed as an independent gene for ASTRAL.
mapfile -t GENES < <(run_tool awk 'BEGIN {FS="\t"} !/^#/ && NF >= 5 {print $1}' "$PCG_BED")
(( ${#GENES[@]} == 13 )) || die "PCG BED must contain exactly 13 genes; found ${#GENES[@]}"

run_tool awk '
    BEGIN { FS="\t" }
    !/^#/ && NF >= 5 {
        if ($3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ || $3 >= $4) {
            printf "Invalid BED coordinates at line %d\n", NR > "/dev/stderr"
            bad=1
        }
        if ($5 != "+" && $5 != "-") {
            printf "Invalid strand at line %d: %s\n", NR, $5 > "/dev/stderr"
            bad=1
        }
        if (++seen[$1] > 1) {
            printf "Duplicate gene in PCG BED: %s\n", $1 > "/dev/stderr"
            bad=1
        }
    }
    END { exit bad }
' "$PCG_BED" || die "Invalid PCG BED"

CONTIG=$(run_tool awk 'BEGIN {FS="\t"} !/^#/ && NF >= 5 {print $2; exit}' "$PCG_BED")
[[ -n "$CONTIG" ]] || die "Could not read the mitochondrial contig from PCG BED"

# Convert a BAM filename such as SRR5949627.rescaled.bam to SRR5949627.
sample_name_from_bam() {
    local bam base sample
    bam="$1"
    base=${bam##*/}
    sample=${base%.bam}
    sample=${sample%.rescaled}
    sample=$(printf '%s' "$sample" | run_tool tr -cs 'A-Za-z0-9_.-' '_')
    printf '%s\n' "$sample"
}

# Read only the first whitespace-delimited column from list_myto.txt.
declare -A REQUESTED_SAMPLES=()
while IFS= read -r sample; do
    [[ -n "$sample" ]] || continue
    REQUESTED_SAMPLES[$sample]=1
done < <(run_tool awk '!/^[[:space:]]*#/ && NF {print $1}' "$SAMPLE_LIST")
(( ${#REQUESTED_SAMPLES[@]} > 0 )) || die "Sample list is empty: $SAMPLE_LIST"

# Collect BAMs recursively and retain only samples listed in SAMPLE_LIST.
# The associative array avoids processing a BAM twice when COMBINED_BAM is
# also inside BAM_ROOT.
declare -a BAMS=()
declare -A SEEN_BAMS=()
declare -A FOUND_SAMPLES=()
add_bam() {
    local bam path sample
    bam="$1"
    [[ -f "$bam" ]] || return 0
    sample=$(sample_name_from_bam "$bam")
    [[ ${REQUESTED_SAMPLES[$sample]+present} ]] || return 0
    path=$(realpath "$bam")
    [[ ${SEEN_BAMS[$path]+present} ]] && return 0
    SEEN_BAMS[$path]=1
    BAMS+=("$path")
    FOUND_SAMPLES[$sample]=1
}

[[ -z "$COMBINED_BAM" ]] || add_bam "$COMBINED_BAM"
for bam in "$BAM_ROOT"/*.bam "$BAM_ROOT"/**/*.bam; do
    add_bam "$bam"
done
for sample in "${!REQUESTED_SAMPLES[@]}"; do
    [[ ${FOUND_SAMPLES[$sample]+present} ]] || die "No *rescaled.bam found for listed sample: $sample"
done
(( ${#BAMS[@]} > 0 )) || die "No requested BAM files found under $BAM_ROOT"

# Initialize per-gene FASTA files so rerunning the script does not append
# duplicate sequences.
for gene in "${GENES[@]}"; do
    : > "$OUT/gene_fastas/${gene}.fa"
done

declare -A SEEN_SAMPLES=()
for bam in "${BAMS[@]}"; do
    run_tool samtools quickcheck -v "$bam" || die "Corrupt or incomplete BAM: $bam"

    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
        run_tool samtools index -@ "$THREADS" "$bam"
    fi

    sample=$(sample_name_from_bam "$bam")
    [[ -n "$sample" ]] || die "Could not derive a sample name from $bam"
    [[ ! ${SEEN_SAMPLES[$sample]+present} ]] || die "Duplicate sample name: $sample"
    SEEN_SAMPLES[$sample]=1

    raw_vcf="$OUT/vcf/${sample}.raw.vcf.gz"
    filtered_vcf="$OUT/vcf/${sample}.filtered.vcf.gz"
    mask_bed="$OUT/vcf/${sample}.low_depth.bed"
    consensus="$OUT/consensus/${sample}.fa"
    mito_consensus="$OUT/consensus/${sample}.mitochondrion.fa"

    printf '[%s] Calling variants\n' "$sample"
    run_tool bcftools mpileup \
        --threads "$THREADS" \
        --fasta-ref "$REF" \
        --min-MQ "$MIN_MQ" \
        --min-BQ "$MIN_BQ" \
        --annotate FORMAT/DP,FORMAT/AD \
        --regions "$CONTIG" \
        --output-type u \
        "$bam" \
      | run_tool bcftools call \
        --threads "$THREADS" \
        --ploidy 1 \
        --multiallelic-caller \
        --output-type z \
        --output "$raw_vcf"
    run_tool bcftools index --force "$raw_vcf"

    run_tool bcftools filter \
        --include "QUAL>=${MIN_QUAL} && FORMAT/DP>=${MIN_DP}" \
        --output-type z \
        --output "$filtered_vcf" \
        "$raw_vcf"
    run_tool bcftools index --force "$filtered_vcf"

    # Positions below MIN_DP become N in consensus. Without this mask,
    # bcftools consensus would silently copy the reference at no-coverage sites.
    run_tool samtools depth \
        -aa -q "$MIN_MQ" -Q "$MIN_BQ" \
        -r "$CONTIG" "$bam" \
      | run_tool awk -v min="$MIN_DP" 'BEGIN {OFS="\t"} $3 < min {print $1, $2-1, $2}' \
      > "$mask_bed"

    # Also mask low-quality variant calls instead of silently turning them
    # into reference alleles in the consensus.
    run_tool bcftools query -f '%CHROM\t%POS\t%QUAL\n' "$raw_vcf" \
      | run_tool awk -v min="$MIN_QUAL" 'BEGIN {OFS="\t"} $3 < min {print $1, $2-1, $2}' \
      >> "$mask_bed"
    run_tool sort -k1,1 -k2,2n -u "$mask_bed" -o "$mask_bed"

    run_tool bcftools consensus \
        --fasta-ref "$REF" \
        --mask "$mask_bed" \
        "$filtered_vcf" \
      > "$consensus"

    # Keep only the mitochondrial contig for downstream extraction.
    run_tool samtools faidx "$consensus" "$CONTIG" > "$mito_consensus"
    run_tool samtools faidx "$mito_consensus"

    printf '[%s] Extracting 13 PCGs\n' "$sample"
    while IFS=$'\t' read -r gene chrom start end strand; do
        sequence=$(run_tool samtools faidx "$mito_consensus" \
            "${chrom}:$((start + 1))-${end}" \
          | run_tool awk 'NR > 1 { printf "%s", $0 }')
        [[ -n "$sequence" ]] || die "No sequence for $gene in $sample"
        if [[ "$strand" == "-" ]]; then
            sequence=$(printf '%s' "$sequence" | run_tool rev | run_tool tr 'ACGTacgt' 'TGCAtgca')
        fi
        printf '>%s\n%s\n' "$sample" "$sequence" >> "$OUT/gene_fastas/${gene}.fa"
    done < <(run_tool awk 'BEGIN {FS=OFS="\t"} !/^#/ && NF >= 5 {print $1, $2, $3, $4, $5}' "$PCG_BED")
done

printf '[alignment] Aligning each PCG with MAFFT\n'
for gene in "${GENES[@]}"; do
    run_tool mafft --auto --thread "$THREADS" \
        "$OUT/gene_fastas/${gene}.fa" \
      > "$OUT/alignments/${gene}.aln.fa"
done

printf '[gene trees] Running IQ-TREE ModelFinder for each PCG\n'
for gene in "${GENES[@]}"; do
    run_tool "$IQTREE_CMD" \
        -s "$OUT/alignments/${gene}.aln.fa" \
        -st DNA \
        -m MFP \
        -B 1000 \
        --alrt 1000 \
        -T "$THREADS" \
        --prefix "$OUT/gene_trees/${gene}"
done

all_gene_trees="$OUT/astral/all_gene_trees.treefile"
: > "$all_gene_trees"
for treefile in "$OUT"/gene_trees/*.treefile; do
    [[ -s "$treefile" ]] || continue
    run_tool awk '1' "$treefile" >> "$all_gene_trees"
done
[[ -s "$all_gene_trees" ]] || die "No gene trees were generated"

run_tool java -jar "$ASTRAL_JAR" \
    -i "$all_gene_trees" \
    -o "$OUT/astral/astral.tree"

printf '\nFinished.\nGene trees: %s\nASTRAL tree: %s\n' \
    "$OUT/gene_trees/" \
    "$OUT/astral/astral.tree"
