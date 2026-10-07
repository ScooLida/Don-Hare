# Mitochondrial Analysis Status

## Analyses

- `scr_31_mito_pcg_ml.sh` calls variants with `DP >= 3`, masks low-depth
  consensus positions, and removes overlapping reference sites from later
  PCGs. The unique-site BED is saved as `pcg_unique_sites.bed`.
- `scr_32_mito_pcg_concat.sh` builds the primary 13-PCG concatenated ML tree
  with one partition per PCG and `MFP+MERGE`, using SH-aLRT 1000 and UFBoot
  1000. ASTRAL and per-PCG gene-tree inference are not run.
- `scr_30_mito_pipeline.sh` runs both scripts in order and never cleans the
  output directory.

## Run commands

```bash
bash scr_30_mito_pipeline.sh
```

Configuration paths can be overridden without editing the scripts, for example:

```bash
OUT=/NatureUsers/ltursunova/hare_work/mito_pcg_analysis \
REF=/NatureUsers/ltursunova/hare_work/myto_hare.fasta \
BAM_ROOT=/NatureUsers/ltursunova/hare_work/MyHare_myto/my_genome \
SAMPLE_LIST=/NatureUsers/ltursunova/hare_work/list_myto.txt \
bash scr_30_mito_pipeline.sh
```

## Current limitation

The source PCG alignments, consensus FASTA files, and BAMs are not present
locally, so the mitochondrial analysis cannot be run correctly from the
current files alone.

The launcher refuses to overwrite an existing partial output. Use
`REBUILD=1 bash scr_30_mito_pipeline.sh` only when an intentional rerun from
BAM is required.

No nuclear tree is included in these routes.
