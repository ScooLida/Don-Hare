# Mitochondrial Analysis Status

## Analyses

- `scr_31_mito_pcg_astral.sh` remains the existing per-PCG analysis.
- `scr_32_mito_pcg_concat.sh` builds the primary 13-PCG concatenated ML tree
  with one partition per PCG. It leaves the ASTRAL tree unchanged.
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

The working copy currently contains only the 13 gene trees under
`PycharmProjects/hares/mito_pcg_analysis/gene_trees/`. The source PCG
alignments, consensus FASTA files, and BAMs are not present locally, so the
concatenated ML analysis cannot be run correctly from the current files alone.

The launcher refuses to overwrite an existing partial output. Use
`REBUILD=1 bash scr_30_mito_pipeline.sh` only when an intentional rerun from
BAM is required.

No nuclear tree is included in these routes.
