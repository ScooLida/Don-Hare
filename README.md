# Don-Hare Analysis Pipelines

The repository contains two independent analysis lines for hare genomic data.
The population-genetic line produces filtered VCF files, PCA, and ADMIXTURE
results. The BUSCO line extracts BUSCO gene sequences, filters and aligns them,
and builds gene and species trees.

The samples `1k`, `3k`, `4k`, and `5kS8` are ancient samples. Modern-only and
all-sample datasets are processed separately. Locus selection for the tree
analysis is performed using modern samples only; ancient sequences are then
added to the same modern-selected loci.

## Line 1: VCF and Population Analyses

Run the scripts in this order:

### `scr_11_bam_to_vcf.sh`

Calls variants jointly from Paleomix BAM files, restricts variant calling to
regions listed in `chr.txt`, and keeps per-sample VCFs for individual use.
Modern samples then define the filtered variant positions. The all-sample VCF
is restricted to exactly the same positions, while ancient samples with no call
remain missing (`./.`).

Filters: `QUAL >= 20` and `MIN(FMT/DP) >= 3`, meaning every modern sample must
have at least three reads at the site. PLINK additionally applies `--geno 0.2`
to modern samples. No LD-pruning filter is applied. The joint call is built
from all sample BAMs so non-variant reference genotypes are not converted to
missing values by merging per-sample VCFs.

After site selection, ancient genotypes with `DP < 2` are set to `./.` while
the site itself is retained for modern samples and other ancient samples.

### `scr_12_pca_admixture.sh`

Runs PCA on modern samples using the modern-selected positions, then projects
all samples onto the modern PCA axes with PLINK 2 allele weights. PLINK applies
`--geno 0.2` only to modern samples. ADMIXTURE tests `K=2..15` on modern
samples and projects all samples for every K using the fixed modern `.P`
matrix; K is selected manually from the CV plot. Ancient-only results are
saved to `population_analysis/ancient_projection_K<K>.tsv` and
`population_analysis/ancient_pca_projection.tsv`.

### `scr_13_population_plots.R`

Creates PCA scatter plots and ADMIXTURE ancestry-proportion plots for both
datasets using the outputs from `scr_12_pca_admixture.sh`.

### `scr_14_dsuite.sh`

Runs `Dsuite Dtrios` on `MyHare_with_all_samples.vcf.gz` using both
`for_data/sets_hare.txt` and `for_data/sets_krol.txt`. The set files define the
population/species assignments and outgroup; samples marked `xxx` are ignored
by Dsuite. Results are written to `dsuite_results/sets_hare_*` and
`dsuite_results/sets_krol_*`.

### `scr_15_fbranch_new_trees.sh`

Uses the grouped ASTRAL tree
`astral_species_tree_with_all_samples_grouped_modern_ancient_separate.tre` as
the topology source. The tree leaves are population/species labels, not
individual sample IDs. Dsuite gets the sample-to-population assignments from
the second column of `sets_hare.txt` or `sets_krol.txt`, so a sample does not
need to occur as an individual leaf: it is included when its assigned
population/species leaf is present in the tree. Samples assigned `xxx` are
ignored. The script verifies that every non-`xxx` population label required by
the set file is present in the grouped tree before running `Dtrios` and
`Fbranch`.

## Line 2: BUSCO Gene Analysis and Trees

Run the scripts in this order:

### `scr_21_bam_to_busco.sh`

Extracts BUSCO gene loci from BAM files and prepares two merged FASTA datasets:
one containing modern samples only, and one containing modern plus ancient
samples. Existing per-sample FASTA files are skipped when they are non-empty
and contain more than `MIN_SEQ_LENGTH` nucleotides.

Filters and settings: `MIN_COVERAGE=2`, discard sequences shorter than
`MIN_SEQ_LENGTH=10`, and use BUSCO coordinates from `BED_FILE`.

### `scr_22_mafft.sh`

Aligns modern FASTA files first with MAFFT using nine parallel jobs, then adds
ancient sequences to the modern alignments with `mafft --addfragments`. This
preserves the modern alignment as the backbone. If an ancient sample has no
sequence for a selected locus, the all-sample alignment contains an `N` sequence
of the same length as the modern alignment. Loci without any modern sequence
are skipped because they cannot define the modern backbone. Modern input
sequences longer than 10000 nt are rejected before MAFFT. An alignment that
exceeds 10000 nt after MAFFT is handled by the later `scr_23` filter. Completed
modern alignments are reused when the script is restarted. The number of
parallel jobs can be supplied at launch, for example `bash scr_22_mafft.sh 2`;
the default is `9`.

### `scr_23_filter_by_gaps.sh`

Filters aligned BUSCO genes using the modern alignments only and writes the
modern-selected locus list. The same list is copied for the all-sample dataset;
absence of an ancient sequence does not remove the locus.

Filters: maximum sequence length `5000`, minimum ATGC count per modern sample
`30`, and maximum non-ATGC proportion `30%`.

### `scr_24_iqtree_astral.sh`

Builds the modern-only IQ-TREE gene trees first with the `GTR+G` model and runs
ASTRAL to produce the modern backbone species tree. It then builds a second set
of gene trees and an all-sample species tree using the same modern-selected
loci after ancient sequences have been added. The second tree is a combined
inference and may change the topology; it is not formal fixed-tree placement.

## Line 3: Mitochondrial PCG Analysis

### `scr_30_mito_pipeline.sh`

Runs the mitochondrial stage in order:

1. `scr_31_mito_pcg_astral.sh` calls variants, creates masked consensus
   sequences, extracts the 13 PCGs, aligns them, and builds the per-PCG trees
   plus the ASTRAL tree.
2. `scr_32_mito_pcg_concat.sh` concatenates the 13 retained PCG alignments
   using one partition per gene and builds the primary ML tree.

All intermediate data are retained under `OUT`: `vcf`, `consensus`,
`gene_fastas`, `alignments`, `gene_trees`, and `astral`. The launcher refuses
to overwrite a partial existing output. Set `REBUILD=1` only for an intentional
rerun from BAM files.

Run on the compute environment with the default `/NatureUsers` paths:

```bash
bash scr_30_mito_pipeline.sh
```

Paths can be overridden through `REF`, `BAM_ROOT`, `SAMPLE_LIST`, `OUT`, and
`PCG_BED` environment variables.

## Other Files

### `scr_98_sex.sh`

An independent sex-analysis script based on X-chromosome and autosomal
coverage. BAM filtering and indexing commands are currently commented out and
must be enabled or revised before use.

### `scr_99_old_ngs_analysis`

Archived older analysis material; it is not part of the two main workflows.
