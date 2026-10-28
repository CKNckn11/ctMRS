# ctMRS

`ctMRS` is an R package for integrating paired single-cell RNA/ATAC data,
transcription-factor motifs, and cohort somatic mutations. It builds a
full-genome TF-peak-gene regulatory network and calculates cell-type
specificity (CTS), mutation gene-regulatory scores (MGRS), and a combined
mutation regulatory score (MRS).

ctMRS assumes that recurrent plasma ctDNA mutations can be interpreted in
their cellular regulatory context by linking them to accessible chromatin
regions and TF-centered gene-regulatory programs. The framework first
characterizes the cohort-level ctDNA mutation landscape, then integrates
paired single-cell RNA-seq and ATAC-seq profiles with TF motif information to
construct a global TF-peak-gene network and resolve cell type-specific
eRegulons. Recurrent ctDNA variants are mapped to accessible regulatory peaks
and connected to candidate target genes through peak-gene associations.
Cell-type specificity (CTS) and mutation-associated gene-regulatory relevance
(MGRS) are subsequently integrated into the mutation-associated regulon score
(MRS). Matched Monte Carlo permutations are used to estimate empirical P
values and false-discovery rates and to prioritize significant
mutation-associated eRegulons.

The required inputs include paired single-cell RNA-seq and ATAC-seq data, TF
motif annotations, a GRCh38 gene annotation, and a cohort-level somatic
mutation table containing genomic coordinates and sample identifiers.

## Workflow

### Plasma ctDNA mutation profiling

Plasma cfDNA is analyzed with Mutect2 to identify somatic variants and
characterize the cohort-level ctDNA mutation landscape. Mutation relevance is
based on `log(1 + number of mutated samples)`. Functional annotations are
retained for description and are not assigned manual impact weights.

<p align="center">
  <img src="example/figure/Fig1A.png" alt="Plasma ctDNA mutation profiling" width="100%">
</p>

### Cell type-specific regulatory network construction

Paired single-cell RNA-seq and ATAC-seq profiles are integrated with TF motif
annotations and peak-gene associations to construct a global TF-peak-gene
regulatory network. The global network is resolved into cell type-specific
TF-centered eRegulons and used to quantify cell-type specificity (CTS).

<p align="center">
  <img src="example/figure/Fig1B.png" alt="Cell type-specific regulatory network construction" width="100%">
</p>

### Mutation-associated regulon scoring and significance testing

Recurrent ctDNA variants are mapped to accessible chromatin peaks and linked
to candidate target genes. Mutation-associated gene-regulatory relevance
(MGRS) is integrated with CTS to calculate the final mutation-associated
regulon score (MRS). Matched Monte Carlo permutations generate empirical P
values and false-discovery rates for prioritizing significant
mutation-associated eRegulons.

<p align="center">
  <img src="example/figure/Fig1C.png" alt="ctMRS scoring and matched-permutation workflow" width="68%">
</p>

```text
RNA + ATAC + motif -> TF-peak-gene global GRN -> cell-type sub-GRNs -> CTS
somatic variants -> cohort recurrence -> mutated ATAC peaks -> MGRS
CTS + MGRS - discordance penalty -> MRS -> matched-control Monte Carlo P/FDR
```

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("CKNckn11/ctMRS")
library(ctMRS)
```

## Required inputs

- A single-cell RNA Seurat object or RDS file.
- A single-cell ATAC Signac/Seurat object or RDS file containing a peak assay
  and motif annotations.
- GRCh38 gene annotation as a `GRanges` object or RDS file.
- A tab-delimited somatic-mutation table containing genomic coordinates,
  alleles, functional annotations, and sample identifiers.

## Example

```r
fit <- run_ctMRS(
  rna = "path/to/rna_object.rds",
  atac = "path/to/atac_object_with_motifs.rds",
  gene_annotation = "path/to/grch38_gene_annotation.rds",
  mutations = "path/to/somatic_mutations.tsv",
  output_dir = "path/to/ctMRS_results",
  celltype_column = "celltype",
  min_recurrent_samples = 1L,
  mutation_peak_buffer = 0L,
  perm_n = 1000L,
  seed = 1234L
)
```

Use `perm_n = 10000L` or more for a final analysis. The primary analysis can
retain `min_recurrent_samples = 1L`, because recurrence enters continuously
through `log1p(n_samples)`.

## Main functions

- `prepare_ctMRS_inputs()` pairs RNA/ATAC cells and extracts matrices.
- `calculate_gene_cts()` calculates gene cell-type specificity.
- `extract_tf_peak()` maps JASPAR motifs to TF-peak pairs.
- `link_peaks_genes()` computes full-genome RNA-ATAC peak-gene correlations.
- `build_global_grn()` builds the TF-peak-gene network.
- `build_celltype_subgrn()` calculates cell-type-specific regulatory edges.
- `calculate_mutation_relevance()` maps mutations to peaks and scores them.
- `calculate_target_mgrs()` calculates edge- and target-level MGRS.
- `score_eregulons()` calculates CTS, MGRS, and MRS.
- `test_eregulons_mc()` calculates empirical P values and FDR.
- `run_ctMRS()` runs the complete workflow.

## Author

Kangnan Chen (`CKNckn11`)  
GitHub: https://github.com/CKNckn11

## License

MIT
