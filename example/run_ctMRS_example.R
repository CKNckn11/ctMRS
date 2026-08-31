library(ctMRS)

# Test run: perm_n = 1000L
# Final publication run: perm_n >= 10000L
fit <- run_ctMRS(
  rna = "path/to/rna_object.rds",
  atac = "path/to/atac_object_with_motifs.rds",
  gene_annotation = "path/to/grch38_gene_annotation.rds",
  mutations = "path/to/somatic_mutations.tsv",
  output_dir = "path/to/ctMRS_results",
  rna_assay = "RNA",
  peak_assay = "ATAC",
  celltype_column = "celltype",
  min_recurrent_samples = 1L,
  mutation_peak_buffer = 0L,
  perm_n = 1000L,
  seed = 1234L
)

head(fit$final_results)
