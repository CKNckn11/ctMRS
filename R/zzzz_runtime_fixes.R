# data.table column symbols are evaluated inside data.table expressions.
utils::globalVariables(c(
  ".", ".N", ":=", "..required", "N", "TF", "TF_raw", "TF_CTS",
  "Celltype", "gene", "peak", "peak_id", "motif_id", "distance", "cor",
  "pvalue", "FDR", "PeakGeneCor", "PeakGeneFDR", "PeakGeneStrength",
  "PeakAccessibility", "PeakNorm", "Peak_CTS", "Target_CTS",
  "CellEdgeScore", "EdgeQuantile", "chr", "start", "end", "strand",
  "type", "width", "gene_start", "gene_end", "tss", "peak_center",
  "variant_id", "sample_id", "ref", "alt", "func", "exonic_func",
  "gene_anno", "impact_class", "n_samples", "max_variant_samples",
  "n_mutated_samples", "VariantMutationRelevance", "mutation_peak_raw",
  "mutation_peak_raw_capped", "MutationPeakScore", "MGRS_edge",
  "Target_MGRS", "CTS_raw", "MGRS_raw", "CTS_z", "MGRS_z", "MRS_raw",
  "Penalty", "CTS_mean_celltype", "CTS_sd_celltype",
  "MGRS_mean_celltype", "MGRS_sd_celltype", "N_targets", "RegulonID",
  "MeanPeakCTS", "MeanTargetCTS", "SDTargetCTS", "MatchDistance",
  "CTS_sd", "MGRS_sd", "CTS_P", "MGRS_P", "MRS_P", "MRS_FDR",
  "MRS_FDR_Celltype", "Significant", "ConcordantPositive", "Candidate",
  "HighConfidence"
))

# Explicit bindings for functions from recommended packages used unqualified in
# data.table-heavy expressions.
head <- utils::head
capture.output <- utils::capture.output
setNames <- stats::setNames

# Some data.table releases do not resolve ..required reliably after package
# byte compilation. Route in-memory mutation tables through the already-tested
# file-input branch without changing any biological calculation.
.ctmrs_calculate_mutation_relevance_impl <- calculate_mutation_relevance

calculate_mutation_relevance <- function(mutations, peak_ids,
                                         min_recurrent_samples = 1L,
                                         buffer = 0L,
                                         output_dir = NULL) {
  temporary_file <- NULL
  if (!is.character(mutations) || length(mutations) != 1L) {
    temporary_file <- tempfile(fileext = ".tsv")
    data.table::fwrite(mutations, temporary_file, sep = "\t", na = "NA")
    mutations <- temporary_file
    on.exit(unlink(temporary_file), add = TRUE)
  }
  .ctmrs_calculate_mutation_relevance_impl(
    mutations = mutations, peak_ids = peak_ids,
    min_recurrent_samples = min_recurrent_samples,
    buffer = buffer, output_dir = output_dir
  )
}
