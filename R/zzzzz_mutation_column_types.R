# Normalize mutation annotation columns before the core parser. This also
# handles small cohorts where an annotation column is entirely missing/NA and
# would otherwise be inferred as logical by fread.
.ctmrs_calculate_mutation_relevance_typed_impl <- calculate_mutation_relevance

calculate_mutation_relevance <- function(mutations, peak_ids,
                                         min_recurrent_samples = 1L,
                                         buffer = 0L,
                                         output_dir = NULL) {
  required <- c(
    "Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
    "Gene.refGene", "ExonicFunc.refGene", "Tumor_Sample_Barcode"
  )
  if (is.character(mutations) && length(mutations) == 1L) {
    if (!file.exists(mutations)) stop("Mutation file not found: ", mutations)
    mutations <- data.table::fread(
      mutations, sep = "\t", select = required,
      na.strings = c("", ".", "NA"), showProgress = TRUE
    )
  } else {
    mutations <- data.table::as.data.table(data.table::copy(mutations))
  }
  .ctmrs_require_columns(mutations, required, "mutations")
  character_columns <- c(
    "Chr", "Ref", "Alt", "Func.refGene", "Gene.refGene",
    "ExonicFunc.refGene", "Tumor_Sample_Barcode"
  )
  for (column in character_columns) {
    data.table::set(mutations, j = column, value = as.character(mutations[[column]]))
  }
  for (column in c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene")) {
    values <- mutations[[column]]
    values[is.na(values) | !nzchar(values)] <- "unknown"
    data.table::set(mutations, j = column, value = values)
  }
  .ctmrs_calculate_mutation_relevance_typed_impl(
    mutations = mutations, peak_ids = peak_ids,
    min_recurrent_samples = min_recurrent_samples,
    buffer = buffer, output_dir = output_dir
  )
}
