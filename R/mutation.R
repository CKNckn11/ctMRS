#' Calculate cohort-recurrence mutation relevance for ATAC peaks
#'
#' Mutation relevance uses only exact-variant cohort recurrence:
#' `log(1 + n_samples)`. Functional classes are descriptive and are never
#' manually weighted.
#'
#' @param mutations Path to an ANNOVAR multianno TSV or an equivalent data frame.
#' @param peak_ids ATAC peak identifiers such as `chr1-100-200`.
#' @param min_recurrent_samples Minimum samples carrying an exact variant.
#'   The recommended primary value is 1.
#' @param buffer Number of bases added on each side of a peak.
#' @param output_dir Optional output directory.
#' @return A list with variants, mutation_to_peak, peak_scores and diagnostics.
#' @export
calculate_mutation_relevance <- function(mutations, peak_ids,
                                         min_recurrent_samples = 1L,
                                         buffer = 0L,
                                         output_dir = NULL) {
  required <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
                "Gene.refGene", "ExonicFunc.refGene", "Tumor_Sample_Barcode")
  if (is.character(mutations) && length(mutations) == 1L) {
    if (!file.exists(mutations)) stop("Mutation file not found: ", mutations)
    header <- names(data.table::fread(mutations, nrows = 0L, sep = "\t"))
    missing <- setdiff(required, header)
    if (length(missing)) stop("Mutation table lacks columns: ", paste(missing, collapse = ", "))
    mut <- data.table::fread(mutations, sep = "\t", select = required,
                             na.strings = c("", ".", "NA"), showProgress = TRUE)
  } else {
    mut <- .ctmrs_as_dt(mutations)
    .ctmrs_require_columns(mut, required, "mutations")
    mut <- mut[, ..required]
  }
  data.table::setnames(mut, required, c(
    "chr", "start", "end", "ref", "alt", "func",
    "gene_anno", "exonic_func", "sample_id"
  ))
  mut[, chr := as.character(chr)]
  mut[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  mut[chr == "chr23", chr := "chrX"]
  mut[chr == "chr24", chr := "chrY"]
  mut[chr %in% c("chrM", "chrMT"), chr := "chrM"]
  mut[, `:=`(start = as.integer(start), end = as.integer(end))]
  mut <- mut[grepl("^chr(?:[0-9]+|X|Y|M)$", chr) & is.finite(start) &
               is.finite(end) & start >= 1L & end >= start &
               !is.na(sample_id) & nzchar(as.character(sample_id))]
  mut[, variant_id := paste(chr, start, ref, alt, sep = ":")]
  mut <- unique(mut, by = c("variant_id", "sample_id"))
  mut[, impact_class := "other"]
  mut[grepl("intergenic", tolower(data.table::fifelse(is.na(func), "", func))), impact_class := "intergenic"]
  mut[grepl("intronic", tolower(data.table::fifelse(is.na(func), "", func))), impact_class := "intronic"]
  mut[grepl("utr|upstream|downstream", tolower(data.table::fifelse(is.na(func), "", func))), impact_class := "regulatory"]
  mut[grepl("splic", tolower(data.table::fifelse(is.na(func), "", func))) |
        grepl("splic", tolower(data.table::fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "splice"]
  mut[grepl("synonymous", tolower(data.table::fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "synonymous"]
  mut[grepl("nonsynonymous|nonframeshift", tolower(data.table::fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "missense"]
  mut[grepl("stopgain|stoploss|frameshift", tolower(data.table::fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "high_impact"]
  variants <- mut[, .(
    chr = chr[1L], start = start[1L], end = end[1L], ref = ref[1L], alt = alt[1L],
    n_samples = data.table::uniqueN(sample_id), impact_class = impact_class[1L],
    func = func[1L], exonic_func = exonic_func[1L], gene_anno = gene_anno[1L]
  ), by = variant_id]
  n_before <- nrow(variants)
  variants <- variants[n_samples >= as.integer(min_recurrent_samples)]
  if (!nrow(variants)) stop("No variant remains after recurrence filtering.")
  variants[, VariantMutationRelevance := log1p(n_samples)]
  match_peak <- stringr::str_match(
    peak_ids, "^(chr(?:[0-9]+|X|Y|M|MT))[:-]([0-9]+)[:-]([0-9]+)$"
  )
  valid <- stats::complete.cases(match_peak[, 2:4, drop = FALSE])
  if (!any(valid)) stop("No valid ATAC peak ID was detected.")
  peak_gr <- GenomicRanges::GRanges(
    seqnames = match_peak[valid, 2],
    ranges = IRanges::IRanges(as.integer(match_peak[valid, 3]), as.integer(match_peak[valid, 4])),
    peak_id = peak_ids[valid]
  )
  query_gr <- peak_gr
  GenomicRanges::start(query_gr) <- pmax(1L, GenomicRanges::start(query_gr) - as.integer(buffer))
  GenomicRanges::end(query_gr) <- GenomicRanges::end(query_gr) + as.integer(buffer)
  mutation_gr <- GenomicRanges::GRanges(
    seqnames = variants$chr,
    ranges = IRanges::IRanges(variants$start, variants$end),
    variant_id = variants$variant_id
  )
  hits <- GenomicRanges::findOverlaps(mutation_gr, query_gr, ignore.strand = TRUE)
  if (!length(hits)) stop("No mutation overlaps an ATAC peak; check genome build and chromosome names.")
  mapping <- data.table::data.table(
    variant_id = S4Vectors::mcols(mutation_gr)$variant_id[S4Vectors::queryHits(hits)],
    peak_id = S4Vectors::mcols(query_gr)$peak_id[S4Vectors::subjectHits(hits)]
  )
  mapping <- unique(merge(mapping, variants, by = "variant_id", all.x = TRUE),
                    by = c("variant_id", "peak_id"))
  scores <- mapping[, .(
    n_variants = data.table::uniqueN(variant_id),
    max_variant_samples = max(n_samples, na.rm = TRUE),
    mutation_peak_raw = sum(VariantMutationRelevance, na.rm = TRUE)
  ), by = peak_id]
  scores[, n_mutated_samples := max_variant_samples]
  cap <- as.numeric(stats::quantile(scores$mutation_peak_raw, 0.99, na.rm = TRUE))
  if (!is.finite(cap)) cap <- max(scores$mutation_peak_raw, na.rm = TRUE)
  scores[, mutation_peak_raw_capped := pmin(mutation_peak_raw, cap)]
  denominator <- log1p(max(scores$mutation_peak_raw_capped, na.rm = TRUE))
  if (is.finite(denominator) && denominator > 0) {
    scores[, MutationPeakScore := log1p(mutation_peak_raw_capped) / denominator]
  } else {
    scores[, MutationPeakScore := as.numeric(mutation_peak_raw_capped > 0)]
  }
  scores[, `:=`(
    MutationRelevanceRaw = mutation_peak_raw,
    MutationRelevanceCapped = mutation_peak_raw_capped,
    MutationRelevance = MutationPeakScore
  )]
  .ctmrs_write(mapping, output_dir, "07_mutation_peak.tsv.gz")
  .ctmrs_write(scores, output_dir, "08_peak_mutation_scores.tsv")
  diagnostics <- data.table::data.table(
    Metric = c("variants_before_filter", "variants_retained", "variants_overlapping_peaks", "mutated_peaks"),
    Value = c(n_before, nrow(variants), data.table::uniqueN(mapping$variant_id), nrow(scores))
  )
  list(variants = variants, mutation_to_peak = mapping,
       peak_scores = scores, diagnostics = diagnostics)
}

#' Calculate mutation-weighted regulatory edges and target MGRS
#' @param celltype_subgrn Cell-type sub-GRN with peak_id and CellEdgeScore.
#' @param peak_scores Peak mutation scores.
#' @param output_dir Optional output directory.
#' @return A list containing mutation_edges, target_mutation_status and diagnostics.
#' @export
calculate_target_mgrs <- function(celltype_subgrn, peak_scores, output_dir = NULL) {
  subgrn <- .ctmrs_as_dt(celltype_subgrn); scores <- .ctmrs_as_dt(peak_scores)
  .ctmrs_require_columns(subgrn, c("TF", "Celltype", "gene", "peak_id", "CellEdgeScore"), "celltype_subgrn")
  .ctmrs_require_columns(scores, c("peak_id", "MutationPeakScore"), "peak_scores")
  if (anyDuplicated(scores$peak_id)) stop("peak_scores contains duplicated peak_id values.")
  if (any(!is.finite(subgrn$CellEdgeScore) | subgrn$CellEdgeScore < 0)) {
    stop("CellEdgeScore contains invalid or negative values.")
  }
  edges <- merge(subgrn, scores, by = "peak_id", all.x = TRUE)
  if (nrow(edges) != nrow(subgrn)) stop("Mutation merge changed the number of sub-GRN edges.")
  for (column in c("MutationPeakScore", "mutation_peak_raw", "mutation_peak_raw_capped")) {
    if (!column %in% names(edges)) edges[, (column) := 0]
    edges[is.na(get(column)), (column) := 0]
  }
  for (column in c("n_variants", "n_mutated_samples", "max_variant_samples")) {
    if (!column %in% names(edges)) edges[, (column) := 0L]
    edges[is.na(get(column)), (column) := 0L]
  }
  edges[, MGRS_edge := MutationPeakScore * CellEdgeScore]
  target <- edges[, .(Target_MGRS = max(MGRS_edge, na.rm = TRUE)), by = .(TF, Celltype, gene)]
  target[!is.finite(Target_MGRS), Target_MGRS := 0]
  diag <- data.table::data.table(
    Metric = c("subgrn_edges", "mutation_linked_edges", "subgrn_peaks",
               "mutation_linked_subgrn_peaks", "positive_target_mgrs"),
    Value = c(nrow(edges), sum(edges$MutationPeakScore > 0),
              data.table::uniqueN(edges$peak_id),
              data.table::uniqueN(edges[MutationPeakScore > 0, peak_id]),
              sum(target$Target_MGRS > 0))
  )
  .ctmrs_write(edges, output_dir, "09_TF_peak_gene_mutation_edges.tsv.gz")
  .ctmrs_write(target, output_dir, "10_TF_celltype_target_MGRS.tsv")
  list(mutation_edges = edges, target_mutation_status = target, diagnostics = diag)
}
