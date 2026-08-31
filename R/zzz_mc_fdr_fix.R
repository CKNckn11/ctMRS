# Keep global and cell-type-specific multiple-testing corrections distinct.
.ctmrs_test_eregulons_mc_impl <- test_eregulons_mc

test_eregulons_mc <- function(target_info, score_table, perm_n = 1000L,
                              theta = 0.5, top_n = 5L, seed = 1234L,
                              control_pool_n = 50L, min_control_pool = 20L,
                              nominal_p = 0.05, fdr_cutoff = 0.05,
                              output_dir = NULL) {
  ans <- .ctmrs_test_eregulons_mc_impl(
    target_info = target_info, score_table = score_table,
    perm_n = perm_n, theta = theta, top_n = top_n, seed = seed,
    control_pool_n = control_pool_n, min_control_pool = min_control_pool,
    nominal_p = nominal_p, fdr_cutoff = fdr_cutoff,
    output_dir = NULL
  )
  final <- data.table::as.data.table(data.table::copy(ans$final_results))
  final[, `:=`(
    CTS_FDR = stats::p.adjust(CTS_P, method = "BH"),
    MGRS_FDR = stats::p.adjust(MGRS_P, method = "BH"),
    MRS_FDR = stats::p.adjust(MRS_P, method = "BH")
  )]
  final[, MRS_FDR_Celltype := stats::p.adjust(MRS_P, method = "BH"), by = Celltype]
  final[, Significant := !is.na(MRS_P) & !is.na(MRS_FDR) &
          MRS_P < nominal_p & MRS_FDR < fdr_cutoff]
  final[, HighConfidence := Significant & ConcordantPositive]
  data.table::setorder(final, MRS_FDR, MRS_P, -MRS_raw)
  ans$final_results <- final
  if (!is.null(output_dir)) {
    .ctmrs_write(ans$monte_carlo, output_dir, "12_matched_MonteCarlo_statistics.tsv")
    .ctmrs_write(final, output_dir, "12_final_results.tsv")
    .ctmrs_write(final[HighConfidence == TRUE], output_dir, "12_final_high_confidence.tsv")
  }
  ans
}
