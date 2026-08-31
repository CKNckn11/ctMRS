#' Aggregate eRegulon CTS, MGRS and MRS
#' @param celltype_subgrn Cell-type-specific sub-GRN.
#' @param target_mutation_status Target-level MGRS table.
#' @param theta Target-node contribution parameter.
#' @param top_n Number of strongest target MGRS values averaged.
#' @param output_dir Optional output directory.
#' @return A list with target_info and score_table.
#' @export
score_eregulons <- function(celltype_subgrn, target_mutation_status,
                            theta = 0.5, top_n = 5L, output_dir = NULL) {
  subgrn <- .ctmrs_as_dt(celltype_subgrn)
  status <- .ctmrs_as_dt(target_mutation_status)
  info <- subgrn[, .(
    TF_CTS = TF_CTS[1L], Target_CTS = max(Target_CTS, na.rm = TRUE),
    MeanPeakCTS = mean(Peak_CTS, na.rm = TRUE), MaxPeakCTS = max(Peak_CTS, na.rm = TRUE),
    MaxCellEdgeScore = max(CellEdgeScore, na.rm = TRUE)
  ), by = .(TF, Celltype, gene)]
  info <- merge(info, status, by = c("TF", "Celltype", "gene"), all.x = TRUE)
  info[is.na(Target_MGRS), Target_MGRS := 0]
  scores <- info[, {
    x <- Target_MGRS; x[!is.finite(x)] <- 0
    top <- mean(head(sort(x, decreasing = TRUE), min(as.integer(top_n), length(x))))
    if (!is.finite(top)) top <- 0
    n <- data.table::uniqueN(gene)
    .(N_targets = n,
      CTS_raw = TF_CTS[1L] + theta * sum(Target_CTS, na.rm = TRUE) / sqrt(n),
      N_mutation_linked_targets = sum(x > 0), MutationLinkedFraction = mean(x > 0),
      MGRS_raw = top + theta * sum(x, na.rm = TRUE) / sqrt(n))
  }, by = .(TF, Celltype)]
  scores[, `:=`(
    CTS_mean_celltype = mean(CTS_raw, na.rm = TRUE),
    CTS_sd_celltype = stats::sd(CTS_raw, na.rm = TRUE),
    MGRS_mean_celltype = mean(MGRS_raw, na.rm = TRUE),
    MGRS_sd_celltype = stats::sd(MGRS_raw, na.rm = TRUE)
  ), by = Celltype]
  scores[!is.finite(CTS_sd_celltype) | CTS_sd_celltype == 0, CTS_sd_celltype := 1]
  scores[!is.finite(MGRS_sd_celltype) | MGRS_sd_celltype == 0, MGRS_sd_celltype := 1]
  scores[, `:=`(
    CTS_z = (CTS_raw - CTS_mean_celltype) / CTS_sd_celltype,
    MGRS_z = (MGRS_raw - MGRS_mean_celltype) / MGRS_sd_celltype
  )]
  scores[, Penalty := abs(CTS_z - MGRS_z) / sqrt(2)]
  scores[, MRS_raw := CTS_z + MGRS_z - Penalty]
  scores[, c("CTS_mean_celltype", "CTS_sd_celltype", "MGRS_mean_celltype", "MGRS_sd_celltype") := NULL]
  .ctmrs_write(scores, output_dir, "11_eRegulon_CTS_MGRS_MRS.tsv")
  list(target_info = info, score_table = scores)
}

#' Matched-control Monte Carlo tests for eRegulons
#' @param target_info Target-level CTS and MGRS information.
#' @param score_table Observed eRegulon scores.
#' @param perm_n Number of permutations.
#' @param theta Target-node contribution parameter.
#' @param top_n Number of strongest targets averaged.
#' @param seed Random seed.
#' @param control_pool_n Maximum matched control TFs.
#' @param min_control_pool Minimum size-matched pool before relaxing size matching.
#' @param nominal_p Nominal MRS threshold.
#' @param fdr_cutoff Global MRS FDR threshold.
#' @param output_dir Optional output directory.
#' @return A list with Monte Carlo statistics and final results.
#' @export
test_eregulons_mc <- function(target_info, score_table, perm_n = 1000L,
                              theta = 0.5, top_n = 5L, seed = 1234L,
                              control_pool_n = 50L, min_control_pool = 20L,
                              nominal_p = 0.05, fdr_cutoff = 0.05,
                              output_dir = NULL) {
  info <- .ctmrs_as_dt(target_info); scores <- .ctmrs_as_dt(score_table)
  meta <- info[, .(
    N_targets = data.table::uniqueN(gene), TF_CTS = TF_CTS[1L],
    MeanTargetCTS = mean(Target_CTS, na.rm = TRUE),
    SDTargetCTS = stats::sd(Target_CTS, na.rm = TRUE),
    MeanPeakCTS = mean(MeanPeakCTS, na.rm = TRUE)
  ), by = .(TF, Celltype)]
  meta[!is.finite(SDTargetCTS), SDTargetCTS := 0]
  scales <- scores[, .(
    CTS_mean = mean(CTS_raw), CTS_sd = stats::sd(CTS_raw),
    MGRS_mean = mean(MGRS_raw), MGRS_sd = stats::sd(MGRS_raw)
  ), by = Celltype]
  scales[!is.finite(CTS_sd) | CTS_sd == 0, CTS_sd := 1]
  scales[!is.finite(MGRS_sd) | MGRS_sd == 0, MGRS_sd := 1]
  info[, RegulonID := paste(TF, Celltype, sep = "___")]
  targets <- split(info, info$RegulonID)
  set.seed(seed)
  out <- vector("list", nrow(scores))
  for (i in seq_len(nrow(scores))) {
    tf <- scores$TF[i]; ct <- scores$Celltype[i]; n <- scores$N_targets[i]
    observed <- scores[i]
    obs_meta <- meta[TF == tf & Celltype == ct]
    candidates <- data.table::copy(meta[Celltype == ct & TF != tf])
    if (nrow(obs_meta) != 1L || !nrow(candidates)) next
    size_pool <- candidates[N_targets >= max(5L, floor(n * 0.7)) & N_targets <= ceiling(n * 1.3)]
    if (nrow(size_pool) >= min_control_pool) candidates <- size_pool
    sds <- c(stats::sd(log1p(candidates$N_targets)), stats::sd(candidates$TF_CTS),
             stats::sd(candidates$MeanTargetCTS), stats::sd(candidates$MeanPeakCTS))
    sds[!is.finite(sds) | sds == 0] <- 1
    candidates[, MatchDistance := sqrt(
      ((log1p(N_targets) - log1p(obs_meta$N_targets)) / sds[1])^2 +
      ((TF_CTS - obs_meta$TF_CTS) / sds[2])^2 +
      ((MeanTargetCTS - obs_meta$MeanTargetCTS) / sds[3])^2 +
      ((MeanPeakCTS - obs_meta$MeanPeakCTS) / sds[4])^2)]
    data.table::setorder(candidates, MatchDistance)
    candidates <- head(candidates[is.finite(MatchDistance)], min(control_pool_n, nrow(candidates)))
    if (!nrow(candidates)) next
    scale_now <- scales[Celltype == ct]
    null_cts <- null_mgrs <- null_mrs <- rep(NA_real_, as.integer(perm_n))
    for (j in seq_len(as.integer(perm_n))) {
      random_tf <- sample(candidates$TF, 1L)
      random_targets <- targets[[paste(random_tf, ct, sep = "___")]]
      if (is.null(random_targets) || !nrow(random_targets)) next
      rows <- sample(seq_len(nrow(random_targets)), n, replace = nrow(random_targets) < n)
      rt <- random_targets[rows]
      cts <- rt$TF_CTS[1L] + theta * sum(rt$Target_CTS, na.rm = TRUE) / sqrt(n)
      x <- rt$Target_MGRS; x[!is.finite(x)] <- 0
      top <- mean(head(sort(x, decreasing = TRUE), min(top_n, length(x))))
      if (!is.finite(top)) top <- 0
      mgrs <- top + theta * sum(x) / sqrt(n)
      cz <- (cts - scale_now$CTS_mean) / scale_now$CTS_sd
      mz <- (mgrs - scale_now$MGRS_mean) / scale_now$MGRS_sd
      null_cts[j] <- cts; null_mgrs[j] <- mgrs
      null_mrs[j] <- cz + mz - abs(cz - mz) / sqrt(2)
    }
    valid <- is.finite(null_cts) & is.finite(null_mgrs) & is.finite(null_mrs)
    null_cts <- null_cts[valid]; null_mgrs <- null_mgrs[valid]; null_mrs <- null_mrs[valid]
    nv <- length(null_mrs); if (nv < 10L) next
    zfun <- function(obs, null) if (stats::sd(null) > 0) (obs - mean(null)) / stats::sd(null) else NA_real_
    out[[i]] <- data.table::data.table(
      TF = tf, Celltype = ct, N_controls = nrow(candidates), N_valid_permutations = nv,
      MeanMatchDistance = mean(candidates$MatchDistance), MinMatchDistance = min(candidates$MatchDistance),
      CTS_P = (1 + sum(null_cts >= observed$CTS_raw)) / (1 + nv),
      MGRS_P = (1 + sum(null_mgrs >= observed$MGRS_raw)) / (1 + nv),
      MRS_P = (1 + sum(null_mrs >= observed$MRS_raw)) / (1 + nv),
      CTS_null_z = zfun(observed$CTS_raw, null_cts),
      MGRS_null_z = zfun(observed$MGRS_raw, null_mgrs),
      MRS_null_z = zfun(observed$MRS_raw, null_mrs),
      Null_CTS_mean = mean(null_cts), Null_CTS_sd = stats::sd(null_cts),
      Null_MGRS_mean = mean(null_mgrs), Null_MGRS_sd = stats::sd(null_mgrs),
      Null_MRS_mean = mean(null_mrs), Null_MRS_sd = stats::sd(null_mrs))
    if (i %% 100L == 0L || i == nrow(scores)) message("Monte Carlo eRegulon ", i, "/", nrow(scores))
  }
  mc <- data.table::rbindlist(out, fill = TRUE)
  final <- merge(scores, mc, by = c("TF", "Celltype"), all.x = TRUE)
  final[, `:=`(
    CTS_FDR = stats::p.adjust(CTS_P, "BH"),
    MGRS_FDR = stats::p.adjust(MGRS_P, "BH"),
    MRS_FDR = stats::p.adjust(MRS_P, "BH"),
    MRS_FDR_Celltype = stats::p.adjust(MRS_P, "BH")
  ), by = Celltype]
  final[, `:=`(
    Significant = !is.na(MRS_P) & !is.na(MRS_FDR) & MRS_P < nominal_p & MRS_FDR < fdr_cutoff,
    ConcordantPositive = CTS_z > 0 & MGRS_z > 0,
    Candidate = !is.na(MRS_P) & MRS_P < nominal_p
  )]
  final[, `:=`(HighConfidence = Significant & ConcordantPositive,
               StrongCandidate = Candidate & ConcordantPositive)]
  data.table::setorder(final, MRS_FDR, MRS_P, -MRS_raw)
  .ctmrs_write(mc, output_dir, "12_matched_MonteCarlo_statistics.tsv")
  .ctmrs_write(final, output_dir, "12_final_results.tsv")
  .ctmrs_write(final[HighConfidence == TRUE], output_dir, "12_final_high_confidence.tsv")
  list(monte_carlo = mc, final_results = final)
}
