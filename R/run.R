#' Run the complete ctMRS pipeline
#'
#' @param rna Seurat RNA object or RDS path.
#' @param atac Seurat ATAC object or RDS path with stored motif matches.
#' @param gene_annotation GRCh38 gene annotation GRanges or RDS path.
#' @param mutations ANNOVAR multianno table or file path.
#' @param output_dir Output directory.
#' @param rna_assay RNA assay name.
#' @param peak_assay ATAC assay name.
#' @param celltype_column Metadata column containing cell types.
#' @param variable_genes Number of variable genes.
#' @param peak_gene_distance Maximum peak-to-TSS distance.
#' @param correlation_chunk_size Peak-gene pairs per chunk.
#' @param min_peak_gene_cor Minimum positive correlation.
#' @param max_peak_gene_fdr Maximum peak-gene BH FDR.
#' @param min_targets_per_tf Minimum targets per eRegulon.
#' @param cell_edge_quantile Within-eRegulon edge quantile retained.
#' @param mutation_peak_buffer Peak extension in bases.
#' @param min_recurrent_samples Minimum exact-variant recurrence.
#' @param theta Target-node contribution parameter.
#' @param top_n Number of strongest target MGRS values averaged.
#' @param perm_n Monte Carlo permutations.
#' @param seed Random seed.
#' @param control_pool_n Maximum matched controls.
#' @param min_control_pool Minimum size-matched control pool.
#' @param nominal_p Nominal MRS threshold.
#' @param fdr_cutoff Global MRS FDR threshold.
#' @return A named list containing intermediate and final objects.
#' @export
run_ctMRS <- function(rna, atac, gene_annotation, mutations, output_dir,
                      rna_assay = "RNA", peak_assay = "ATAC",
                      celltype_column = "celltype", variable_genes = 5000L,
                      peak_gene_distance = 500000L,
                      correlation_chunk_size = 2000L,
                      min_peak_gene_cor = 0.05, max_peak_gene_fdr = 0.05,
                      min_targets_per_tf = 5L, cell_edge_quantile = 0.5,
                      mutation_peak_buffer = 0L, min_recurrent_samples = 1L,
                      theta = 0.5, top_n = 5L, perm_n = 1000L,
                      seed = 1234L, control_pool_n = 50L,
                      min_control_pool = 20L, nominal_p = 0.05,
                      fdr_cutoff = 0.05) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Step 1/10: preparing paired RNA/ATAC inputs")
  input <- prepare_ctMRS_inputs(rna, atac, gene_annotation, rna_assay,
                                peak_assay, celltype_column)
  message("Step 2/10: calculating gene CTS")
  target_scores <- calculate_gene_cts(input$rna_data, input$celltype, output_dir)
  message("Step 3/10: extracting TF-peak motif matches")
  tf_peak <- extract_tf_peak(input$single_cell, peak_assay,
                             rownames(input$rna_data), output_dir)
  message("Step 4/10: calculating full-genome peak-gene links")
  links <- link_peaks_genes(
    input$single_cell, input$rna_data, input$atac_counts,
    input$gene_annotation, rna_assay, peak_assay, variable_genes,
    peak_gene_distance, correlation_chunk_size, min_peak_gene_cor,
    max_peak_gene_fdr, output_dir
  )
  message("Step 5/10: building global GRN")
  grn <- build_global_grn(links$positive, tf_peak, min_targets_per_tf, output_dir)
  message("Step 6/10: building cell-type sub-GRNs")
  sub <- build_celltype_subgrn(grn, target_scores, input$atac_counts,
                               input$celltype, cell_edge_quantile,
                               min_targets_per_tf, output_dir)
  message("Step 7/10: mapping somatic mutations to ATAC peaks")
  mutation <- calculate_mutation_relevance(
    mutations, input$peak_ids, min_recurrent_samples,
    mutation_peak_buffer, output_dir
  )
  overlap_n <- length(intersect(unique(grn$peak_id), unique(mutation$peak_scores$peak_id)))
  message("Mutated peaks entering global GRN = ", overlap_n, "/",
          data.table::uniqueN(grn$peak_id), " GRN peaks")
  message("Step 8/10: calculating mutation-weighted target MGRS")
  mgrs <- calculate_target_mgrs(sub$celltype_subgrn, mutation$peak_scores, output_dir)
  message("Step 9/10: aggregating eRegulon CTS, MGRS and MRS")
  observed <- score_eregulons(sub$celltype_subgrn,
                              mgrs$target_mutation_status,
                              theta, top_n, output_dir)
  message("Step 10/10: matched-control Monte Carlo testing")
  tested <- test_eregulons_mc(
    observed$target_info, observed$score_table, perm_n, theta, top_n,
    seed, control_pool_n, min_control_pool, nominal_p, fdr_cutoff,
    output_dir
  )
  result <- list(
    settings = list(
      variable_genes = variable_genes, peak_gene_distance = peak_gene_distance,
      min_peak_gene_cor = min_peak_gene_cor, max_peak_gene_fdr = max_peak_gene_fdr,
      min_targets_per_tf = min_targets_per_tf, cell_edge_quantile = cell_edge_quantile,
      mutation_peak_buffer = mutation_peak_buffer,
      min_recurrent_samples = min_recurrent_samples, theta = theta,
      top_n = top_n, perm_n = perm_n, seed = seed
    ),
    grn = grn, peak_accessibility = sub$peak_accessibility,
    celltype_subgrn = sub$celltype_subgrn,
    mutation_to_peak = mutation$mutation_to_peak,
    peak_scores = mutation$peak_scores,
    mutation_edges = mgrs$mutation_edges,
    target_mutation_status = mgrs$target_mutation_status,
    score_table = observed$score_table,
    monte_carlo = tested$monte_carlo,
    final_results = tested$final_results
  )
  saveRDS(result, file.path(output_dir, "ctMRS_complete_objects.rds"))
  summary <- data.table::data.table(
    Metric = c("paired_cells", "ATAC_peaks", "global_GRN_edges", "global_GRN_peaks",
               "celltype_subGRN_edges", "mutated_peaks", "mutated_GRN_peaks",
               "eRegulons", "significant_MRS", "high_confidence_MRS"),
    Value = c(ncol(input$rna_data), nrow(input$atac_counts), nrow(grn),
              data.table::uniqueN(grn$peak_id), nrow(sub$celltype_subgrn),
              nrow(mutation$peak_scores), overlap_n, nrow(tested$final_results),
              sum(tested$final_results$Significant, na.rm = TRUE),
              sum(tested$final_results$HighConfidence, na.rm = TRUE))
  )
  data.table::fwrite(summary, file.path(output_dir, "00_run_summary.tsv"), sep = "\t")
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "00_sessionInfo.txt"))
  message("ctMRS completed. Results: ", output_dir)
  result
}
