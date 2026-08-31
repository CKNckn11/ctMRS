#' Prepare paired RNA and ATAC inputs
#'
#' @param rna A Seurat object or path to an RNA RDS file.
#' @param atac A Seurat object or path to an ATAC RDS file.
#' @param gene_annotation A GRanges object or path to its RDS file.
#' @param rna_assay RNA assay name.
#' @param peak_assay ChromatinAssay name.
#' @param celltype_column RNA metadata column containing cell types.
#' @return A list containing the combined Seurat object, matrices, cell types,
#'   peak identifiers and gene annotation.
#' @export
prepare_ctMRS_inputs <- function(rna, atac, gene_annotation,
                                 rna_assay = "RNA", peak_assay = "ATAC",
                                 celltype_column = "celltype") {
  rna <- .ctmrs_read(rna, "RNA")
  atac <- .ctmrs_read(atac, "ATAC")
  gene_annotation <- .ctmrs_read(gene_annotation, "Gene annotation")
  if (!rna_assay %in% SeuratObject::Assays(rna)) stop("RNA assay not found: ", rna_assay)
  if (!peak_assay %in% SeuratObject::Assays(atac)) stop("ATAC assay not found: ", peak_assay)
  if (!inherits(atac[[peak_assay]], "ChromatinAssay")) {
    stop("ATAC assay is not a Signac ChromatinAssay.")
  }
  if (!celltype_column %in% colnames(rna[[]])) {
    stop("RNA metadata lacks cell-type column: ", celltype_column)
  }
  common_cells <- intersect(colnames(rna), colnames(atac))
  if (!length(common_cells)) stop("RNA and ATAC have no shared cell barcodes.")
  rna <- rna[, common_cells]
  atac <- atac[, common_cells]
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    rna <- SeuratObject::JoinLayers(rna, assay = rna_assay)
  }
  single_cell <- atac
  single_cell[[rna_assay]] <- rna[[rna_assay]]
  single_cell[[celltype_column]] <- rna[[celltype_column, drop = TRUE]]
  celltype <- as.character(single_cell[[celltype_column, drop = TRUE]])
  if (anyNA(celltype) || any(!nzchar(celltype))) stop("Cell-type annotations contain NA or empty values.")
  rna_data <- SeuratObject::LayerData(single_cell, assay = rna_assay, layer = "data")
  atac_counts <- SeuratObject::LayerData(single_cell, assay = peak_assay, layer = "counts")
  if (!identical(colnames(rna_data), colnames(atac_counts))) {
    atac_counts <- atac_counts[, colnames(rna_data), drop = FALSE]
  }
  if (!identical(colnames(rna_data), colnames(atac_counts))) stop("RNA and ATAC cell order differs.")
  peak_ids <- rownames(atac_counts)
  if (is.null(peak_ids) || anyDuplicated(peak_ids)) stop("ATAC peak IDs are missing or duplicated.")
  list(single_cell = single_cell, rna_data = rna_data,
       atac_counts = atac_counts, celltype = celltype,
       peak_ids = peak_ids, gene_annotation = gene_annotation,
       rna_assay = rna_assay, peak_assay = peak_assay)
}

#' Calculate gene cell-type specificity by cosine normalization
#' @param rna_data Gene-by-cell normalized expression matrix.
#' @param celltype Cell-type label for every matrix column.
#' @param output_dir Optional output directory.
#' @return A data.table with gene, Celltype and Gene_CTS.
#' @export
calculate_gene_cts <- function(rna_data, celltype, output_dir = NULL) {
  if (ncol(rna_data) != length(celltype)) stop("celltype length must equal the number of cells.")
  groups <- factor(celltype)
  membership <- Matrix::sparseMatrix(
    i = seq_along(groups), j = as.integer(groups), x = 1,
    dims = c(length(groups), nlevels(groups)),
    dimnames = list(colnames(rna_data), levels(groups))
  )
  average <- rna_data %*% membership
  average <- sweep(as.matrix(average), 2L, as.numeric(table(groups)), "/")
  denominator <- sqrt(rowSums(average^2))
  denominator[!is.finite(denominator) | denominator == 0] <- 1
  cts <- average / denominator
  ans <- data.table::data.table(
    gene = rep(rownames(cts), times = ncol(cts)),
    Celltype = rep(colnames(cts), each = nrow(cts)),
    Gene_CTS = as.vector(cts)
  )
  .ctmrs_write(ans, output_dir, "00_gene_celltype_CTS.tsv")
  ans
}

#' Extract TF--peak motif matches from a Signac Motif object
#' @param single_cell Combined Seurat object.
#' @param peak_assay ChromatinAssay name.
#' @param rna_genes RNA feature names used to retain expressed TFs.
#' @param output_dir Optional output directory.
#' @return A data.table with motif_id, peak_id and TF.
#' @export
extract_tf_peak <- function(single_cell, peak_assay = "ATAC", rna_genes,
                            output_dir = NULL) {
  motif_obj <- Signac::Motifs(single_cell[[peak_assay]])
  if (is.null(motif_obj)) stop("ATAC assay has no stored Motif object.")
  motif_matrix <- motif_obj@data
  motif_names <- motif_obj@motif.names
  if (!length(motif_names) || is.null(names(motif_names))) stop("Motif object lacks motif-to-TF names.")
  motif_map <- data.table::rbindlist(lapply(names(motif_names), function(id) {
    data.table::data.table(motif_id = id, TF_raw = as.character(motif_names[[id]]))
  }))
  motif_map <- motif_map[, .(TF = unlist(strsplit(TF_raw, "::", fixed = TRUE))), by = motif_id]
  motif_map[, TF := trimws(sub("\\(.*$", "", TF))]
  lookup <- setNames(rna_genes, toupper(rna_genes))
  motif_map[, TF := unname(lookup[toupper(TF)])]
  motif_map <- unique(motif_map[!is.na(TF) & nzchar(TF), .(motif_id, TF)])
  sm <- Matrix::summary(motif_matrix)
  tf_peak <- data.table::data.table(
    motif_id = colnames(motif_matrix)[sm$j],
    peak_id = rownames(motif_matrix)[sm$i]
  )
  tf_peak <- unique(merge(tf_peak, motif_map, by = "motif_id", allow.cartesian = TRUE))
  .ctmrs_write(tf_peak, output_dir, "01_TF_peak.tsv.gz")
  tf_peak
}

#' Calculate full-genome peak--gene correlations
#' @param single_cell Combined Seurat object.
#' @param rna_data Gene-by-cell normalized RNA matrix.
#' @param atac_counts Peak-by-cell ATAC counts.
#' @param gene_annotation GRCh38 gene annotation as GRanges.
#' @param rna_assay RNA assay name.
#' @param peak_assay ATAC assay name.
#' @param variable_genes Number of variable genes.
#' @param distance Maximum peak-to-TSS distance.
#' @param chunk_size Number of pairs evaluated per chunk.
#' @param min_cor Minimum positive correlation.
#' @param max_fdr Maximum BH FDR.
#' @param output_dir Optional output directory.
#' @return A list containing all, significant and positive links.
#' @export
link_peaks_genes <- function(single_cell, rna_data, atac_counts, gene_annotation,
                             rna_assay = "RNA", peak_assay = "ATAC",
                             variable_genes = 5000L, distance = 500000L,
                             chunk_size = 2000L, min_cor = 0.05,
                             max_fdr = 0.05, output_dir = NULL) {
  SeuratObject::DefaultAssay(single_cell) <- rna_assay
  single_cell <- Seurat::FindVariableFeatures(
    single_cell, assay = rna_assay, selection.method = "vst",
    nfeatures = as.integer(variable_genes), verbose = FALSE
  )
  genes_use <- intersect(SeuratObject::VariableFeatures(single_cell[[rna_assay]]), rownames(rna_data))
  if (!length(genes_use)) stop("No variable gene is present in RNA data.")
  depth <- Matrix::colSums(atac_counts)
  if (any(depth <= 0)) stop("ATAC matrix contains cells with zero depth.")
  atac_norm <- t(t(atac_counts) / depth) * 1e4
  atac_norm@x <- log1p(atac_norm@x)
  peak_ids <- rownames(atac_counts)
  peak_gr <- Signac::granges(single_cell[[peak_assay]])
  if (length(peak_gr) != nrow(atac_norm)) stop("Peak ranges and matrix differ in length.")
  m <- S4Vectors::mcols(gene_annotation)
  anno <- data.table::data.table(
    chr = as.character(GenomicRanges::seqnames(gene_annotation)),
    start = GenomicRanges::start(gene_annotation), end = GenomicRanges::end(gene_annotation),
    strand = as.character(GenomicRanges::strand(gene_annotation)),
    gene = as.character(m$gene_name), type = as.character(m$type)
  )
  anno <- anno[type == "exon" & !is.na(gene) & nzchar(gene) & gene %in% genes_use & strand %in% c("+", "-")]
  if (!nrow(anno)) stop("No variable gene has a valid exon annotation.")
  genes <- anno[, .(gene_start = min(start), gene_end = max(end)), by = .(gene, chr, strand)]
  genes[, width := gene_end - gene_start + 1L]
  data.table::setorder(genes, gene, -width)
  genes <- genes[!duplicated(gene)]
  genes[, `:=`(width = NULL, tss = ifelse(strand == "+", gene_start, gene_end))]
  windows <- genes[, .(chr, start = pmax(1L, tss - distance), end = tss + distance, gene, tss)]
  peaks <- data.table::data.table(
    peak_idx = seq_along(peak_gr), peak = peak_ids,
    chr = as.character(GenomicRanges::seqnames(peak_gr)),
    start = GenomicRanges::start(peak_gr), end = GenomicRanges::end(peak_gr)
  )
  peaks[, peak_center := as.integer(round((start + end) / 2))]
  data.table::setkey(windows, chr, start, end)
  pairs <- data.table::foverlaps(peaks, windows,
    by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"),
    type = "any", nomatch = 0L)
  pairs[, distance := peak_center - tss]
  pairs <- unique(pairs[, .(peak, gene, distance)])
  pairs <- pairs[peak %in% rownames(atac_norm) & gene %in% rownames(rna_data)]
  if (!nrow(pairs)) stop("No candidate peak-gene pair was found.")
  chunks <- split(seq_len(nrow(pairs)), ceiling(seq_len(nrow(pairs)) / as.integer(chunk_size)))
  n_cells <- ncol(rna_data)
  results <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    p <- pairs[chunks[[i]]]
    a <- atac_norm[match(p$peak, rownames(atac_norm)), , drop = FALSE]
    r <- rna_data[match(p$gene, rownames(rna_data)), , drop = FALSE]
    a2 <- a; a2@x <- a2@x^2
    r2 <- r; r2@x <- r2@x^2
    sa <- Matrix::rowSums(a); sr <- Matrix::rowSums(r)
    cov_num <- Matrix::rowSums(a * r) - sa * sr / n_cells
    ass <- pmax(Matrix::rowSums(a2) - sa^2 / n_cells, 0)
    rss <- pmax(Matrix::rowSums(r2) - sr^2 / n_cells, 0)
    correlation <- cov_num / sqrt(ass * rss)
    correlation[!is.finite(correlation)] <- NA_real_
    correlation <- pmax(pmin(correlation, 0.999999), -0.999999)
    t_stat <- correlation * sqrt((n_cells - 2) / pmax(1 - correlation^2, .Machine$double.eps))
    results[[i]] <- data.table::data.table(
      peak = p$peak, gene = p$gene, distance = p$distance,
      cor = correlation, pvalue = 2 * stats::pt(-abs(t_stat), df = n_cells - 2)
    )[is.finite(cor) & is.finite(pvalue)]
    if (i %% 10L == 0L || i == length(chunks)) message("Peak-gene chunk ", i, "/", length(chunks))
  }
  all_links <- data.table::rbindlist(results)
  all_links[, FDR := stats::p.adjust(pvalue, method = "BH")]
  significant <- all_links[abs(cor) >= min_cor & FDR <= max_fdr]
  positive <- all_links[cor >= min_cor & FDR <= max_fdr]
  if (!nrow(positive)) stop("No positive peak-gene link passed thresholds.")
  .ctmrs_write(all_links, output_dir, "02_peak_gene_all.tsv.gz")
  .ctmrs_write(significant, output_dir, "02_peak_gene_significant.tsv")
  .ctmrs_write(positive, output_dir, "02_peak_gene_positive.tsv")
  list(all = all_links, significant = significant, positive = positive)
}

#' Construct a TF--peak--gene global GRN
#' @param peak_gene_links Positive significant peak-gene links.
#' @param tf_peak TF-peak motif table.
#' @param min_targets Minimum distinct targets per TF.
#' @param output_dir Optional output directory.
#' @return Global GRN data.table using the canonical peak_id schema.
#' @export
build_global_grn <- function(peak_gene_links, tf_peak, min_targets = 5L,
                             output_dir = NULL) {
  links <- .ctmrs_as_dt(peak_gene_links)
  tf_peak <- .ctmrs_as_dt(tf_peak)
  .ctmrs_require_columns(links, c("peak", "gene", "distance", "cor", "pvalue", "FDR"), "peak_gene_links")
  .ctmrs_require_columns(tf_peak, c("peak_id", "motif_id", "TF"), "tf_peak")
  raw <- merge(links, tf_peak, by.x = "peak", by.y = "peak_id", allow.cartesian = TRUE)
  raw <- raw[TF != gene]
  grn <- raw[, .(
    motif_id = paste(sort(unique(motif_id)), collapse = ";"),
    distance = distance[1L], PeakGeneCor = cor[1L], PeakGeneP = pvalue[1L],
    PeakGeneFDR = FDR[1L]
  ), by = .(TF, gene, peak_id = peak)]
  grn[, PeakGeneStrength := PeakGeneCor *
        (pmin(-log10(pmax(PeakGeneFDR, .Machine$double.xmin)), 50) / 50)]
  counts <- unique(grn[, .(TF, gene)])[, .N, by = TF]
  grn <- grn[TF %in% counts[N >= min_targets, TF]]
  if (!nrow(grn)) stop("No TF retained the configured number of targets.")
  data.table::setcolorder(grn, c("TF", "gene", "peak_id", "motif_id", "distance",
                                 "PeakGeneCor", "PeakGeneP", "PeakGeneFDR", "PeakGeneStrength"))
  .ctmrs_write(grn, output_dir, "03_TF_peak_gene_global_GRN.tsv")
  grn
}

#' Build cell-type-specific sub-GRNs
#' @param grn Global GRN from `build_global_grn()`.
#' @param target_scores Gene CTS table.
#' @param atac_counts Peak-by-cell count matrix.
#' @param celltype Cell-type labels.
#' @param edge_quantile Within-eRegulon edge quantile retained.
#' @param min_targets Minimum targets in each TF-cell-type eRegulon.
#' @param output_dir Optional output directory.
#' @return A list containing accessibility, complete cell-type GRN and sub-GRN.
#' @export
build_celltype_subgrn <- function(grn, target_scores, atac_counts, celltype,
                                  edge_quantile = 0.5, min_targets = 5L,
                                  output_dir = NULL) {
  .ctmrs_validate_probability(edge_quantile, "edge_quantile")
  grn <- .ctmrs_as_dt(grn); target_scores <- .ctmrs_as_dt(target_scores)
  groups <- factor(celltype)
  access <- data.table::rbindlist(lapply(levels(groups), function(ct) {
    cells <- colnames(atac_counts)[celltype == ct]
    data.table::data.table(peak_id = rownames(atac_counts), Celltype = ct,
      PeakAccessibility = as.numeric(Matrix::rowMeans(atac_counts[, cells, drop = FALSE] > 0)))
  }))
  access[, PeakNorm := sqrt(sum(PeakAccessibility^2)), by = peak_id]
  access[!is.finite(PeakNorm) | PeakNorm == 0, PeakNorm := 1]
  access[, `:=`(Peak_CTS = PeakAccessibility / PeakNorm, PeakNorm = NULL)]
  ctgrn <- merge(grn, target_scores, by = "gene", allow.cartesian = TRUE)
  data.table::setnames(ctgrn, "Gene_CTS", "Target_CTS")
  tf_cts <- data.table::copy(target_scores)
  data.table::setnames(tf_cts, c("gene", "Gene_CTS"), c("TF", "TF_CTS"))
  ctgrn <- merge(ctgrn, tf_cts, by = c("TF", "Celltype"), all.x = TRUE)
  ctgrn[is.na(TF_CTS), TF_CTS := 0]
  ctgrn <- merge(ctgrn, access, by = c("peak_id", "Celltype"), all.x = TRUE)
  ctgrn[is.na(PeakAccessibility), PeakAccessibility := 0]
  ctgrn[is.na(Peak_CTS), Peak_CTS := 0]
  ctgrn[, CellEdgeScore := PeakGeneStrength * Peak_CTS *
          sqrt(pmax(TF_CTS, 0) * pmax(Target_CTS, 0))]
  ctgrn <- ctgrn[is.finite(CellEdgeScore) & CellEdgeScore > 0]
  ctgrn[, EdgeQuantile := stats::quantile(CellEdgeScore, edge_quantile, na.rm = TRUE),
        by = .(TF, Celltype)]
  subgrn <- ctgrn[CellEdgeScore >= EdgeQuantile]
  valid <- unique(subgrn[, .(TF, Celltype, gene)])[, .N, by = .(TF, Celltype)][N >= min_targets]
  subgrn <- merge(subgrn, valid[, .(TF, Celltype)], by = c("TF", "Celltype"))
  .ctmrs_write(access, output_dir, "04_peak_celltype_accessibility.tsv.gz")
  .ctmrs_write(ctgrn, output_dir, "05_celltype_GRN.tsv.gz")
  .ctmrs_write(subgrn, output_dir, "06_celltype_subGRN.tsv.gz")
  list(peak_accessibility = access, celltype_grn = ctgrn, celltype_subgrn = subgrn)
}
