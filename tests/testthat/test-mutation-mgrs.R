test_that("mutation relevance and MGRS retain peak IDs", {
  mutations <- data.frame(
    Chr = c("1", "1", "1"), Start = c(150L, 150L, 350L),
    End = c(150L, 150L, 350L), Ref = c("A", "A", "G"),
    Alt = c("T", "T", "C"), Func.refGene = c("intronic", "intronic", "intergenic"),
    Gene.refGene = c("A", "A", "B"), ExonicFunc.refGene = NA_character_,
    Tumor_Sample_Barcode = c("S1", "S2", "S1"), check.names = FALSE
  )
  x <- calculate_mutation_relevance(mutations, c("chr1-100-200", "chr1-300-400"))
  expect_equal(nrow(x$peak_scores), 2L)
  expect_equal(x$peak_scores[peak_id == "chr1-100-200", max_variant_samples], 2L)
  expect_true(all(x$peak_scores$MutationPeakScore > 0))
  subgrn <- data.table::data.table(
    TF = c("TF1", "TF1"), Celltype = "AT1", gene = "G1",
    peak_id = c("chr1-100-200", "chr1-300-400"), CellEdgeScore = c(0.5, 0.2)
  )
  y <- calculate_target_mgrs(subgrn, x$peak_scores)
  expect_equal(nrow(y$mutation_edges), nrow(subgrn))
  expect_gt(y$target_mutation_status$Target_MGRS, 0)
})
