.ctmrs_as_dt <- function(x) {
  data.table::as.data.table(data.table::copy(x))
}

.ctmrs_require_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(object_name, " lacks columns: ", paste(missing, collapse = ", "),
         ". Available columns: ", paste(names(x), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.ctmrs_write <- function(x, output_dir, filename) {
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(x, file.path(output_dir, filename), sep = "\t")
  }
  invisible(x)
}

.ctmrs_read <- function(x, label) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop(label, " file not found: ", x, call. = FALSE)
    return(readRDS(x))
  }
  x
}

.ctmrs_add_grn_aliases <- function(grn) {
  grn <- .ctmrs_as_dt(grn)
  aliases <- list(
    Regions = "peak_id",
    cor = "PeakGeneCor",
    pvalue = "PeakGeneP",
    Strength = "PeakGeneStrength"
  )
  for (new_name in names(aliases)) {
    old_name <- aliases[[new_name]]
    if (!new_name %in% names(grn) && old_name %in% names(grn)) {
      grn[, (new_name) := get(old_name)]
    }
  }
  grn
}

.ctmrs_validate_probability <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop(name, " must be a single number between 0 and 1.", call. = FALSE)
  }
}

