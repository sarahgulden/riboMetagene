#' Extract 5' end count tables for several BAM files
#'
#' @param samples Data frame from list_bam_samples() (needs filepath and sample_name)
#' @param unique_only Passed on to bam_to_5prime_counts()
#' @param out_dir Optional folder. If given, every table is saved there as <sample_name>.rds,
#'   and samples that already have a file are loaded instead of re-reading the BAM
#' @return A named list with one count table per sample (named by sample_name)
#' @export
extract_5prime <- function(samples, unique_only = FALSE, out_dir = NULL) {

  if (!is.null(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  result <- list()

  for (i in seq_len(nrow(samples))) {

    name <- samples$sample_name[i]
    cat("=== Sample", i, "of", nrow(samples), ":", name, "===\n")

    rds_file <- if (!is.null(out_dir)) file.path(out_dir, paste0(name, ".rds")) else NULL

    # already done in an earlier run? then just load it
    if (!is.null(rds_file) && file.exists(rds_file)) {
      cat("Found a saved table, loading it\n")
      result[[name]] <- readRDS(rds_file)
      next
    }

    # read the BAM; if it is broken (e.g. a truncated download), report it and carry on
    tab <- tryCatch(
      bam_to_5prime_counts(samples$filepath[i], unique_only = unique_only),
      error = function(e) {
        cat("FAILED:", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(tab)) next

    if (!is.null(rds_file)) saveRDS(tab, rds_file)
    result[[name]] <- tab
  }

  cat("Done:", length(result), "of", nrow(samples), "samples\n")
  result
}
