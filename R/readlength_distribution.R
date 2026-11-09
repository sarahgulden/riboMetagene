#' Read length distribution per sample
#'
#' @param counts Named list of count tables from extract_5prime()
#' @param samples Optional data frame with sample_name plus experiment columns
#'   (tissue, timepoint, ...) that are added to the result
#' @return A data.frame: sample_name, readlength, n (reads), pct (% of that sample's reads),
#'   plus the experiment columns if samples was given
#' @export
readlength_distribution <- function(counts, samples = NULL) {

  per_sample <- lapply(names(counts), function(name) {
    tab <- counts[[name]]

    # sum the counts per read length (same as the tapply check you just ran)
    n <- tapply(tab$count, tab$readlength, sum)

    data.frame(
      sample_name = name,
      readlength  = as.integer(names(n)),
      n           = as.numeric(n),
      pct         = 100 * as.numeric(n) / sum(n)
    )
  })
  df <- do.call(rbind, per_sample)

  # add the experiment columns (tissue, timepoint, ...) if given
  if (!is.null(samples)) {
    meta <- samples[, setdiff(colnames(samples), c("filepath", "filename")), drop = FALSE]
    df <- merge(df, meta, by = "sample_name")
  }
  df
}
