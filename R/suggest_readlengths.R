#' Suggest read lengths for the analysis
#'
#' @param qc Data frame from readlength_distribution()
#' @param method "cumulative" = the most abundant lengths until they make up `cutoff` of the reads;
#'   "top_n" = the `n` most abundant lengths
#' @param cutoff Fraction of reads to cover for method "cumulative" (default 0.9)
#' @param n Number of lengths for method "top_n" (default 5)
#' @param contiguous If TRUE, fill gaps so the result is one continuous range (e.g. 27:31)
#' @return Integer vector of suggested read lengths
#' @export
suggest_readlengths <- function(qc, method = c("cumulative", "top_n"),
                                cutoff = 0.9, n = 5, contiguous = TRUE) {

  method    <- match.arg(method)
  n_samples <- length(unique(qc$sample_name))

  # average % of reads per read length; every sample counts equally
  avg <- tapply(qc$pct, qc$readlength, sum) / n_samples
  avg <- setNames(as.numeric(avg), names(avg))
  avg <- sort(avg, decreasing = TRUE)

  if (method == "top_n") {
    chosen <- as.integer(names(avg)[seq_len(min(n, length(avg)))])
  } else {
    # take the most abundant lengths until together they reach the cutoff
    n_needed <- which(cumsum(avg) >= cutoff * 100 - 1e-8)[1]
    if (is.na(n_needed)) n_needed <- length(avg)
    chosen <- as.integer(names(avg)[seq_len(n_needed)])
  }

  # fill gaps so the result is one continuous range
  if (contiguous) chosen <- seq(min(chosen), max(chosen))
  chosen <- sort(chosen)

  # how much of each sample's reads do the chosen lengths cover?
  in_range <- qc$readlength %in% chosen
  covered  <- tapply(qc$pct[in_range], qc$sample_name[in_range], sum)

  cat("Suggested read lengths:", min(chosen), "-", max(chosen), "nt\n")
  cat("They cover on average", round(mean(covered), 1), "% of the reads (lowest sample:",
      round(min(covered), 1), "%)\n")

  chosen
}
