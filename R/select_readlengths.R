#' Select read lengths with a consistent offset and clear periodicity
#'
#' A read length is kept if (1) at least `min_agree` of the samples have their start-codon
#' peak at the consensus offset, and (2) in at least `min_frame0` of the samples frame 0 is
#' the dominant frame of the P-sites in the CDS.
#'
#' @param offsets Data frame from find_psite_offsets()
#' @param score Data frame from periodicity_score()
#' @param min_agree Minimum share of samples at the consensus offset (default 0.9)
#' @param min_frame0 Minimum share of samples with frame 0 dominant (default 0.8)
#' @return The selected read lengths (the table of both criteria is printed)
#' @export
select_readlengths <- function(offsets, score, min_agree = 0.9, min_frame0 = 0.8) {

  peaks <- attr(offsets, "peaks")
  if (is.null(peaks)) stop("'offsets' has to come from find_psite_offsets()", call. = FALSE)

  res <- offsets
  res$agree  <- NA_real_
  res$frame0 <- NA_real_

  for (i in seq_len(nrow(res))) {
    rl <- res$readlength[i]
    # share of samples whose start-codon peak is at the consensus offset
    res$agree[i] <- mean(peaks$position[peaks$readlength == rl] == res$peak_position[i])
    # share of samples in which frame 0 is the dominant frame
    s <- score[score$readlength == rl, ]
    res$frame0[i] <- if (nrow(s) > 0) mean(s$dominant_frame == 0) else NA_real_
  }

  res$keep <- !is.na(res$frame0) & res$agree >= min_agree & res$frame0 >= min_frame0

  print(data.frame(readlength = res$readlength, offset = res$peak_position,
                   agree = round(res$agree, 2), frame0 = round(res$frame0, 2), keep = res$keep))
  cat("Selected read lengths:", paste(res$readlength[res$keep], collapse = ", "), "\n")
  res$readlength[res$keep]
}
