#' Extract signal matrix around codon ranges for one readlength
#'
#' @param cov A list with $plus and $minus Rle coverage
#'   (one element of the output from readlength_to_coverage(), e.g. cov1[["28"]])
#' @param windows_plus GRanges of codon windows on the plus strand
#' @param windows_minus GRanges of codon windows on the minus strand
#' @return A matrix: rows = individual codon sites, columns = position in window
#' @export
process_sample_rl <- function(cov, windows_plus, windows_minus) {

  rle_plus  <- cov$plus
  rle_minus <- cov$minus

  # extract signal matrices
  mat_plus  <- as.matrix(rle_plus[windows_plus])

  mat_minus <- as.matrix(rle_minus[windows_minus])
  mat_minus <- mat_minus[, ncol(mat_minus):1]  # reverse columns (minus strand reads 3'->5')

  # combine strands
  mat <- rbind(mat_plus, mat_minus)

  return(mat)
}
