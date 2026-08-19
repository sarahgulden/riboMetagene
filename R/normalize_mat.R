#' Filter and normalize a signal matrix
#'
#' Removes rows (codon sites) with too many zero-signal positions,
#' then applies min-max normalization to each remaining row.
#'
#' @param mat A signal matrix (rows = sites, columns = position in window)
#' @param threshold Maximum allowed proportion of zero values per row (default 0.8)
#' @param label Optional text describing this matrix (e.g. sample/readlength/feature),
#'   printed alongside the filtering summary
#' @return A normalized matrix, or NULL if no rows survive filtering
#' @export
normalize_mat <- function(mat, threshold = 0.8, label = NULL) {

  # filter rows with too many zeros
  zero_prop    <- rowSums(mat == 0) / ncol(mat)
  mat_filtered <- mat[zero_prop < threshold, , drop = FALSE]

  prefix <- if (!is.null(label)) paste0("[", label, "] ") else ""
  cat(prefix, "Kept", nrow(mat_filtered), "of", nrow(mat), "rows after filtering\n")

  if (is.null(nrow(mat_filtered)) || nrow(mat_filtered) == 0) {
    cat(prefix, "No rows survived filtering — skipping\n")
    return(NULL)
  }

  # min-max normalization per row
  minmax_norm <- function(x) {
    if (max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }

  mat_norm <- t(apply(mat_filtered, 1, minmax_norm))

  return(mat_norm)
}
