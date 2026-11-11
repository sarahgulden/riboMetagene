#' Find the P-site offset (peak position) per readlength
#'
#' Uses start-codon signal only, since that's the biologically meaningful
#' P-site reference. Finds the most common peak position (mode) across
#' all groups (e.g. samples/tissues/timepoints) for each readlength.
#'
#' @param signal_df Long dataframe from build_signal_df() (needs feature,
#'   readlength, position, signal, and any grouping columns like sample_name)
#' @param group_cols Character vector of columns identifying a "group"
#'   (e.g. c("sample_name") or c("tissue", "timepoint")) - the peak is
#'   found once per group per readlength, then the most common value is kept
#' @return A dataframe: readlength, peak_position
#' @export
find_psite_offsets <- function(signal_df, group_cols = "sample_name") {

  start_df <- signal_df[signal_df$feature == "start" & signal_df$position < 0, ]

  # find the peak (max signal) per group per readlength
  split_cols <- c(group_cols, "readlength")
  peak_list <- split(start_df, start_df[split_cols])

  peaks <- do.call(rbind, lapply(peak_list, function(d) {
    if (nrow(d) == 0) return(NULL)
    d[which.max(d$signal), ]
  }))

  # show how the groups voted: rows = read length, columns = peak position, values = number of groups
  cat("Votes per read length (columns = peak position, values = number of samples):\n")
  print(table(peaks$readlength, peaks$position))

  # count how often each position wins, per readlength
  peak_counts <- as.data.frame(table(peaks$readlength, peaks$position))
  names(peak_counts) <- c("readlength", "position", "n")
  peak_counts <- peak_counts[peak_counts$n > 0, ]

  # keep the most common (mode) position per readlength
  offsets <- do.call(rbind, lapply(split(peak_counts, peak_counts$readlength), function(d) {
    d[which.max(d$n), ]
  }))

  out <- data.frame(
    readlength    = as.numeric(as.character(offsets$readlength)),
    peak_position = as.numeric(as.character(offsets$position)),
    stringsAsFactors = FALSE
  )

  # keep the peak of every group, so select_readlengths() can check how well the groups agree
  attr(out, "peaks") <- peaks
  out

}
