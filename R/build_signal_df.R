#' Build a long-format signal dataframe across samples, readlengths, and features
#'
#' Runs the full pipeline: for each sample, extracts read lengths from the BAM,
#' builds per-readlength coverage, extracts signal around start and/or stop
#' codons, normalizes, and averages into a single mean signal per position.
#'
#' @param samples Dataframe from list_bam_samples() (needs filepath, sample_name)
#' @param readlengths Which readlengths to process, e.g. 27:30
#' @param chrom_sizes Named vector of chromosome lengths
#' @param start_ranges Output of get_codon_ranges(..., codon = "start_codon")
#' @param stop_ranges Output of get_codon_ranges(..., codon = "stop_codon")
#' @param window The window size used when building start_ranges/stop_ranges
#'   (needed to correctly label positions, e.g. window=25 means positions -25:27)
#' @param threshold Zero-filtering threshold passed to normalize_mat()
#' @return A long dataframe: sample_name, feature, readlength, position, signal
#' @export
build_signal_df <- function(samples, readlengths, chrom_sizes,
                            start_ranges = NULL, stop_ranges = NULL,
                            window = 25, threshold = 0.8) {

  # figure out which features we're actually running (start, stop, or both)
  features <- list()
  if (!is.null(start_ranges)) features[["start"]] <- start_ranges
  if (!is.null(stop_ranges))  features[["stop"]]  <- stop_ranges

  all_rows <- list()

  for (i in 1:nrow(samples)) {

    cat("=== Sample", i, "of", nrow(samples), ":", samples$sample_name[i], "===\n")

    # step 1: get every read's length for this sample
    reads <- get_readlengths(samples$filepath[i])

    # step 2: build per-readlength, per-strand coverage
    cov <- readlength_to_coverage(reads, readlengths, chrom_sizes)

    # step 3: loop over readlengths
    for (rl in readlengths) {

      rl_char <- as.character(rl)

      # step 4: loop over features (start and/or stop)
      for (feat_name in names(features)) {

        feat_ranges <- features[[feat_name]]

        mat <- process_sample_rl(
          cov           = cov[[rl_char]],
          windows_plus  = feat_ranges$plus,
          windows_minus = feat_ranges$minus
        )

        mat_norm <- normalize_mat(
          mat,
          threshold = threshold,
          label = paste(samples$sample_name[i], feat_name, rl_char, sep = " | ")
        )
        if (is.null(mat_norm)) next  # skip if nothing survived filtering

        mean_signal <- colMeans(mat_norm, na.rm = TRUE)

        # build the position labels: window on each side + 3nt codon
        n_pos <- length(mean_signal)
        positions <- seq(-window, -window + n_pos - 1)

        row_id <- paste(samples$sample_name[i], feat_name, rl_char, sep = "_")

        all_rows[[row_id]] <- data.frame(
          sample_name = samples$sample_name[i],
          feature     = feat_name,
          readlength  = rl,
          position    = positions,
          signal      = mean_signal,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  final_df <- do.call(rbind, all_rows)
  rownames(final_df) <- NULL
  return(final_df)
}
