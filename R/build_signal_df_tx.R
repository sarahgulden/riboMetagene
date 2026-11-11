#' Signal profiles around start and stop codons, in transcript coordinates
#'
#' For every sample and read length: cuts a window around the start and/or stop codon of
#' each selected transcript out of the mapped 5' ends, filters and normalizes the sites
#' (normalize_mat()), and averages them into one profile. Position 0 is the first
#' nucleotide of the codon; negative positions are towards the 5' end of the transcript.
#'
#' @param mapped Named list of tables from map_to_transcripts() (names = sample_name)
#' @param tx_info The tx_info element of prepare_annotation()
#' @param readlengths Read lengths to use, e.g. 27:30
#' @param features "start", "stop" or both
#' @param upstream Nucleotides before the codon (towards the 5' end of the transcript)
#' @param downstream Nucleotides after the codon (towards the 3' end of the transcript)
#' @param threshold Passed to normalize_mat(): maximum share of zero positions per site, threshold of 0.95 means at least 3 positions not 0
#' @param common_set If TRUE, use only transcripts that are valid for all requested features,
#'   so start and stop profiles come from the same genes
#' @param samples Optional data frame with sample_name plus experiment columns to add
#' @return A data.frame: sample_name, feature, readlength, position, signal (+ experiment columns)
#' @export
build_signal_df_tx <- function(mapped, tx_info, readlengths,
                               features = c("start", "stop"),
                               upstream = 25, downstream = 25,
                               threshold = 0.95, common_set = FALSE, samples = NULL) {

  features <- match.arg(features, several.ok = TRUE)
  if (is.data.frame(mapped)) stop("'mapped' has to be a named list with one table per sample")

  win_len   <- upstream + 3 + downstream           # the codon itself is 3 nt
  positions <- seq(-upstream, length.out = win_len)

  # 1. sites whose whole window lies inside the transcript structure
  #    start: needs `upstream` nt of 5' UTR and `downstream` nt of CDS after the codon
  #    stop:  needs `upstream` nt of CDS before the codon and `downstream` nt of 3' UTR
  #    (cds_len runs to the end of the stop codon)
  sites <- list()
  if ("start" %in% features) {
    keep <- tx_info$utr5_len >= upstream & tx_info$cds_len >= 3 + downstream
    sites$start <- data.frame(transcript_id = tx_info$transcript_id[keep],
                              win_start = tx_info$utr5_len[keep] + 1L - upstream)
  }
  if ("stop" %in% features) {
    stop_first <- tx_info$utr5_len + tx_info$cds_len - 2L          # first nt of the stop codon
    keep <- tx_info$cds_len >= upstream + 3 & tx_info$utr3_len >= downstream
    sites$stop <- data.frame(transcript_id = tx_info$transcript_id[keep],
                             win_start = stop_first[keep] - upstream)
  }

  if (common_set && length(sites) > 1) {
    common <- Reduce(intersect, lapply(sites, function(s) s$transcript_id))
    sites  <- lapply(sites, function(s) s[s$transcript_id %in% common, ])
  }
  for (feat in names(sites)) {
    cat(feat, "codon:", nrow(sites[[feat]]), "transcripts with a complete window\n")
  }

  # 2. one profile per sample, feature and read length
  all_rows <- list()
  for (name in names(mapped)) {
    cat("=== Sample:", name, "===\n")
    tab <- mapped[[name]]
    tab <- tab[tab$readlength %in% readlengths, ]

    for (feat in names(sites)) {
      s <- sites[[feat]]

      # row of every mapped 5' end in the site matrix, and its column in the window
      row <- match(tab$transcript_id, s$transcript_id)       # NA = transcript has no valid window
      col <- tab$tx_pos - s$win_start[row] + 1
      ok  <- !is.na(row) & col >= 1 & col <= win_len
      t_ok <- data.frame(row = row[ok], col = col[ok],
                         readlength = tab$readlength[ok], count = tab$count[ok])

      for (rl in readlengths) {
        t_rl <- t_ok[t_ok$readlength == rl, ]

        # matrix: rows = sites, columns = positions in the window, values = 5' end counts
        mat <- matrix(0, nrow = nrow(s), ncol = win_len)
        mat[cbind(t_rl$row, t_rl$col)] <- t_rl$count

        mat_norm <- normalize_mat(mat, threshold = threshold,
                                  label = paste(name, feat, rl, sep = " | "))
        if (is.null(mat_norm)) next

        all_rows[[paste(name, feat, rl)]] <- data.frame(
          sample_name = name, feature = feat, readlength = rl,
          position = positions, signal = colMeans(mat_norm),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- do.call(rbind, all_rows)
  rownames(df) <- NULL

  # 3. add the experiment columns (tissue, timepoint, ...) if given
  if (!is.null(samples)) {
    meta <- samples[, setdiff(colnames(samples), c("filepath", "filename")), drop = FALSE]
    df <- merge(df, meta, by = "sample_name")
  }
  df
}
