#' Frame of every P-site in the CDS
#'
#' Shifts every mapped 5' end to its P-site (5' end - offset) and counts, per sample and
#' read length, in which frame of the CDS the P-sites fall. Frame 0 = first nucleotide of a codon.
#' With correct offsets most P-sites fall in frame 0.
#'
#' @param mapped Named list of tables from map_to_transcripts()
#' @param tx_info The tx_info element of prepare_annotation()
#' @param offsets Data frame from find_psite_offsets() (readlength, peak_position).
#'   Read lengths without an offset are left out.
#' @param skip_codons Codons to leave out at the start of the CDS, and before the stop codon (default 5)
#' @param samples Optional data frame with sample_name plus experiment columns to add
#' @return A data.frame: sample_name, readlength, frame, n (P-sites), fraction (within the read length)
#' @export
psite_frames <- function(mapped, tx_info, offsets, skip_codons = 5, samples = NULL) {

  if (is.data.frame(mapped)) stop("'mapped' has to be a named list with one table per sample")

  all_rows <- list()
  for (name in names(mapped)) {
    cat("=== Sample:", name, "===\n")
    tab <- mapped[[name]]

    # 1. offset of the read length and transcript info of every row; drop rows without them
    off <- offsets$peak_position[match(tab$readlength, offsets$readlength)]
    i   <- match(tab$transcript_id, tx_info$transcript_id)
    ok  <- !is.na(off) & !is.na(i)
    tab <- tab[ok, ]
    off <- off[ok]
    i   <- i[ok]

    # 2. P-site position on the transcript (the offset is negative, e.g. -12,
    #    so subtracting it moves the 5' end 12 nt downstream)
    psite   <- tab$tx_pos - off
    cds_pos <- psite - tx_info$utr5_len[i]        # 1 = first nucleotide of the start codon

    # 3. keep P-sites inside the CDS, without the first and last codons and the stop codon
    #    (cds_len runs to the end of the stop codon)
    keep    <- cds_pos > 3 * skip_codons & cds_pos <= tx_info$cds_len[i] - 3 - 3 * skip_codons
    tab     <- tab[keep, ]
    cds_pos <- cds_pos[keep]
    if (nrow(tab) == 0) { cat("No P-sites left, skipping\n"); next }

    # 4. frame of every P-site: 0 = first nucleotide of a codon
    frame <- factor((cds_pos - 1) %% 3, levels = 0:2)

    # 5. count the P-sites per read length (rows) and frame (columns)
    n <- tapply(tab$count, list(tab$readlength, frame), sum)
    n[is.na(n)] <- 0

    res <- data.frame(
      sample_name = name,
      readlength  = rep(as.integer(rownames(n)), times = 3),
      frame       = rep(0:2, each = nrow(n)),
      n           = as.numeric(n),
      stringsAsFactors = FALSE
    )
    res$fraction <- res$n / ave(res$n, res$readlength, FUN = sum)
    all_rows[[name]] <- res
  }

  df <- do.call(rbind, all_rows)
  rownames(df) <- NULL

  # add the experiment columns (tissue, timepoint, ...) if given
  if (!is.null(samples)) {
    meta <- samples[, setdiff(colnames(samples), c("filepath", "filename")), drop = FALSE]
    df <- merge(df, meta, by = "sample_name")
  }
  df
}

#' Periodicity score per sample and read length
#'
#' From the frame fractions (f0, f1, f2) it calculates the amplitude of the signal at period 3:
#' sqrt(f0^2 + f1^2 + f2^2 - f0*f1 - f0*f2 - f1*f2). It is 0 if the three frames are equal
#' and 1 if all P-sites are in one frame.
#'
#' @param frames Data frame from psite_frames()
#' @return A data.frame: sample_name, readlength, n_psites, frame0 (fraction in frame 0),
#'   dominant_frame, score (plus the experiment columns if frames has them)
#' @export
periodicity_score <- function(frames) {

  groups <- split(frames, paste(frames$sample_name, frames$readlength))

  rows <- lapply(groups, function(d) {
    if (sum(d$n) == 0) return(NULL)
    f <- d$fraction[order(d$frame)]                 # fractions of frame 0, 1, 2
    data.frame(
      sample_name    = d$sample_name[1],
      readlength     = d$readlength[1],
      n_psites       = sum(d$n),
      frame0         = f[1],
      dominant_frame = d$frame[which.max(d$fraction)],
      score          = sqrt(sum(f^2) - f[1] * f[2] - f[1] * f[3] - f[2] * f[3]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # carry over the experiment columns
  meta <- unique(frames[, setdiff(colnames(frames), c("readlength", "frame", "n", "fraction")),
                        drop = FALSE])
  merge(out, meta, by = "sample_name")
}
