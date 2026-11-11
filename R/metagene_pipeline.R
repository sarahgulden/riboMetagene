#' Run the full Ribo-seq metagene pipeline
#'
#' From 5' end tables to P-site-corrected profiles around start and stop codons:
#' maps the 5' ends to transcripts, builds the profiles, finds the P-site offset of every read
#' length, checks the periodicity in the CDS, selects the read lengths with consistent offsets
#' and clear periodicity, and applies the correction.
#'
#' @param counts Named list of 5' end tables (from extract_5prime() or read_5prime_tables()),
#'   or tables that were already mapped with map_to_transcripts(). Already mapped tables skip
#'   the slowest step, which helps when you rerun with other settings.
#' @param annot Output of prepare_annotation() (a list with exons_tx and tx_info)
#' @param samples Optional data frame with sample_name plus your experiment columns
#'   (tissue, timepoint, ...). They are added to the results.
#' @param readlengths Read lengths to consider, e.g. 26:32 (the QC step helps to choose)
#' @param select If TRUE, keep only read lengths that pass select_readlengths(); if FALSE,
#'   keep all of `readlengths`
#' @param min_agree,min_frame0 Passed to select_readlengths()
#' @param features "start", "stop" or both. "start" is always needed, because the offsets
#'   are determined from the start codon.
#' @param upstream,downstream Nucleotides before and after the codon in the profiles
#' @param threshold Passed to build_signal_df_tx(): maximum share of zero positions per site
#' @param common_set Passed to build_signal_df_tx()
#' @param skip_codons Passed to psite_frames(): codons skipped at both ends of the CDS
#' @param trim_edges If TRUE, keep only corrected positions that all selected read lengths cover
#' @param keep_mapped If TRUE, the mapped tables are returned too (they are large)
#' @return A list with
#'   \item{corrected}{profiles of the selected read lengths with position_corrected}
#'   \item{offsets}{P-site offset per read length (with the votes in attr "peaks")}
#'   \item{periodicity}{periodicity score per sample and read length}
#'   \item{frames}{P-sites per frame, sample and read length}
#'   \item{signal}{profiles of all considered read lengths, before correction}
#'   \item{readlengths}{the read lengths that were used}
#'   \item{mapped}{only if keep_mapped = TRUE}
#' @export
metagene_pipeline <- function(counts, annot, samples = NULL,
                              readlengths = 26:32,
                              select = TRUE, min_agree = 0.9, min_frame0 = 0.8,
                              features = c("start", "stop"),
                              upstream = 25, downstream = 25,
                              threshold = 0.95, common_set = FALSE,
                              skip_codons = 5, trim_edges = TRUE,
                              keep_mapped = FALSE) {

  # 0. check the input before the slow steps start
  if (is.data.frame(counts)) {
    stop("'counts' has to be a named list with one table per sample", call. = FALSE)
  }
  if (!all(c("exons_tx", "tx_info") %in% names(annot))) {
    stop("'annot' has to come from prepare_annotation()", call. = FALSE)
  }
  features <- match.arg(features, c("start", "stop"), several.ok = TRUE)
  if (!"start" %in% features) {
    stop("The offsets come from the start codon, so 'features' has to include \"start\"",
         call. = FALSE)
  }
  if (!is.null(samples)) {
    missing <- setdiff(names(counts), samples$sample_name)
    if (length(missing) > 0) {
      stop("These samples are missing in 'samples': ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
  }

  # 1. 5' ends -> transcript coordinates (skipped if the tables are already mapped)
  cat("##### Step 1/5: mapping 5' ends to transcripts #####\n")
  already_mapped <- all(c("transcript_id", "tx_pos") %in% colnames(counts[[1]]))
  if (already_mapped) {
    cat("The tables are already mapped, skipping\n")
    mapped <- counts
  } else {
    mapped <- map_to_transcripts(counts, annot$exons_tx)
  }

  # 2. profiles around the codons, for all considered read lengths
  cat("##### Step 2/5: signal profiles #####\n")
  signal_df <- build_signal_df_tx(mapped, annot$tx_info, readlengths = readlengths,
                                  features = features, upstream = upstream,
                                  downstream = downstream, threshold = threshold,
                                  common_set = common_set, samples = samples)

  # 3. P-site offset per read length (the samples vote)
  cat("##### Step 3/5: P-site offsets #####\n")
  offsets <- find_psite_offsets(signal_df, group_cols = "sample_name")

  # 4. periodicity: in which frame of the CDS do the P-sites fall?
  cat("##### Step 4/5: periodicity in the CDS #####\n")
  frames <- psite_frames(mapped, annot$tx_info, offsets, skip_codons = skip_codons,
                         samples = samples)
  score  <- periodicity_score(frames)

  # 5. choose the read lengths and apply the correction
  cat("##### Step 5/5: read lengths and correction #####\n")
  if (select) {
    use_rl <- select_readlengths(offsets, score, min_agree = min_agree, min_frame0 = min_frame0)
    if (length(use_rl) == 0) {
      stop("No read length passed the checks. Try lower values for min_agree / min_frame0, ",
           "or select = FALSE.", call. = FALSE)
    }
  } else {
    use_rl <- intersect(readlengths, offsets$readlength)
    cat("Using all requested read lengths:", paste(use_rl, collapse = ", "), "\n")
  }

  offsets_use  <- offsets[offsets$readlength %in% use_rl, ]
  corrected_df <- apply_psite_correction(signal_df[signal_df$readlength %in% use_rl, ],
                                         offsets_use)

  # every read length is shifted by a different amount, so at the edges some read
  # lengths have no data: keep only the positions that all selected read lengths cover
  if (trim_edges) {
    lo <- max(tapply(corrected_df$position_corrected, corrected_df$readlength, min))
    hi <- min(tapply(corrected_df$position_corrected, corrected_df$readlength, max))
    corrected_df <- corrected_df[corrected_df$position_corrected >= lo &
                                   corrected_df$position_corrected <= hi, ]
    cat("Kept corrected positions", lo, "to", hi, "(covered by all selected read lengths)\n")
  }
  rownames(corrected_df) <- NULL

  result <- list(corrected = corrected_df, offsets = offsets, periodicity = score,
                 frames = frames, signal = signal_df, readlengths = use_rl)
  if (keep_mapped) result$mapped <- mapped

  cat("Done.\n")
  result
}
