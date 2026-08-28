#' Run the full ribosome profiling metagene pipeline
#'
#' Wraps the complete workflow: for each BAM file, extracts read lengths,
#' builds per-readlength coverage, extracts signal around start and/or
#' stop codons, normalizes, finds P-site offsets (from start-codon signal),
#' and applies P-site correction. Returns one long dataframe ready for
#' plotting/comparison across samples.
#'
#' @param samples Dataframe from list_bam_samples() (needs filepath, sample_name).
#'   Any extra columns you've added (e.g. tissue, timepoint) are kept and
#'   merged into the final result automatically.
#' @param readlengths Which readlengths to process, e.g. 27:35
#' @param chrom_sizes Named vector of chromosome lengths
#' @param start_ranges Output of get_codon_ranges(..., codon = "start_codon").
#'   Required, since P-site offsets are always calculated from start-codon signal.
#' @param stop_ranges Output of get_codon_ranges(..., codon = "stop_codon").
#'   Optional - omit if you only want start-codon analysis.
#' @param window Window size used when building start_ranges/stop_ranges
#' @param threshold Zero-filtering threshold passed to normalize_mat()
#' @param offset_group_cols Columns identifying a "group" for P-site offset
#'   consistency checking, e.g. "sample_name" or c("tissue", "timepoint")
#' @return A dataframe with sample metadata, feature, readlength, position,
#'   signal, peak_position, and position_corrected
#' @export
metagene_pipeline <- function(samples,
                              readlengths,
                              chrom_sizes,
                              start_ranges,
                              stop_ranges = NULL,
                              window = 25,
                              threshold = 0.8,
                              offset_group_cols = "sample_name") {

  # step 1: build the raw signal dataframe (loops over all samples/readlengths/features)
  cat("=== Step 1/3: extracting signal from all samples ===\n")
  signal_df <- build_signal_df(
    samples      = samples,
    readlengths  = readlengths,
    chrom_sizes  = chrom_sizes,
    start_ranges = start_ranges,
    stop_ranges  = stop_ranges,
    window       = window,
    threshold    = threshold
  )

  # step 2: bring back any extra metadata columns from `samples`
  # (e.g. tissue, timepoint) - merge on whatever columns are shared
  extra_cols <- setdiff(colnames(samples), c("filepath", "filename"))
  meta <- samples[, extra_cols, drop = FALSE]
  signal_df <- merge(signal_df, meta, by = "sample_name")

  # step 3: find P-site offsets from start-codon signal, then apply correction
  cat("=== Step 2/3: finding P-site offsets ===\n")
  offsets <- find_psite_offsets(signal_df, group_cols = offset_group_cols)
  print(offsets)

  cat("=== Step 3/3: applying P-site correction ===\n")
  corrected_df <- apply_psite_correction(signal_df, offsets)

  cat("Done!", nrow(corrected_df), "rows in final dataframe\n")

  return(corrected_df)
}
