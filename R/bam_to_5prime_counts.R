#' Extract the 5' ends of all reads in a BAM file as a compact count table
#'
#' Reads the BAM once and keeps, for every read, only its 5' end (1 nt, strand-aware)
#' and its read length. Identical entries are merged and counted, so the result
#' is much smaller than the BAM and does not depend on any annotation.
#'
#' @param bam_file Path to a BAM file
#' @param unique_only If TRUE, keep only uniquely mapped alignments (NH tag == 1)
#' @return A data.frame with columns: chr, pos, strand, readlength, count
#' @export
bam_to_5prime_counts <- function(bam_file, unique_only = FALSE, check_complete = TRUE) {

  # stop early on truncated files; the wrapper reports the failure and moves on
  if (check_complete && !bam_is_complete(bam_file)) {
    stop("BAM file looks truncated (end-of-file marker missing): ", basename(bam_file))
  }

  # 1. read all alignments, plus the NH tag (= number of places a read maps to)
  aln <- GenomicAlignments::readGAlignments(
    bam_file,
    param = Rsamtools::ScanBamParam(tag = "NH")
  )
  cat("Read", length(aln), "alignments from", basename(bam_file), "\n")

  # 2. multi-mappers: report them, and optionally remove them
  nh <- S4Vectors::mcols(aln)$NH
  if (!is.null(nh)) {
    n_multi <- sum(nh > 1, na.rm = TRUE)
    cat(n_multi, "alignments have NH > 1 (multi-mappers)\n")
    if (unique_only && n_multi > 0) aln <- aln[which(nh == 1)]
  }

  # 3. read length of every alignment, and its 5' end (strand-aware: leftmost on +, rightmost on -)
  readlength <- GenomicAlignments::qwidth(aln)
  five_prime <- GenomicRanges::resize(GenomicRanges::granges(aln), width = 1, fix = "start")

  # 4. plain table: one row per alignment
  df <- data.frame(
    chr        = as.character(GenomeInfoDb::seqnames(five_prime)),
    pos        = BiocGenerics::start(five_prime),
    strand     = as.character(BiocGenerics::strand(five_prime)),
    readlength = readlength
  )

  # 5. merge identical rows and count them: one row per position / strand / read length
  dplyr::count(df, dplyr::across(dplyr::everything()), name = "count")
}
