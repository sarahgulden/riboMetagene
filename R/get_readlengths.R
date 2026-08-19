#' Read a BAM file and extract read length for each alignment
#'
#' @param bam_file Path to a BAM file
#' @return A GRanges object where each read has a "readlength" column
#' @export
get_readlengths <- function(bam_file) {

  # read the BAM file (GAlignments = one row per aligned read)
  aln <- GenomicAlignments::readGAlignments(bam_file)

  # convert to GRanges (easier to work with downstream)
  gr <- GenomicRanges::granges(aln)

  # add read length as a column
  # qwidth = length of the read as it aligned to the genome
  gr$readlength <- GenomicAlignments::qwidth(aln)

  return(gr)
}
