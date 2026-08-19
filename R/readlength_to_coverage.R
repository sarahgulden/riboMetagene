#' Convert tagged reads into per-readlength, per-strand coverage (Rle)
#'
#' Uses only the 5' end of each read (not the full footprint width),
#' since P-site offset analysis is based on where reads start, not
#' the entire region they cover.
#'
#' @param gr GRanges with a "readlength" column (from get_readlengths())
#' @param readlengths Which readlengths to keep, e.g. 27:30
#' @param chrom_sizes Named vector of chromosome lengths (names = chromosome names)
#' @return A named list: one element per readlength, each containing
#'   a list with $plus and $minus Rle coverage
#' @export
readlength_to_coverage <- function(gr, readlengths, chrom_sizes) {

  gr <- gr[gr$readlength %in% readlengths]

  GenomeInfoDb::seqlengths(gr) <- chrom_sizes[GenomeInfoDb::seqlevels(gr)]

  # collapse each read down to just its 5' end (1 nt), strand-aware
  gr <- GenomicRanges::resize(gr, width = 1, fix = "start")

  result <- list()

  for (rl in readlengths) {

    gr_rl <- gr[gr$readlength == rl]

    gr_plus  <- gr_rl[GenomicRanges::strand(gr_rl) == "+"]
    gr_minus <- gr_rl[GenomicRanges::strand(gr_rl) == "-"]

    rle_plus  <- GenomicRanges::coverage(gr_plus)
    rle_minus <- GenomicRanges::coverage(gr_minus)

    result[[as.character(rl)]] <- list(plus = rle_plus, minus = rle_minus)
  }

  return(result)
}
