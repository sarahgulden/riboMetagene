#' Get start or stop codon ranges, split by strand, extended by a window
#'
#' @param annotation A GRanges annotation object (e.g. from a GFF3/GTF)
#'   with a "type" column that includes "start_codon"/"stop_codon"
#' @param codon Which codon type to extract: "start_codon" or "stop_codon"
#' @param window How many nt to extend on each side of the codon
#' @return A list with $plus and $minus GRanges, each extended by window
#' @export
get_codon_ranges <- function(annotation, codon = "stop_codon", window = 25) {

  annotation_type <- as.character(annotation$type)

  if (!codon %in% annotation_type) {
    stop("'", codon, "' not found in annotation$type. Available types: ",
         paste(unique(annotation_type), collapse = ", "))
  }

  codons <- annotation[annotation_type == codon]
  codons_3 <- codons[BiocGenerics::width(codons) == 3]

  plus  <- codons_3[BiocGenerics::strand(codons_3) == "+"]
  minus <- codons_3[BiocGenerics::strand(codons_3) == "-"]

  # explicit as.character() avoids an Rle/%in% dispatch issue that only
  # shows up when this function runs inside a package (not interactively)
  plus  <- plus[as.character(GenomeInfoDb::seqnames(plus))   %in% GenomeInfoDb::standardChromosomes(plus)]
  minus <- minus[as.character(GenomeInfoDb::seqnames(minus)) %in% GenomeInfoDb::standardChromosomes(minus)]

  plus_range  <- plus  + window
  minus_range <- minus + window

  plus_range  <- GenomicRanges::trim(plus_range)
  minus_range <- GenomicRanges::trim(minus_range)

  return(list(plus = plus_range, minus = minus_range))
}
