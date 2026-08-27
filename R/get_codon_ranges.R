#' Get start or stop codon ranges, split by strand, extended by a window
#'
#' @param annotation A GRanges annotation object (e.g. from a GFF3/GTF)
#'   with a "type" column that includes "start_codon"/"stop_codon"
#' @param codon Which codon type to extract: "start_codon" or "stop_codon"
#' @param window How many nt to extend on each side of the codon
#' @return A list with $plus and $minus GRanges, each extended by window
#' @export
get_codon_ranges <- function(annotation, codon = "stop_codon", window = 25) {

  # make sure type is character, not factor - avoids comparison issues
  annotation_type <- as.character(annotation$type)

  # check that the requested codon type actually exists in this annotation
  # (catches typos or wrong file early, instead of silently returning 0 rows)
  if (!codon %in% annotation_type) {
    stop("'", codon, "' not found in annotation$type. Available types: ",
         paste(unique(annotation_type), collapse = ", "))
  }

  # extract codon type
  codons <- annotation[annotation_type == codon]

  # keep only complete 3nt codons
  codons_3 <- codons[BiocGenerics::width(codons) == 3]

  # split by strand
  plus  <- codons_3[BiocGenerics::strand(codons_3) == "+"]
  minus <- codons_3[BiocGenerics::strand(codons_3) == "-"]

  # remove pseudochromosomes
  plus  <- plus[GenomeInfoDb::seqnames(plus)  %in% GenomeInfoDb::standardChromosomes(plus)]
  minus <- minus[GenomeInfoDb::seqnames(minus) %in% GenomeInfoDb::standardChromosomes(minus)]

  # extend by window
  plus_range  <- plus  + window
  minus_range <- minus + window

  # make sure windows don't run off the start/end of a chromosome
  # (requires seqlengths to be set on the annotation - if not set, this
  # step is silently skipped, since trim() has nothing to trim against)
  plus_range  <- GenomicRanges::trim(plus_range)
  minus_range <- GenomicRanges::trim(minus_range)

  return(list(plus = plus_range, minus = minus_range))
}

