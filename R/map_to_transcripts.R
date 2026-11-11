#' Map 5' end count tables to transcript coordinates
#'
#' Converts the genomic 5' end positions from extract_5prime() into positions along the
#' selected transcripts (introns removed, position 1 = first nucleotide of the transcript).
#' Only 5' ends inside an exon of a selected transcript, on the same strand, are kept.
#'
#' How it works: every exon gets a number = how many transcript nucleotides come before it.
#' A 5' end inside an exon then converts to a transcript position by simple arithmetic.
#'
#' @param counts One count table (from bam_to_5prime_counts()) or a named list of them
#'   (from extract_5prime() or read_5prime_tables())
#' @param exons_tx The exons_tx element of prepare_annotation()
#' @return Same structure as `counts`: each table has the columns
#'   transcript_id, tx_pos, readlength, count
#' @export
map_to_transcripts <- function(counts, exons_tx) {

  # 1. one GRanges with all exons (named by transcript ID, in 5'->3' order within a transcript)
  ex <- unlist(exons_tx)

  # 2. for every exon: nucleotides of its transcript that come before it
  #    (cumulative exon width within the transcript, minus the exon's own width)
  tx_offset <- as.integer(ave(BiocGenerics::width(ex), names(ex),
                              FUN = function(w) cumsum(w) - w))

  # a list of samples: map one after the other
  if (!is.data.frame(counts)) {
    result <- list()
    for (name in names(counts)) {
      cat("=== Mapping", name, "===\n")
      result[[name]] <- map_one_sample(counts[[name]], ex, tx_offset)
      invisible(gc())   # free memory before the next sample
    }
    return(result)
  }

  map_one_sample(counts, ex, tx_offset)
}

# internal: map the 5' end table of one sample
map_one_sample <- function(tab, ex, tx_offset) {

  total <- sum(tab$count)

  # keep only 5' ends on chromosomes/scaffolds that the annotation knows
  known <- tab$chr %in% GenomeInfoDb::seqlevels(ex)
  if (!any(known)) {
    stop("None of the chromosome names in the count table are in the annotation. ",
         "Do both use the same naming (e.g. 'chr1' vs '1')?", call. = FALSE)
  }
  cat("Skipping", sum(tab$count[!known]), "reads on sequences that are not in the annotation\n")
  tab <- tab[known, ]

  # one 1 nt range per row, with exactly the same chromosome list as the exons
  gr <- GenomicRanges::GRanges(
    seqnames = factor(tab$chr, levels = GenomeInfoDb::seqlevels(ex)),
    ranges   = IRanges::IRanges(start = tab$pos, width = 1),
    strand   = tab$strand
  )

  # which exon contains each 5' end? (strand-aware: only exons on the same strand match)
  hits <- IRanges::findOverlaps(gr, ex)
  q <- S4Vectors::queryHits(hits)      # row of the count table
  s <- S4Vectors::subjectHits(hits)    # exon

  # position inside the exon, counted from the exon's 5' end:
  # the left end on the + strand, the right end on the - strand
  on_minus    <- as.character(BiocGenerics::strand(ex))[s] == "-"
  pos_in_exon <- ifelse(on_minus,
                        BiocGenerics::end(ex)[s] - tab$pos[q] + 1L,
                        tab$pos[q] - BiocGenerics::start(ex)[s] + 1L)

  out <- data.frame(
    transcript_id = names(ex)[s],
    tx_pos        = tx_offset[s] + pos_in_exon,
    readlength    = tab$readlength[q],
    count         = tab$count[q],
    stringsAsFactors = FALSE
  )

  cat("Mapped", sum(out$count), "of", total, "reads (",
      round(100 * sum(out$count) / total, 1), "%) to the selected transcripts\n")
  out
}
