#' Transcripts that carry a given tag
#'
#' Helper for select_transcripts(): the IDs of all transcripts with a certain tag
#' in the annotation, e.g. "Ensembl_canonical" in GENCODE.
#'
#' @param annotation GRanges annotation with the columns "tag" and "transcript_id"
#' @param tag The tag to look for
#' @return A character vector of transcript IDs
#' @export
transcripts_with_tag <- function(annotation, tag) {

  tags <- annotation$tag
  if (is.null(tags)) stop("The annotation has no 'tag' column", call. = FALSE)

  # GENCODE stores several tags per row (a list); other annotations have one tag per row
  has_tag <- if (is.character(tags) || is.factor(tags)) {
    as.character(tags) == tag
  } else {
    any(tags == tag)
  }

  unique(as.character(annotation$transcript_id[which(has_tag)]))
}

#' Build a table of transcript structures (5' UTR, CDS and 3' UTR lengths)
#'
#' For every transcript that has exons, a start codon and a stop codon, this gets the
#' length of the transcript and of its 5' UTR, CDS and 3' UTR. Everything is in transcript
#' coordinates (introns removed, position 1 = first nucleotide at the 5' end).
#'
#' @param annotation GRanges (e.g. from rtracklayer::import() of a GENCODE GFF3/GTF) with
#'   the columns type, gene_id and transcript_id, and rows of type exon, start_codon, stop_codon
#' @return A list with
#'   \item{exons_tx}{GRangesList: the exons of every transcript, in 5'->3' order}
#'   \item{tx_info}{data.frame: transcript_id, gene_id, tx_length, utr5_len, cds_len, utr3_len.
#'     cds_len runs from the start codon to the end of the stop codon.}
#' @export
build_transcript_table <- function(annotation) {

  # 1. check that the annotation has what we need
  needed  <- c("type", "gene_id", "transcript_id")
  missing <- setdiff(needed, names(S4Vectors::mcols(annotation)))
  if (length(missing) > 0) {
    stop("The annotation needs the column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  type <- as.character(annotation$type)
  for (t in c("exon", "start_codon", "stop_codon")) {
    if (!t %in% type) stop("No rows of type '", t, "' in the annotation", call. = FALSE)
  }

  # 2. exons ordered from the 5' to the 3' end of each transcript:
  #    ascending coordinates on the + strand, descending on the - strand
  exons     <- annotation[type == "exon"]
  direction <- ifelse(as.character(BiocGenerics::strand(exons)) == "-", -1, 1)
  exons     <- exons[order(as.character(exons$transcript_id),
                           direction * BiocGenerics::start(exons))]
  exons_tx  <- S4Vectors::split(exons, as.character(exons$transcript_id))

  # 3. start and stop codon positions in transcript coordinates
  codon_in_tx <- function(codon) {
    cod    <- annotation[type == codon]
    cod_tx <- GenomicFeatures::mapToTranscripts(cod, exons_tx)   # genomic -> transcript coordinates

    # keep only matches where the codon landed in its own transcript
    same_tx <- names(exons_tx)[cod_tx$transcriptsHits] ==
      as.character(cod$transcript_id[cod_tx$xHits])
    cod_tx  <- cod_tx[same_tx]

    # transcript coordinates always run 5'->3': rebuild as "+" ranges, merge the pieces of a
    # codon that an intron splits (e.g. 1 nt + 2 nt), and keep the complete 3 nt codons
    cod_tx <- GenomicRanges::GRanges(
      seqnames = as.character(GenomeInfoDb::seqnames(cod_tx)),
      ranges   = IRanges::IRanges(BiocGenerics::start(cod_tx), BiocGenerics::end(cod_tx)),
      strand   = "+"
    )
    cod_tx <- IRanges::reduce(cod_tx)
    cod_tx[BiocGenerics::width(cod_tx) == 3]
  }
  start_tx <- codon_in_tx("start_codon")
  stop_tx  <- codon_in_tx("stop_codon")

  # 4. one row per transcript that has both codons
  start_ids <- as.character(GenomeInfoDb::seqnames(start_tx))
  stop_ids  <- as.character(GenomeInfoDb::seqnames(stop_tx))
  tx_ids    <- intersect(start_ids, stop_ids)

  start_first <- BiocGenerics::start(start_tx)[match(tx_ids, start_ids)]   # first nt of the start codon
  stop_last   <- BiocGenerics::end(stop_tx)[match(tx_ids, stop_ids)]       # last nt of the stop codon
  tx_length   <- as.integer(sum(BiocGenerics::width(exons_tx)))[match(tx_ids, names(exons_tx))]
  gene_id     <- as.character(exons$gene_id)[match(tx_ids, as.character(exons$transcript_id))]

  tx_info <- data.frame(
    transcript_id = tx_ids,
    gene_id       = gene_id,
    tx_length     = tx_length,
    utr5_len      = start_first - 1L,               # nt before the start codon
    cds_len       = stop_last - start_first + 1L,   # start codon to end of the stop codon
    utr3_len      = tx_length - stop_last,          # nt after the stop codon
    stringsAsFactors = FALSE
  )

  # 5. sanity check: start codon to end of stop codon must be a multiple of 3
  in_frame <- tx_info$cds_len > 0 & tx_info$cds_len %% 3 == 0
  cat("Transcripts with exons, start codon and stop codon:", nrow(tx_info), "\n")
  cat("Dropped", sum(!in_frame), "whose ORF length is not a multiple of 3\n")
  tx_info <- tx_info[in_frame, ]
  rownames(tx_info) <- NULL

  list(exons_tx = exons_tx[tx_info$transcript_id], tx_info = tx_info)
}

#' Choose one transcript per gene
#'
#' Rules, in this order: (1) transcripts in `prefer` come first (e.g. canonical ones);
#' (2) among the remaining candidates of a gene the fallback decides:
#' "longest_cds" (ties: longest transcript) or "random".
#'
#' @param tx_info The tx_info table from build_transcript_table()
#' @param prefer Optional vector of preferred transcript IDs, e.g. from transcripts_with_tag()
#' @param fallback "longest_cds" (default) or "random"
#' @param seed Seed for fallback = "random". Same annotation + same seed = same result.
#'   Note that this sets the random seed of your R session.
#' @return tx_info with one row per gene
#' @export
select_transcripts <- function(tx_info, prefer = NULL,
                               fallback = c("longest_cds", "random"), seed = 1) {

  fallback <- match.arg(fallback)

  # TRUE for the transcripts the user prefers
  tx_info$preferred <- tx_info$transcript_id %in% prefer

  # sort by transcript_id first, so the result does not depend on the order of the input
  tx_info <- tx_info[order(tx_info$transcript_id), ]

  # order within every gene: preferred transcripts first, then by the fallback rule
  if (fallback == "longest_cds") {
    ord <- order(tx_info$gene_id, !tx_info$preferred, -tx_info$cds_len, -tx_info$tx_length)
  } else {
    set.seed(seed)
    ord <- order(tx_info$gene_id, !tx_info$preferred, stats::runif(nrow(tx_info)))
  }
  tx_info <- tx_info[ord, ]

  # the first transcript of every gene is the chosen one
  chosen <- tx_info[!duplicated(tx_info$gene_id), ]

  cat("Selected", nrow(chosen), "transcripts, one per gene:",
      sum(chosen$preferred), "from the preferred set,",
      sum(!chosen$preferred), "by the fallback rule (", fallback, ")\n")

  chosen$preferred <- NULL
  rownames(chosen) <- NULL
  chosen
}

#' Prepare the annotation: one transcript per gene, with exons and lengths
#'
#' Runs build_transcript_table() and select_transcripts().
#'
#' @inheritParams select_transcripts
#' @param annotation GRanges annotation, see build_transcript_table()
#' @return A list with exons_tx (GRangesList) and tx_info (data.frame) of the chosen transcripts
#' @export
prepare_annotation <- function(annotation, prefer = NULL,
                               fallback = c("longest_cds", "random"), seed = 1) {

  fallback <- match.arg(fallback)
  tx       <- build_transcript_table(annotation)
  tx_info  <- select_transcripts(tx$tx_info, prefer = prefer, fallback = fallback, seed = seed)

  list(exons_tx = tx$exons_tx[tx_info$transcript_id], tx_info = tx_info)
}
