#' Plot the periodicity of the P-sites in the CDS
#'
#' @param frames Data frame from psite_frames() (res$frames of metagene_pipeline())
#' @param type "boxplot" = fraction of P-sites per frame for every read length (one box per
#'   sample); "heatmap" = share of P-sites in frame 0 for every sample and read length;
#'   "bar" = P-sites per frame, summed over the read lengths of each sample and averaged per group
#' @param readlengths Only use these read lengths, e.g. res$readlengths. NULL = all in `frames`
#' @param color_by Column that defines the groups of the "bar" plot (e.g. "timepoint").
#'   It has to be the same for all rows of a sample. Default "sample_name" (one group per sample)
#' @param label_by Column for the row labels of the "heatmap". It has to identify every sample
#'   uniquely (default "sample_name")
#' @param order_by Column to sort the rows of the "heatmap" by (e.g. "timepoint")
#' @param low_colour,high_colour Colours of the "heatmap" for less and more than 1/3 of the
#'   P-sites in frame 0
#' @return A ggplot object
#' @export
plot_periodicity <- function(frames, type = c("boxplot", "heatmap", "bar"),
                             readlengths = NULL, color_by = "sample_name",
                             label_by = "sample_name", order_by = NULL,
                             low_colour = "#B2182B", high_colour = "#2166AC") {
  check_plot_packages("ggplot2")
  type <- match.arg(type)

  if (!is.null(readlengths)) frames <- frames[frames$readlength %in% readlengths, ]
  if (nrow(frames) == 0) stop("No rows left after selecting the read lengths", call. = FALSE)

  # 1. boxplot: one box per read length and frame; dashed line = 1/3 = no periodicity
  if (type == "boxplot") {
    p <- ggplot2::ggplot(frames, ggplot2::aes(x = factor(.data$readlength), y = .data$fraction,
                                              fill = factor(.data$frame))) +
      ggplot2::geom_boxplot(outlier.size = 0.8) +
      ggplot2::geom_hline(yintercept = 1/3, linetype = "dashed") +
      ggplot2::labs(x = "Read length (nt)", y = "Fraction of P-sites in the CDS", fill = "Frame") +
      ggplot2::theme_classic(base_size = 12)
    return(p)
  }

  # 2. heatmap: frame 0 share per sample (rows) and read length (columns)
  if (type == "heatmap") {
    f0 <- frames[frames$frame == 0, ]
    if (anyDuplicated(paste(f0[[label_by]], f0$readlength)) > 0) {
      stop("'", label_by, "' does not identify every sample uniquely", call. = FALSE)
    }

    # row order: by order_by (if given), then by label; the first sample is drawn on top
    ord <- if (is.null(order_by)) order(f0[[label_by]]) else order(f0[[order_by]], f0[[label_by]])
    row_levels <- unique(as.character(f0[[label_by]])[ord])
    f0$row_label <- factor(as.character(f0[[label_by]]), levels = rev(row_levels))

    p <- ggplot2::ggplot(f0, ggplot2::aes(x = factor(.data$readlength), y = .data$row_label,
                                          fill = .data$fraction)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(midpoint = 1/3, low = low_colour, mid = "white",
                                    high = high_colour, name = "Frame 0") +
      ggplot2::labs(x = "Read length (nt)", y = NULL) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    return(p)
  }

  # 3. bar plot: per sample, sum the P-sites over the read lengths and calculate the fraction per frame
  per_sample <- stats::aggregate(n ~ sample_name + frame, data = frames, FUN = sum)
  per_sample$fraction <- per_sample$n / stats::ave(per_sample$n, per_sample$sample_name, FUN = sum)

  # attach the grouping column, then average the samples of each group
  meta <- unique(frames[, unique(c("sample_name", color_by)), drop = FALSE])
  per_sample <- merge(per_sample, meta, by = "sample_name")
  bars <- stats::aggregate(per_sample$fraction,
                           by = list(group = per_sample[[color_by]], frame = per_sample$frame),
                           FUN = mean)
  names(bars)[3] <- "fraction"

  ggplot2::ggplot(bars, ggplot2::aes(x = factor(.data$frame), y = .data$fraction,
                                     fill = .data$group)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_hline(yintercept = 1/3, linetype = "dashed") +
    ggplot2::labs(x = "Frame", y = "Fraction of P-sites in the CDS", fill = color_by) +
    ggplot2::theme_classic(base_size = 12)
}
