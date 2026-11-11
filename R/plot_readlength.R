#' @importFrom rlang .data
NULL

# internal helper: stop with a clear message if a plotting package is missing
check_plot_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Please install: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

#' Line plot of the read length distribution
#'
#' One line per sample, showing the % of reads at each read length.
#'
#' @param qc Data frame from readlength_distribution()
#' @param color_by Column used to colour the lines, e.g. "timepoint" or "tissue"
#'   (default "sample_name"). For a custom order of the groups, make that column a factor first.
#' @param highlight Read lengths to mark, e.g. the output of suggest_readlengths().
#'   Dashed lines are drawn at the start and end of the range. NULL = no lines.
#' @param length_range Optional c(min, max) to show only part of the read lengths,
#'   e.g. c(24, 40). The percentages stay the same.
#' @return A ggplot object
#' @export
plot_readlength_line <- function(qc, color_by = "sample_name", highlight = NULL,
                                 length_range = NULL) {
  check_plot_packages("ggplot2")

  if (!is.null(length_range)) {
    qc <- qc[qc$readlength >= length_range[1] & qc$readlength <= length_range[2], ]
  }

  # one x-axis tick every 2 nt
  rl <- range(qc$readlength)

  p <- ggplot2::ggplot(qc, ggplot2::aes(x = .data$readlength, y = .data$pct,
                                        colour = .data[[color_by]],
                                        group = .data$sample_name)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_continuous(breaks = seq(rl[1], rl[2], by = 2)) +
    ggplot2::scale_colour_viridis_d(end = 0.9) +
    ggplot2::labs(x = "Read length (nt)", y = "% of reads", colour = color_by) +
    ggplot2::theme_classic(base_size = 12)

  if (!is.null(highlight)) {
    p <- p + ggplot2::geom_vline(xintercept = range(highlight), linetype = "dashed",
                                 colour = "grey40", alpha = 0.6)
  }
  p
}

#' Heatmap of the read length distribution
#'
#' One row per sample, one column per read length, fill = % of reads.
#'
#' @param qc Data frame from readlength_distribution()
#' @param label_by Column used for the row labels. It has to identify every sample
#'   uniquely (default "sample_name").
#' @param order_by Optional column to order the rows by (e.g. "timepoint"), then by label
#' @param highlight Read lengths to mark with dashed lines, e.g. from suggest_readlengths()
#' @param length_range Optional c(min, max) to show only part of the read lengths
#' @param high_colour Fill colour for the highest percentage (low is always very light)
#' @return A ggplot object
#' @export
plot_readlength_heatmap <- function(qc, label_by = "sample_name", order_by = NULL,
                                    highlight = NULL, length_range = NULL,
                                    high_colour = "#17848E") {
  check_plot_packages("ggplot2")

  if (!is.null(length_range)) {
    qc <- qc[qc$readlength >= length_range[1] & qc$readlength <= length_range[2], ]
  }

  # row order: by order_by (if given), then by label; the first sample is drawn on top
  ord <- if (is.null(order_by)) order(qc[[label_by]]) else order(qc[[order_by]], qc[[label_by]])
  row_levels <- unique(as.character(qc[[label_by]])[ord])
  qc$row_label <- factor(as.character(qc[[label_by]]), levels = rev(row_levels))

  rl <- range(qc$readlength)

  p <- ggplot2::ggplot(qc, ggplot2::aes(x = .data$readlength, y = .data$row_label,
                                        fill = .data$pct)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#F0F9F9", high = high_colour, name = "% of reads") +
    ggplot2::scale_x_continuous(breaks = seq(rl[1], rl[2], by = 2)) +
    ggplot2::labs(x = "Read length (nt)", y = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  if (!is.null(highlight)) {
    p <- p + ggplot2::geom_vline(xintercept = range(highlight), linetype = "dashed",
                                 colour = "grey40", alpha = 0.6)
  }
  p
}

#' Read length QC figure: line plot on top, heatmap below
#'
#' @inheritParams plot_readlength_line
#' @param label_by Column used for the heatmap row labels
#' @param order_by Column to order the heatmap rows by. Default: the colour column
#'   (if it is not "sample_name").
#' @param title,subtitle Text above the figure
#' @return A patchwork object (can be saved with ggplot2::ggsave())
#' @export
plot_readlength_qc <- function(qc, color_by = "sample_name", label_by = "sample_name",
                               order_by = NULL, highlight = NULL, length_range = NULL,
                               title = "Read length distribution", subtitle = NULL) {
  check_plot_packages(c("ggplot2", "patchwork"))

  if (is.null(order_by) && color_by != "sample_name") order_by <- color_by

  # line plot: the x-axis is shown only on the heatmap below
  p_line <- plot_readlength_line(qc, color_by, highlight, length_range) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  # colouring by sample: the heatmap rows already name the samples, so no legend
  if (color_by == "sample_name") {
    p_line <- p_line + ggplot2::theme(legend.position = "none")
  }

  p_heat <- plot_readlength_heatmap(qc, label_by, order_by, highlight, length_range)

  # stack vertically, the line plot a bit smaller than the heatmap
  patchwork::wrap_plots(p_line, p_heat, ncol = 1, heights = c(0.6, 1)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}
