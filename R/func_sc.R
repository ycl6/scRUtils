#' Plot G2M against G1 scores from `cyclone` results
#'
#' This function produces a scatter plot to show the predicted G1 and
#' G2M scores and cell cycle phases of single cells after performing
#' classification using `cyclone()` from the \pkg{scran} package.
#'
#' @param x A `list` returned by `cyclone()` containing `phases`, `scores`
#' and `normalized.scores` slots.
#' @param phase_color A character vector of color codes indicating the colour
#' of the cell cycle phases. Default is `NULL`, and use [choosePalette()] to
#' select a palette.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.8.
#' @param text_size A numeric scalar indicating the size of cell cycle phase
#' labels. Default is 6.
#' @param xlab The title of the x-axis. Default is "G1 score".
#' @param ylab The title of the y-axis. Default is "G2/M score".
#' @param title Plot title. Default is `NULL`.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#'
#' @return A `ggplot` object
#'
#' @details NULL
#'
#' @author I-Hsuan Lin
#'
#' @name plotCyclone
#'
#' @seealso [scran::cyclone()]
#'
#' @export
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(phases_assignments)
#'
#' # Create plot, minimal input
#' plotCyclone(phases_assignments)
#'
#' # Assign colours
#' my.cols <- c("red", "blue", "black")
#' plotCyclone(phases_assignments, phase_color = my.cols)
#'
#' # Assign named colours
#' my.cols <- setNames(c("red", "blue", "black"), c("G1", "S", "G2M"))
#' plotCyclone(phases_assignments, phase_color = my.cols)
plotCyclone <- function(x, phase_color = NULL, point_size = 2, point_alpha = 0.8, text_size = 6,
                        xlab = "G1 score", ylab = "G2/M score", title = NULL, theme_size = 18) {
  .check_phase(x)
  dat <- data.frame(x = x$score$G1, y = x$score$G2M, phases = x$phases, stringsAsFactors = TRUE)

  # Set up plot
  aes <- aes_string(x = "x", y = "y", color = "phases", shape = "phases")
  phase_color <- choosePalette(dat$phases, phase_color)
  phase_shape <- setNames(c(0, 2, 8), c("G1", "S", "G2M"))[names(phase_color)]

  p <- ggplot(dat, aes) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = phase_color) +
    scale_shape_manual(values = phase_shape) +
    guides(
      color = guide_legend("Phase", override.aes = list(size = point_size * 2, alpha = 1)),
      shape = guide_legend("Phase")
    ) +
    scale_x_continuous(xlab, limits = c(0, 1)) +
    scale_y_continuous(ylab, limits = c(0, 1)) +
    geom_segment(aes(x = 1 / 2, y = 0, xend = 1 / 2, yend = 1 / 2), size = 0.5, colour = "black", linetype = "dashed") +
    geom_segment(aes(x = 0, y = 1 / 2, xend = 1 / 2, yend = 1 / 2), size = 0.5, colour = "black", linetype = "dashed") +
    geom_segment(aes(x = 1 / 2, y = 1 / 2, xend = 1, yend = 1), size = 0.5, colour = "black", linetype = "dashed") +
    annotate("label", x = 0.25, y = 0.02, size = text_size, fill = "white", alpha = 0.5, label = "S") +
    annotate("label", x = 0.95, y = 0.50, size = text_size, fill = "white", alpha = 0.5, label = "G1") +
    annotate("label", x = 0.50, y = 0.98, size = text_size, fill = "white", alpha = 0.5, label = "G2/M") +
    theme_cowplot(theme_size)

  # Add plot title
  if (!is.null(title)) p + ggtitle(title) else p
}

#' Plot expression frequency against mean counts
#'
#' This function uses `plotRowData()` from the \pkg{scater} package to
#' produce a scatter plot of expression frequency against mean counts by
#' using per-feature quality control metrics produced by \pkg{scuttle}.
#'
#' @param sce A `SingleCellExperiment` object containing the required QC
#' metrics for features ('mean' and 'detected') that were returned by
#' `addPerCellQC()` or `perFeatureQCMetrics()` in its `rowData` slot.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.8.
#' @param anno_size A numeric scalar indicating the size of annotation text.
#' Default is 6.
#' @param text_size A numeric scalar indicating the size of the label. This
#' is passed to `geom_text_repel()`. Default is 4.
#' @param text_color A string indicating the colour of the label. This is
#' passed to `geom_text_repel()`. Default is "black".
#' @param trend_color A string indicating the colour of the smoothed curve.
#' Default is "firebrick",
#' @param trend_size A numeric scalar indicating the size of the smoothed
#' curve. Default is 1.
#' @param trend_se Logical scalar indicating whether to display confidence
#' interval around smooth. Default is `TRUE`.
#' @param box.padding A scalar indicating the amount of padding around
#' bounding box, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param max.overlaps Exclude text labels that overlap too many things.
#' This is passed to `geom_text_repel()`. Default is Inf.
#' @param seed Random seed passed to [set.seed()]. This is passed to
#' `geom_text_repel()`. Default is 12321.
#' @param xlab The title of the x-axis. Default is "log2 Mean expression".
#' @param ylab The title of the y-axis. Default is "Percentage of expressing
#' cells".
#' @param title Plot title. Default is `NULL`.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#' @param ... Other arguments passed on to `geom_text_repel()`.
#'
#' @return A `ggplot` object
#'
#' @details The original function of the same name can be found in the legacy
#' \pkg{scater} package and is now deprecated.
#'
#' @author I-Hsuan Lin
#'
#' @name plotExprsFreqVsMean
#'
#' @seealso [scater::plotRowData()], [scuttle::addPerFeatureQCMetrics()],
#' [scuttle::perFeatureQCMetrics()], [ggrepel::geom_text_repel()]
#'
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData
#' @importFrom scater plotRowData
#' @importFrom rlang abort
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' plotExprsFreqVsMean(sce, title = "Expression frequency vs. mean expression level")
plotExprsFreqVsMean <- function(sce, point_size = 2, point_alpha = 0.8, anno_size = 6,
                                text_size = 4, text_color = "black", trend_color = "firebrick",
                                trend_size = 1, trend_se = TRUE, box.padding = 0.5,
                                max.overlaps = Inf, seed = 12321, xlab = "log2 Mean expression",
                                ylab = "Percentage of expressing cells", title = NULL,
                                theme_size = 18, ...) {
  .is.sce(sce)
  if (!all(c("mean", "detected") %in% colnames(rowData(sce)))) {
    abort(sprintf("Unable to find the required columns 'mean' and/or 'detected' in `rowData(%s)`.", deparse(substitute(sce))))
  }

  freqs <- rowData(sce)$detected
  means <- log2(rowData(sce)$mean + 1)
  pct_dropout <- 100 - freqs

  ## data frame with expression mean and frequency.
  mn_vs_fq <- data.frame(mn = means, fq = freqs / 100)
  text_x_loc <- min(mn_vs_fq$mn) + 0.6 * diff(range(mn_vs_fq$mn))

  p <- plotRowData(sce,
    x = data.frame("X" = means, check.names = FALSE),
    y = data.frame("Y" = freqs, check.names = FALSE),
    colour_by = I(pct_dropout), point_size = point_size,
    point_alpha = point_alpha, theme_size = theme_size
  ) +
    geom_smooth(
      data = mn_vs_fq, aes_string(x = "mn", y = "100*fq"),
      alpha = 1, color = trend_color, size = trend_size, se = trend_se
    ) +
    geom_hline(yintercept = 50, linetype = 2) + # 50% dropout
    scale_y_continuous(limits = c(0, 110), breaks = c(seq(0, 100, by = 25))) +
    annotate("text",
      x = text_x_loc, y = 40, size = anno_size,
      label = paste0(sum(mn_vs_fq$fq >= 0.50), " genes are expressed\nin at least 50% of cells")
    ) +
    annotate("text",
      x = text_x_loc, y = 20, size = anno_size,
      label = paste0(sum(mn_vs_fq$fq >= 0.25), " genes are expressed\nin at least 25% of cells")
    ) +
    labs(x = xlab, y = ylab, color = "Dropouts (%)")

  # Add labels
  aes <- aes_string(label = "ifelse(X > quantile(means, 0.999) & Y > quantile(freqs, 0.999), names(sce), '')")
  p <- p + geom_text_repel(aes,
    alpha = 1,
    color = text_color, box.padding = box.padding, size = text_size,
    max.overlaps = max.overlaps, seed = seed, min.segment.length = unit(0, "lines"), ...
  )

  # Add plot title
  if (!is.null(title)) p + ggtitle(title) else p
}

#' Plot logcounts variance against logcounts mean
#'
#' This function calculates the per-feature mean and variance using normalised
#' counts (accessible via `logcounts()`) from a `SingleCellExperiment` object,
#' and returns a logcounts mean-variance scatter plot.
#'
#' @param sce A `SingleCellExperiment` object containing the required QC
#' metrics returned by `addPerCellQC()`/`perFeatureQCMetrics()`.
#' @param top_n An integer scalar indicating the first `top_n` genes with
#' the highest calculated mean to label. Default is 5.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.8.
#' @param text_size A numeric scalar indicating the size of the label. This
#' is passed to `geom_text_repel()`. Default is 4.
#' @param text_color A string indicating the colour of the label. This is
#' passed to `geom_text_repel()`. Default is "black".
#' @param box.padding A scalar indicating the amount of padding around
#' bounding box, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param max.overlaps Exclude text labels that overlap too many things.
#' This is passed to `geom_text_repel()`. Default is Inf.
#' @param seed Random seed passed to [set.seed()]. This is passed to
#' `geom_text_repel()`. Default is 12321.
#' @param xlab The title of the x-axis. Default is "Mean log-counts".
#' @param ylab The title of the y-axis. Default is "Variance of log-counts".
#' @param title Plot title. Default is `NULL`.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#' @param ... Other arguments passed on to `geom_text_repel()`.
#'
#' @return A `ggplot` object
#'
#' @details When the `rowData` slot of the `SingleCellExperiment` input
#' contains QC metrics for features that were returned by `addPerCellQC()` or
#' `perFeatureQCMetrics()`, this function will use `detected` (percentage of
#' expressed features above the detection limit) to calculate `pct_dropout`
#' (percentage of dropouts) and colour the points accordingly.
#'
#' @author I-Hsuan Lin
#'
#' @name plotVarianceVsMean
#'
#' @seealso [scater::plotRowData()], [scuttle::logNormCounts()],
#' [scuttle::calculateCPM()], [ggrepel::geom_text_repel()]
#'
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment logcounts
#' @importFrom scater plotRowData
#' @importFrom utils head
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' plotVarianceVsMean(sce, title = "logcounts mean-variance plot")
plotVarianceVsMean <- function(sce, top_n = 5, point_size = 2, point_alpha = 0.8,
                               text_size = 4, text_color = "black",
                               box.padding = 0.5, max.overlaps = Inf, seed = 12321,
                               xlab = "Mean log-counts", ylab = "Variance of log-counts",
                               title = NULL, theme_size = 18, ...) {
  .is.sce(sce)
  .check_assayname(sce, "logcounts")

  mean <- .clac_rowMeans(logcounts(sce))
  variance <- .clac_rowVars(logcounts(sce))

  rowData(sce)$rowname <- if (is.null(rownames(sce))) NULL else rownames(sce)

  if ("detected" %in% colnames(rowData(sce))) {
    pct_dropout <- 100 - rowData(sce)$detected
    p <- plotRowData(sce,
      x = data.frame("X" = mean, check.names = FALSE),
      y = data.frame("Y" = variance, check.names = FALSE),
      point_size = point_size, point_alpha = point_alpha, theme_size = theme_size,
      other_fields = "rowname", colour_by = I(pct_dropout)
    ) +
      guides(color = guide_colourbar(title = "Dropouts (%)"))
  } else {
    p <- plotRowData(sce,
      x = data.frame("X" = mean, check.names = FALSE),
      y = data.frame("Y" = variance, check.names = FALSE),
      point_size, point_alpha = point_alpha, theme_size = theme_size,
      other_fields = "rowname"
    )
  }

  p <- p + labs(x = xlab, y = ylab)

  # Add labels
  sel <- names(head(sort(mean, decreasing = TRUE), top_n))
  aes_repel <- aes_string(label = 'ifelse(rowname %in% sel, rowname, "")')
  p <- p + geom_text_repel(aes_repel,
    alpha = 1, color = text_color, box.padding = box.padding, size = text_size,
    max.overlaps = max.overlaps, seed = seed, min.segment.length = unit(0, "lines"), ...
  )

  # Add plot title
  if (!is.null(title)) p + ggtitle(title) else p
}

#' Plot variance modelling results
#'
#' This function requires a `SingleCellExperiment` object and the resulting
#' DataFrame after performing variance modelling, such as `modelGeneVar()`
#' from the \pkg{scran} package, and produces a mean-variance scatter plot.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param var A DataFrame returned by variance modelling functions such as
#' `modelGeneVar()`, `modelGeneVarByPoisson()`, etc.
#' @param hvg A character vector containing a set of highly variable genes
#' to highlight in red. Default is `NULL`.
#' @param top_n An integer scalar indicating the first `top_n` genes to label.
#' When `hvg = NULL`, the first `top_n` genes with the highest `mean` are
#' selected. When HVGs are provided, the first `top_n` HVGs are selected.
#' Default is 10.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param text_size A numeric scalar indicating the size of the label. This
#' is passed to `geom_text_repel()`. Default is 4.
#' @param text_color A string indicating the colour of the label. This is
#' passed to `geom_text_repel()`. Default is "black".
#' @param trend_color A string indicating the colour of the smoothed curve.
#' Default is "gold",
#' @param trend_size A numeric scalar indicating the size of the smoothed
#' curve. Default is 2.
#' @param box.padding A scalar indicating the amount of padding around
#' bounding box, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param max.overlaps Exclude text labels that overlap too many things.
#' This is passed to `geom_text_repel()`. Default is Inf.
#' @param seed Random seed passed to [set.seed()]. This is passed to
#' `geom_text_repel()`. Default is 12321.
#' @param title Plot title. Default is `NULL`.
#' @param xlab The title of the x-axis. Default is "Mean log-counts".
#' @param ylab The title of the y-axis. Default is "Variance of log-counts".
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#' @param ... Other arguments passed on to `geom_text_repel()`.
#'
#' @return A `ggplot` object
#'
#' @details
#' In addition to the two required inputs, if `rowData(sce)` has `is_pass`
#' and/or `is_ambient` columns containing logical values (i.e. `TRUE` and
#' `FALSE`), denoting if a gene passed QC or if its expression is attributable
#' to ambient contamination, the resulting figure will show the genes in
#' different shapes. Also, genes that are HVGs are coloured in red when a gene
#' list is provided using the `hvg` argument.
#'
#' @author I-Hsuan Lin
#'
#' @name plotVariableFeature
#'
#' @seealso [scran::modelGeneVar()], [scran::modelGeneVarByPoisson()],
#' [scran::modelGeneVarWithSpikes()], [ggrepel::geom_text_repel()]
#'
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr pull
#' @importFrom utils head
#' @importFrom stats setNames
#' @importFrom stats quantile
#' @importFrom tibble add_column
#' @importFrom tibble rownames_to_column
#' @importFrom S4Vectors metadata
#' @importFrom rlang abort
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' var <- metadata(sce)[["modelGeneVar"]]
#' hvg <- metadata(sce)[["HVG"]]
#'
#' # View modelGeneVar result
#' plotVariableFeature(sce, var, top_n = 5, title = "modelGeneVar")
#'
#' # Highlight HVGs
#' plotVariableFeature(sce, var, hvg, title = "modelGeneVar (highlight HVGs)")
plotVariableFeature <- function(sce, var, hvg = NULL, top_n = 10, point_size = 2,
                                text_size = 4, text_color = "black",
                                trend_color = "gold", trend_size = 2,
                                box.padding = 0.5, max.overlaps = Inf, seed = 12321,
                                xlab = "Mean log-counts", ylab = "Variance of log-counts",
                                title = NULL, theme_size = 18, ...) {
  .is.sce(sce)
  if (!all(c("mean", "total", "p.value", "FDR") %in% colnames(var))) abort("`var` must be an object returned by one of the `modelGene` functions by `scran`.")
  if (!all(rownames(var) %in% rownames(sce))) abort(sprintf("Not all the genes in `var` are found in `%s`.", deparse(substitute(sce))))

  dat <- as.data.frame(var) %>% rownames_to_column()
  if (!is.null(hvg)) {
    dat <- dat %>%
      add_column(HVG = FALSE) %>%
      mutate(HVG = replace(.data$HVG, .data$rowname %in% hvg, TRUE))
  }

  dat$status <- "Passed"
  is_fail <- is_ambient <- c()
  if ("is_pass" %in% names(rowData(sce))) {
    is_fail <- rownames(sce)[rowData(sce)$is_pass == FALSE]
    dat$status[which(dat$rowname %in% is_fail)] <- "Failed"
  }

  if ("is_ambient" %in% names(rowData(sce))) {
    is_ambient <- rownames(sce)[rowData(sce)$is_ambient == TRUE]
    dat$status[which(dat$rowname %in% is_ambient)] <- "Ambient Contamination"
  }

  dat$status <- as.factor(dat$status)

  # Set up plot
  aes <- aes_string(x = "mean", y = "total", shape = "status", alpha = "status")

  p <- ggplot(dat, aes) +
    theme_cowplot(theme_size) +
    scale_shape_manual(values = setNames(c(4, 1, 16), c("Failed", "Ambient Contamination", "Passed"))) +
    scale_alpha_manual(values = setNames(c(1, 1, 0.4), c("Failed", "Ambient Contamination", "Passed"))) +
    guides(
      shape = guide_legend("Gene status", override.aes = list(size = point_size * 2, alpha = 1)),
      alpha = guide_legend("Gene status")
    ) +
    labs(x = xlab, y = ylab)

  # Deal with additional colours and labels
  if (is.null(hvg)) {
    # Get first N genes with highes mean
    sel <- arrange(dat, desc(mean)) %>%
      pull(.data$rowname) %>%
      head(top_n)

    # Add point
    p <- p + geom_point(size = point_size)
  } else {
    # Get first N HVGs
    sel <- head(hvg, top_n)

    # Add point
    p <- p + geom_point(aes_string(color = "HVG"), size = point_size) +
      scale_color_manual(
        labels = paste0(c("Not HVG", "HVG"), ":", table(dat$HVG)),
        values = c("black", "red")
      ) +
      guides(color = guide_legend(override.aes = list(size = point_size * 2, alpha = 1)))
  }

  # Add trend line
  p <- p + stat_function(
    fun = function(x) metadata(var)$trend(x), geom = "line",
    alpha = 1, color = trend_color, size = trend_size
  )

  # Add labels
  aes_repel <- aes_string(label = 'ifelse(rowname %in% sel, rowname, "")')
  p <- p + geom_text_repel(aes_repel,
    alpha = 1, color = text_color,
    box.padding = box.padding, size = text_size,
    max.overlaps = max.overlaps, seed = seed, ...
  )

  # Add plot title
  if (!is.null(title)) p + ggtitle(title) else p
}

#' Plot distribution of approximate silhouette widths
#'
#' This function creates a violin scatter plot to show the distribution of
#' the approximate silhouette width across cells in each cluster.
#'
#' @param object A numeric matrix-like object containing observations in rows
#' and variables in columns.
#' @param clusters A vector of length equal to `ncol(object)`, indicating the
#' cluster assigned to each observation.
#' @param printDiff Logical scalar indicating whether to print the precentage
#' of cells assigned to different cluster(s). Default is `TRUE`.
#' @param plot Logical scalar indicating whether to draw the plot. Default is
#' `TRUE`.
#' @param cluster_color A character vector of color codes indicating the
#' colour of the clusters. Default is `NULL`, and use [choosePalette()] to
#' select a palette.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.8.
#' @param point_shape An integer scalar (between 0 and 25) indicating the
#' shape aesthetics. Default is 16.
#' @param swarm_method A string indicating the method for arranging points.
#' Method recongised by \pkg{ggbeeswarm} are "quasirandom", "pseudorandom",
#' "smiley" and "frowney". Default is "quasirandom".
#' @param add_mean Logical scalar indicating whether add a horizontal
#' reference line at mean silhouette width. Default is `TRUE`.
#' @param mean_color A string indicating the colour of the reference line.
#' Default is "firebrick".
#' @param mean_size A numeric scalar indicating the size of the reference
#' line. Default is 1.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#'
#' @return A `ggplot` object when `plot = TRUE`
#'
#' @details
#' The function uses `approxSilhouette()` from the \pkg{bluster} package' to
#' calculate the silhouette widths using an approximate approach. Using the
#' information stored in the returned DataFrame, it:
#'   - prints the descriptive statistics of approximate silhouette widths.
#'   - calculates the percentage of cells assigned to the closest cluster and
#' prints the result when `printDiff = TRUE`.
#'   - creates a violin scatter plot to show thedistribution of the approximate
#' silhouette width when `plot = TRUE`.
#'
#' These information can allow users to identify poorly separate clusters
#' quickly.
#'
#' @author I-Hsuan Lin
#'
#' @name plotSilhouette
#'
#' @seealso [bluster::approxSilhouette()]
#'
#' @export
#' @import ggplot2
#' @importFrom bluster approxSilhouette
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom rlang inform
#' @importFrom rlang .data
#' @importFrom tidyr complete
#' @importFrom tibble tibble
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom cowplot theme_cowplot
#' @examples
#' # Load demo dataset
#' data(sce)
#'
#' mat <- SingleCellExperiment::reducedDim(sce, "PCA")
#' plotSilhouette(mat, sce$label)
plotSilhouette <- function(object, clusters, printDiff = TRUE, plot = TRUE,
                           cluster_color = NULL, point_size = 2, point_alpha = 0.8,
                           point_shape = 16, swarm_method = "quasirandom", add_mean = TRUE,
                           mean_color = "firebrick", mean_size = 1, theme_size = 18) {
  sil.approx <- approxSilhouette(object, clusters)
  sil.data <- as.data.frame(sil.approx)
  sil.data$closest <- factor(ifelse(sil.data$width > 0, sil.data$cluster, sil.data$other))
  sil.data$isSame <- as.character(sil.data$closest) == as.character(sil.data$cluster)

  inform("Silhouette width summary:")
  print(summary(sil.data$width))

  if (printDiff) {
    grouped <- sil.data %>%
      group_by(.data$cluster, .data$isSame) %>%
      summarise(N = n())

    # Check if 'grouped$isSame' is all TRUE
    if (!FALSE %in% unique(grouped$isSame)) grouped <- rbind(grouped, tibble(cluster = grouped$cluster, isSame = FALSE, N = 0))

    print(grouped %>% mutate("% Different" = round(.data$N / sum(.data$N) * 100, 2)) %>%
      ungroup() %>% complete(.data$cluster, .data$isSame, fill = list(N = 0, "% Different" = 0)) %>%
      dplyr::filter(.data$isSame == FALSE) %>% dplyr::select(-c(.data$N, .data$isSame)), n = Inf)
  }

  if (plot) {
    cluster_color <- choosePalette(clusters, cluster_color)
    mean.width <- mean(sil.data$width)
    aes <- aes_string(x = "cluster", y = "width", colour = "closest")
    p <- ggplot(sil.data, aes) +
      geom_quasirandom(size = point_size, alpha = point_alpha, shape = point_shape, method = swarm_method) +
      geom_hline(yintercept = 0, size = 0.5, color = "black") +
      scale_color_manual(values = cluster_color) +
      theme_cowplot(theme_size) +
      ylab("Silhouette width Si")

    p <- if (add_mean) p + geom_hline(yintercept = mean.width, color = mean_color, size = mean_size, linetype = "dashed") else p

    ncol <- ceiling(length(table(clusters)) / 20) # show 20 clusters in a column

    # Return plot
    p + guides(color = guide_legend(
      title = "Cluster", ncol = ncol,
      override.aes = list(size = point_size * 4, alpha = 1)
    ))
  }
}

#' Create QC plots from `findDoubletClusters` output
#'
#' This function produces QC plots to aid deciding potential doublet clusters
#' after running `findDoubletClusters()` from the \pkg{scDblFinder} package.
#'
#' @param dbl A DataFrame returned by `findDoubletClusters()` that identifies
#' potential clusters of doublet cells.
#' @param clusters A vector of length equal to number of rows in `dbl`,
#' indicating the cluster levels. Default is `NULL` and use the row names in
#' `dbl`.
#' @param cluster_color A character vector of color codes indicating the
#' colour of the clusters. Default is `NULL`, and use [choosePalette()] to
#' select a palette.
#' @param qc_plot A numeric scalar indicating which QC plot to create.
#' Default is 0 to produce a compound figure, or:
#'   - 1: A scatter plot of median number of significant genes (`median.de`)
#' against number of significant genes (`num.de`).
#'   - 2: A barplot showing the proportion of cells in each of the query
#' cluster (`prop`).
#'   - 3: A scatter plot of ratio of the median library sizes for the second
#' source cluster (`lib.size2`) against first source cluster (`lib.size1`).
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param text_size A numeric scalar indicating the size of the label.
#' This is passed to `geom_text_repel()`. Default is 8.
#' @param box.padding A scalar indicating the amount of padding around
#' bounding box, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param point.padding A scalar indicating the amount of padding around
#' labeled point, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param max.overlaps Exclude text labels that overlap too many things.
#' This is passed to `geom_text_repel()`. Default is Inf.
#' @param seed Random seed passed to [set.seed()]. This is passed to
#' `geom_text_repel()`. Default is 12321.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#'
#' @return A `ggplot` object when `qc_plot` is not 0, or a \pkg{ggplot2}
#' plot with an object of class `c("gg", "ggplot")` when `qc_plot = 0`
#'
#' @details NULL
#'
#' @author I-Hsuan Lin
#'
#' @name plotqcDoubletClusters
#'
#' @seealso [scDblFinder::findDoubletClusters()]
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom rlang abort
#' @importFrom rlang .data
#' @importFrom cowplot theme_cowplot
#' @importFrom cowplot plot_grid
#' @examples
#' library(SingleCellExperiment)
#'
#' data(dbl_results)
#'
#' # Create compound figure
#' plotqcDoubletClusters(dbl_results)
#'
#' # Create median.de against num.de
#' plotqcDoubletClusters(dbl_results, qc_plot = 1)
plotqcDoubletClusters <- function(dbl, clusters = NULL, cluster_color = NULL, qc_plot = 0,
                                  point_size = 2, text_size = 6, box.padding = 0.5,
                                  point.padding = 0.5, max.overlaps = Inf, seed = 12321,
                                  theme_size = 18) {
  if (!qc_plot %in% c(0, 1, 2, 3)) abort("Invalid `qc_plot` setting.")
  if (!all(c("num.de", "median.de", "prop", "lib.size1", "lib.size2") %in% colnames(dbl))) abort("Unable to find all required columns in `dbl`.")
  dat <- rownames_to_column(as.data.frame(dbl), var = "Cluster") %>%
    mutate(Cluster = as.factor(.data$Cluster))

  dbl.clu <- row.names(dbl)
  if (is.null(clusters)) {
    if (requireNamespace("gtools")) {
      label <- gtools::mixedsort(dbl.clu)
    } else if (.check_wholenum(dbl.clu)) {
      label <- sort(as.numeric(dbl.clu))
    } else {
      label <- sort(dbl.clu)
    }
  } else {
    if (length(clusters) == length(dbl.clu)) {
      label <- clusters
    } else {
      label <- if (is.factor(clusters) & length(levels(clusters) == length(dbl.clu))) levels(clusters) else abort("Number of clusters in `dbl` and `clusters` are different.")
    }
  }
  dat <- mutate(dat, Cluster = factor(.data$Cluster, levels = label))

  cluster_color <- choosePalette(dat$Cluster, cluster_color)

  if (qc_plot %in% c(0, 1)) {
    aes <- aes_string(x = "num.de", y = "median.de", color = "Cluster")
    p1 <- ggplot(dat, aes) +
      geom_point(size = point_size) +
      guides(color = guide_legend(override.aes = list(size = point_size * 2, alpha = 1))) +
      scale_color_manual(values = cluster_color) +
      theme_cowplot(theme_size) +
      ggtitle("num.de vs. median.de")

    p1 <- p1 + geom_text_repel(aes_string(label = "Cluster"),
      size = text_size,
      box.padding = box.padding, point.padding = point.padding,
      max.overlaps = max.overlaps, seed = seed, show.legend = FALSE
    )

    # Return plot
    if (qc_plot != 0) {
      return(p1)
    }
  }

  if (qc_plot %in% c(0, 2)) {
    aes <- aes_string(x = "Cluster", y = "prop", fill = "Cluster")
    p2 <- ggplot(dat, aes) +
      geom_col(size = 1, width = 0.8) +
      scale_fill_manual(values = cluster_color) +
      theme_cowplot(theme_size) +
      ggtitle("proportion of cells")

    # Return plot
    if (qc_plot != 0) {
      return(p2)
    }
  }

  if (qc_plot %in% c(0, 3)) {
    aes <- aes_string(x = "lib.size1", y = "lib.size2", color = "Cluster")
    p3 <- ggplot(dat, aes) +
      geom_point(size = point_size) +
      guides(color = guide_legend(override.aes = list(size = point_size * 2, alpha = 1))) +
      xlim(0, max(c(dbl$lib.size1, dbl$lib.size2))) +
      ylim(0, max(c(dbl$lib.size1, dbl$lib.size2))) +
      scale_color_manual(values = cluster_color) +
      theme_cowplot(theme_size) +
      ggtitle("lib.size1 vs. lib.size2")

    p3 <- p3 + geom_text_repel(aes_string(label = "Cluster"),
      size = text_size,
      box.padding = box.padding, point.padding = point.padding,
      max.overlaps = max.overlaps, seed = seed, show.legend = FALSE
    )

    # Return plot
    if (qc_plot != 0) {
      return(p3)
    }
  }

  # Return compound figure
  plot_grid(p1, p2, p3, align = "vh", ncol = 2)
}

#' Add labels to reduced dimension plots
#'
#' This function is designed to be used with [scater::plotReducedDim()]
#' from the \pkg{scater} package without specifying text label (i.e. via
#' `text_by`). Then use `add_label()` to place labels centrally.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param dimname A string or integer scalar indicating the reduced dimension
#' result in `reducedDims(sce)` to plot. Default is "TSNE".
#' @param text_by A string indicating the column metadata field with which to
#' add text labels on the plot. Alternatively, a character vector of the same
#' length as `colData(sce)` indicating the labels of each cell. Default is
#' "label".
#' @param text_type A string indicating to add text directly to the plot
#' (`"text"`) or draw a rectangle underneath the text (`"label"`). Default
#' is "text".
#' @param text_size A numeric scalar indicating the size of the label. This
#' is passed to `geom_text_repel()`. Default is 8.
#' @param text_color A string indicating the colour of the label. This is
#' passed to `geom_text_repel()`. Default is "black".
#' @param box.padding A scalar indicating the amount of padding around
#' bounding box, as unit or number. This is passed to `geom_text_repel()`.
#' Default is 0.5.
#' @param max.overlaps Exclude text labels that overlap too many things.
#' This is passed to `geom_text_repel()`. Default is 20.
#'
#' @return A `geom` (geometric object)
#'
#' @details The repel-away-from-center behaviour in `plotReducedDim()` should
#' be fixed in \pkg{scater} v1.23.5.
#'
#' @author I-Hsuan Lin
#'
#' @name add_label
#'
#' @seealso [scater::plotReducedDim()]
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scater retrieveCellInfo
#' @importFrom stats median
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes_string
#' @examples
#' library(scater)
#'
#' # Load demo dataset
#' data(sce)
#'
#' plotReducedDim(sce, "TSNE", colour_by = "label") + add_label(sce, "TSNE")
add_label <- function(sce, dimname = "TSNE", text_by = "label", text_type = "text",
                      text_size = 8, text_color = "black", box.padding = 0.5,
                      max.overlaps = 20) {
  .is.sce(sce)
  .check_dimname(sce, dimname)

  red_dim <- as.matrix(reducedDim(sce, dimname))
  df_to_plot <- data.frame(red_dim[, seq_len(2), drop = FALSE])
  colnames(df_to_plot)[seq_len(2)] <- c("X", "Y")

  if (length(text_by) == 1) {
    text_out <- retrieveCellInfo(sce, text_by, search = "colData")
  } else {
    if (length(text_by) != ncol(sce)) abort(sprintf("`text_by` and `ncol(%s)` have different lengths.", deparse(substitute(sce))))
    text_out <- list(name = "label", value = text_by)
  }

  text_out$val <- as.factor(text_out$val)
  by_text_x <- vapply(split(df_to_plot$X, text_out$val), median, FUN.VALUE = 0)
  by_text_y <- vapply(split(df_to_plot$Y, text_out$val), median, FUN.VALUE = 0)

  if (text_type == "label") {
    geom_label_repel(
      data = data.frame(x = by_text_x, y = by_text_y, label = names(by_text_x)),
      mapping = aes_string(x = "x", y = "y", label = "label"),
      inherit.aes = FALSE, size = text_size, colour = text_color,
      max.overlaps = max.overlaps, force = 0
    )
  } else {
    geom_text_repel(
      data = data.frame(x = by_text_x, y = by_text_y, label = names(by_text_x)),
      mapping = aes_string(x = "x", y = "y", label = "label"),
      inherit.aes = FALSE, size = text_size, colour = text_color,
      max.overlaps = max.overlaps, force = 0
    )
  }
}

#' Colour cells by a feature in a 2-dimension representation
#'
#' This function uses `plotReducedDim()` from the \pkg{scater} package to
#' show cells on a pre-calculated low-dimensional projection (such as UMAP or
#' t-SNE) and colour by a choosen cell-specific feature. It uses `add_label()`
#' to label cells when `text_by` is used.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param feature A string indicating the column name in `colData(sce)`
#' containing the feature to be coloured. Alternatively, a character vector of
#' the same length as `colData(sce)` indicating the feature type of each cell.
#' @param dimname A string or integer scalar indicating the reduced dimension
#' result in `reducedDims(sce)` to plot. Default is "TSNE".
#' @param feat_desc A string that describes the coloured feature. Default is
#' `NULL`.
#' @param feat_color A character vector of colour codes indicating the colours
#' of the features, or a palette function that creates a vector of colours
#' along a colour map. Default is `NULL`, and use [choosePalette()] to
#' select a palette.
#' @param color_breaks One of:
#'   - `NULL` for no breaks
#'   - `waiver()` for the default breaks (the scale limits)
#'   - A character vector of breaks
#'   - A function that takes the limits as input and returns breaks
#'     as output
#' @param exprs_by A string or integer scalar specifying which assay from
#' `sce` to obtain expression values from, for use in point aesthetics. Use
#' `assayNames(sce)` to find all availavle assays in `sce`. Default is
#' "logcounts".
#' @param text_by A string indicating the column metadata field with which
#' to add text labels on the plot. Default is `NULL`.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.5.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#' @param legend_pos The position of legends ("none", "left", "right",
#' "bottom","top". Use "none" to disable plot legend. Default is "right".
#' @param legend_just The anchor point for positioning legend inside plot
#' ("center" or two-element numeric vector) or the justification according
#' to the plot area when positioned outside the plot. Default is "center".
#' @param guides_ncol,guides_nrow An integer scalar indicating the desired
#' number of column and row of discrete legends. Default is NULL.
#' @param guides_barwidth,guides_barheight A numeric or a `grid::unit()`
#' object specifying the width and height of the colourbar. Default is NULL.
#' @param guides_size An integer scalar indicating the desired size of the
#' points in a discrete legend. Default is `point_size * 2`.
#' @param title Plot title. Default is `NULL` and the `dimname` is used as
#' title.
#' @param show_title Logical scalar indicating whether to show plot title.
#' Default is `TRUE`.
#' @param show_subtitle Logical scalar indicating whether to show plot
#' subtitle, Default is `TRUE`.
#' @param other_fields Additional cell-based fields to include in the
#' DataFrame. Default is `list()`.
#' @param ... Other arguments passed on to [add_label()].
#'
#' @return A `ggplot` object
#'
#' @details
#' The function uses `plotReducedDim()` from the \pkg{scater} package to
#' create a base plot and adds aesthetic, such as point colours, legend
#' controls, title and subtitles, to the final figure.
#'
#' @author I-Hsuan Lin
#'
#' @name plotProjection
#'
#' @seealso [scater::plotReducedDim()]
#'
#' @export
#' @import ggplot2
#' @importFrom scater plotReducedDim
#' @examples
#' # Load demo dataset
#' data(sce)
#'
#' # Colour and label cells by cluster
#' plotProjection(sce, "label",
#'   dimname = "TSNE", text_by = "label",
#'   feat_desc = "Cluster"
#' )
#'
#' # Colour cells by DKD62 expression and label by cluster
#' plotProjection(sce, "DKD62",
#'   dimname = "UMAP", text_by = "label",
#'   feat_desc = "DKD62 Expression", guides_barheight = 15
#' )
plotProjection <- function(sce, feature, dimname = "TSNE", feat_desc = NULL, feat_color = NULL,
                           color_breaks = waiver(), exprs_by = "logcounts", text_by = NULL,
                           point_size = 2, point_alpha = 0.5, theme_size = 18,
                           legend_pos = "right", legend_just = "center",
                           guides_ncol = NULL, guides_nrow = NULL,
                           guides_barwidth = NULL, guides_barheight = NULL,
                           guides_size = point_size * 2, title = NULL, show_title = TRUE,
                           show_subtitle = TRUE, other_fields = list(), ...) {
  .is.sce(sce)
  .check_dimname(sce, dimname)

  feature <- .check_feature(sce, feature, exprs_by)
  discrete <- !is.numeric(feature$values)

  if (discrete) {
    color <- choosePalette(feature$values, feat_color)
  } else {
    color <- if (is.null(feat_color)) c("blue", "yellow", "red") else feat_color
  }

  suppressMessages({
    p <- plotReducedDim(sce, dimname,
      colour_by = I(feature$values),
      point_size = point_size, point_alpha = point_alpha,
      theme_size = theme_size, other_fields = other_fields
    ) +
      theme(legend.position = legend_pos, legend.justification = legend_just)

    # add labels
    p <- if (!is.null(text_by)) p + add_label(sce, dimname, text_by = text_by, ...) else p

    # add title & subtitle
    if (show_title) {
      p <- if (!is.null(title)) p + labs(title = title) else p + labs(title = dimname)
    } else {
      p <- p + theme(plot.title = element_blank())
    }

    if (show_subtitle) {
      p <- if (!is.null(feat_desc)) p + labs(subtitle = paste("Coloured by", feat_desc)) else p + labs(subtitle = paste("Coloured by", feature$name))
    }

    # add color
    if (discrete) {
      p <- p + scale_color_manual(values = color, breaks = color_breaks) +
        guides(color = guide_legend(
          title = NULL, ncol = guides_ncol, nrow = guides_nrow,
          override.aes = list(size = guides_size, alpha = 1)
        ))
    } else {
      p <- p + scale_color_gradientn(colours = color, breaks = color_breaks) +
        guides(color = guide_colorbar(
          title = NULL, barwidth = guides_barwidth, barheight = guides_barheight,
          frame.colour = "black", ticks.colour = "black"
        ))
    }
    return(p)
  })
}

#' Colour cells by a feature in two 2-dimension representation side-by-side
#'
#' This function make use of `plotProjection()` to show cells on two
#' pre-calculated low-dimensional projection (such as UMAP and t-SNE) in
#' a compound figure.
#'
#' @param dimnames A character or integer vector of length 2 indicating the
#' reduced dimension result in `reducedDims(sce)` to plot. Default is
#' `c("TSNE", "UMAP")`.
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 1.
#' @param titles A character vector of length 2 indicating the plot title.
#' Default is `NULL` and the `dimnames` is used as titles.
#' @param add_void Logical scalar indicating whether to add a completely empty
#' theme using [ggplot2::theme_void()] in place of the legend. This is
#' applicable when `legend_pos = "none"`. The empty theme is added to the
#' right of the main plot using `plot_grid` Default is `FALSE`.
#' @param rel_widths A numeric vector of relative columns widths. For example,
#' in a two-column grid, `rel_widths = c(2, 1)` would make the first column
#' twice as wide as the second column. This is applicable when `legend_pos`
#' is "right" or "left". Default is `c(15, 1)` to make the joined plot 15 times
#' as wide as the legend column if `legend_pos = "right"`.
#' @param rel_heights A numerical vector of relative rows heights. Works just
#' as `rel_widths` does, but for rows rather than columns. This is applicable
#' when `legend_pos` is "bottom" or "top". Default is `c(15, 1)` to make
#' the joined plot 15 times as tall as the legend row if
#' `legend_pos = "bottom"`.
#' @inheritParams plotProjection
#'
#' @return A \pkg{ggplot2} plot with an object of class `c("gg", "ggplot")`
#'
#' @details
#' The function uses `plotProjection()` to create 2 base plots, then
#' `get_legend()` from the \pkg{cowplot} package to produce a shared legend
#' and `plot_grid()` from \pkg{cowplot} to create a compound figure with
#' legend placed at the desired position as specified by `legend_pos`.
#'
#' @author I-Hsuan Lin
#'
#' @name plotProjections
#'
#' @seealso [plotProjection()], [scater::plotReducedDim()]
#'
#' @export
#' @import ggplot2
#' @importFrom rlang warn
#' @importFrom rlang abort
#' @importFrom cowplot get_legend
#' @importFrom cowplot plot_grid
#' @examples
#' # Load demo dataset
#' data(sce)
#'
#' # Plot TSNE and UMAP side-by-side
#' plotProjections(sce, "label",
#'   dimname = c("TSNE", "UMAP"), text_by = "label",
#'   feat_desc = "Cluster"
#' )
#'
#' # Show DKD62 expression
#' plotProjections(sce, "DKD62",
#'   dimname = c("TSNE", "UMAP"), text_by = "label",
#'   feat_desc = "DKD62 Expression", guides_barwidth = 15
#' )
plotProjections <- function(sce, feature, dimnames = c("TSNE", "UMAP"), feat_desc = NULL,
                            feat_color = NULL, color_breaks = waiver(), exprs_by = "logcounts",
                            text_by = NULL, point_size = 1, point_alpha = 0.5,
                            theme_size = 18, legend_pos = "right", legend_just = "center",
                            guides_ncol = NULL, guides_nrow = NULL,
                            guides_barwidth = NULL, guides_barheight = NULL,
                            guides_size = point_size * 2, titles = NULL, show_title = TRUE,
                            show_subtitle = TRUE, other_fields = list(), add_void = FALSE,
                            rel_widths = c(15, 1), rel_heights = c(15, 1), ...) {
  titles <- if (is.null(titles)) dimnames else titles
  if (length(dimnames) != 2L) abort("Wrong `dimnames` format.")
  if (length(titles) != 2L) abort("Wrong `titles` format.")
  if (length(rel_widths) != 2L) abort("Wrong `rel_widths` format.")
  if (length(rel_heights) != 2L) abort("Wrong `rel_heights` format.")

  p1 <- plotProjection(sce, feature,
    dimname = dimnames[1], feat_desc = feat_desc,
    feat_color = feat_color, color_breaks = color_breaks, exprs_by = exprs_by,
    text_by = text_by, point_size = point_size, point_alpha = point_alpha,
    theme_size = theme_size, legend_pos = "none", legend_just = legend_just,
    guides_ncol = guides_ncol, guides_nrow = guides_nrow,
    guides_barwidth = guides_barwidth, guides_barheight = guides_barheight,
    guides_size = guides_size, title = titles[1], show_title = show_title,
    show_subtitle = show_subtitle, other_fields = other_fields, ...
  )

  p2 <- plotProjection(sce, feature,
    dimname = dimnames[2], feat_desc = feat_desc,
    feat_color = feat_color, color_breaks = color_breaks, exprs_by = exprs_by,
    text_by = text_by, point_size = point_size, point_alpha = point_alpha,
    theme_size = theme_size, legend_pos = "none", legend_just = legend_just,
    guides_ncol = guides_ncol, guides_nrow = guides_nrow,
    guides_barwidth = guides_barwidth, guides_barheight = guides_barheight,
    guides_size = guides_size, title = titles[2], show_title = show_title,
    show_subtitle = show_subtitle, other_fields = other_fields, ...
  )

  if (legend_pos == "none") {
    if (add_void) {
      legend <- ggplot() +
        theme_void()
      legend_pos <- "right"
    } else {
      return(plot_grid(p1, p2, align = "vh", nrow = 1))
    }
  } else {
    discrete <- !is.numeric(.check_feature(sce, feature, exprs_by)$values)
    # add color, use shared legend for both plots
    if (discrete) {
      tmp <- p1 + guides(color = guide_legend(
        title = NULL, ncol = guides_ncol, nrow = guides_nrow,
        override.aes = list(size = guides_size, alpha = 1)
      ))
    } else {
      tmp <- p1 + guides(color = guide_colorbar(
        title = NULL, barwidth = guides_barwidth, barheight = guides_barheight,
        frame.colour = "black", ticks.colour = "black"
      ))
    }
    legend <- get_legend(tmp + theme(legend.position = legend_pos, legend.justification = legend_just))
  }

  if (legend_pos == "right") {
    plot_grid(plot_grid(p1, p2, align = "vh", nrow = 1), legend, nrow = 1, rel_widths = rel_widths)
  } else if (legend_pos == "left") {
    if (rel_widths[1] > rel_widths[2]) {
      warn("Re-order `rel_widths` to place smaller value first for legend's relative width.")
      rel_widths <- sort(rel_widths)
    }
    plot_grid(legend, plot_grid(p1, p2, align = "vh", nrow = 1), nrow = 1, rel_widths = rel_widths)
  } else if (legend_pos == "bottom") {
    plot_grid(plot_grid(p1, p2, align = "vh", nrow = 1), legend, ncol = 1, rel_heights = rel_heights)
  } else {
    if (rel_heights[1] > rel_heights[2]) {
      warn("Re-order `rel_heights` to place smaller value first for legend's relative heights.")
      rel_heights <- sort(rel_heights)
    }
    plot_grid(legend, plot_grid(p1, p2, align = "vh", nrow = 1), ncol = 1, rel_heights = rel_heights)
  }
}

#' Show ligand-receptor gene expression in reduced dimensions
#'
#' This function produces a figure to show the gene expression intensity of a
#' ligand-receptor pair on a pre-calculated low-dimensional projection (such
#' as UMAP or t-SNE).
#'
#' @param sce A `SingleCellExperiment` object.
#' @param dimname A string or integer scalar indicating the reduced dimension
#' result in `reducedDims(sce)` to plot. Default is "TSNE".
#' @param lr_pair A character vector of length 2 containing the ligand and
#' receptor gene symbol.
#' @param lr_desc A character vector of length 2 containing short description
#' to change legend title. Default is `c("Ligand","Receptor")`.
#' @param lr_color A character vector of length 2 containing colour aesthetics.
#' Default is `c("blue","red")`.
#' @param lr_sep A character string to define how the 2 genes terms are
#' separated. Default is "-".
#' @param oneplot Logical scalar indicating whether to overlay expressions
#' and produce a single plot or produces two side-by-side plots. Default is
#' `TRUE`.
#' @param exprs_by A string or integer scalar indicating which assay to obtain
#' expression values from, for use in point aesthetics. Default is "logcounts".
#' @param same_scale Logical scalar indicating whether to use same scale limits
#' on both genes. Default is `TRUE`.
#' @param low_color A string containing the  color code to indicate the colour
#' on the lower-end of the colour scale. Default is "gray90".
#' @param point_size A numeric scalar indicating the size of the points.
#' Default is 2.
#' @param point_alpha A numeric scalar (between 0 and 1) indicating the
#' transparency. Default is 0.5.
#' @param point_shape An integer scalar (between 0 and 25) indicating the
#' shape aesthetics. Also accepts an integer vector length 2 to use different
#' shape on the two plots. Default is 16.
#' @param guides_barwidth,guides_barheight A numeric or a `grid::unit()`
#' object specifying the width and height of the colourbar. Default is NULL.
#' @param text_by A string indicating the column metadata field with which
#' to add text labels on the plot. Default is `NULL`.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#' @param ... Other arguments passed on to [add_label()].
#'
#' @return A `ggplot` object when `oneplot = TRUE` or a \pkg{ggplot2} plot
#' with an object of class `c("gg", "ggplot")` when `oneplot = FALSE`
#'
#' @details
#' The function is based on `plotReducedDim()` from the \pkg{scater} package,
#' and uses `new_scale_colour()` from the \pkg{ggnewscale} package to add an
#' additional `layer` where a second `geom` uses another colour scale to show
#' the expression intensity of a second gene.
#'
#' Even though `plotReducedDimLR()` was designed initially to show ligand-
#' receptor expression, by changing the `lr_desc` and `lr_sep` arguments, one
#' can also use this function to show the expression of two genes, that could
#' be expressed in a mutually exclusive fashion, or specific to certain cell
#' clusters or cell types.
#'
#' If it is too difficult to visualise the expression of two colour scales on
#' the same figure, one can use `oneplot = FALSE` to create a compound figure
#' with 2 sub-plots, each showing their respective colours.
#'
#' @author I-Hsuan Lin
#'
#' @name plotReducedDimLR
#'
#' @seealso[scater::plotReducedDim()], [ggnewscale::new_scale_colour()]
#'
#' @export
#' @import ggplot2
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scater retrieveCellInfo
#' @importFrom stats quantile
#' @importFrom rlang abort
#' @importFrom scales squish
#' @importFrom ggnewscale new_scale_colour
#' @importFrom cowplot theme_cowplot
#' @importFrom cowplot draw_label
#' @importFrom cowplot ggdraw
#' @importFrom cowplot plot_grid
#' @examples
#' # Load demo dataset
#' data(sce)
#'
#' # Plot TSNE on 2 genes
#' plotReducedDimLR(sce, "TSNE", c("VWF23", "VVU55"))
#'
#' # Customise annotation
#' plotReducedDimLR(sce, "TSNE", c("OOX46", "BFP78"),
#'   lr_desc = c("B", "D"), lr_sep = " and ", text_by = "label"
#' )
#'
#' # Show 2 plots
#' plotReducedDimLR(sce, "TSNE", c("OOX46", "BFP78"),
#'   lr_desc = c("Grp B", "Grp D"), lr_sep = " and ", text_by = "label",
#'   guides_barheight = 10, oneplot = FALSE
#' )
plotReducedDimLR <- function(sce, dimname = "TSNE", lr_pair, lr_desc = c("Ligand", "Receptor"),
                             lr_color = c("blue", "red"), lr_sep = "-", oneplot = TRUE,
                             exprs_by = "logcounts", same_scale = TRUE, low_color = "gray90",
                             point_size = 2, point_alpha = 0.5, point_shape = 16,
                             guides_barwidth = NULL, guides_barheight = NULL,
                             text_by = NULL, theme_size = 18, ...) {
  .is.sce(sce)
  .check_dimname(sce, dimname)
  if (length(lr_pair) != 2L) abort("Wrong `lr_pair` format.")
  if (length(lr_color) != 2L) abort("Wrong `lr_color` colour format.")
  if (length(lr_desc) != 2L) abort("Wrong `lr_desc` description format.")

  if (length(point_shape) != 2L) {
    if (length(point_shape) == 1L) {
      point_shape <- c(point_shape, point_shape)
    } else {
      abort("Wrong `point_shape` format.")
    }
  }

  gene1 <- retrieveCellInfo(sce, lr_pair[1], exprs_values = exprs_by)
  gene2 <- retrieveCellInfo(sce, lr_pair[2], exprs_values = exprs_by)

  label1 <- paste0(lr_desc[1], "\n", gene1$name)
  label2 <- paste0(lr_desc[2], "\n", gene2$name)

  df_to_plot <- as.data.frame(reducedDim(sce, dimname))[, 1:2] # first 2 columns
  colnames(df_to_plot) <- c("X", "Y")

  dat <- cbind(df_to_plot, data.frame(L = gene1$value, R = gene2$value))

  # Set up colour scale limits (capped at 5% and 95%)
  limits1 <- if (quantile(gene1$value, 0.05) > 0) c(0, quantile(gene1$value, 0.95)) else c(quantile(gene1$value, 0.05), quantile(gene1$value, 0.95))
  limits2 <- if (quantile(gene2$value, 0.05) > 0) c(0, quantile(gene2$value, 0.95)) else c(quantile(gene2$value, 0.05), quantile(gene2$value, 0.95))

  if (same_scale) {
    limits1 <- limits2 <- c(min(c(limits1, limits2)), max(c(limits1, limits2)))
  }

  # Set up plot
  p <- ggplot(data = dat, mapping = aes_string(x = "X", y = "Y"))
  title <- paste(lr_pair, collapse = lr_sep)

  if (oneplot) {
    p <- p + geom_point(aes_string(color = "L"), shape = point_shape[1], size = point_size, alpha = point_alpha) +
      scale_colour_gradientn(label1,
        colours = c(low_color, lr_color[1]), limits = limits1, oob = squish,
        guide = guide_colorbar(
          order = 1, barwidth = guides_barwidth, barheight = guides_barheight,
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      new_scale_colour() +
      geom_point(aes_string(color = "R"), shape = point_shape[2], size = point_size, alpha = point_alpha) +
      scale_colour_gradientn(label2,
        colours = c(low_color, lr_color[2]), limits = limits2, oob = squish,
        guide = guide_colorbar(
          order = 2, barwidth = guides_barwidth, barheight = guides_barheight,
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      theme_cowplot(theme_size) + labs(
        title = title, x = paste(dimname, "1"),
        y = paste(dimname, "2")
      )

    # Return a single plot
    if (!is.null(text_by)) p + add_label(sce, dimname, text_by = text_by, ...) else p
  } else {
    p1 <- p + geom_point(aes_string(color = "L"), shape = point_shape[1], size = point_size, alpha = point_alpha) +
      scale_colour_gradientn(label1,
        colours = c(low_color, lr_color[1]), limits = limits1, oob = squish,
        guide = guide_colorbar(
          barwidth = guides_barwidth, barheight = guides_barheight,
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      theme_cowplot(theme_size) + labs(x = paste(dimname, "1"), y = paste(dimname, "2"))

    p2 <- p + geom_point(aes_string(color = "R"), shape = point_shape[2], size = point_size, alpha = point_alpha) +
      scale_colour_gradientn(label2,
        colours = c(low_color, lr_color[2]), limits = limits2, oob = squish,
        guide = guide_colorbar(
          barwidth = guides_barwidth, barheight = guides_barheight,
          frame.colour = "black", ticks.colour = "black"
        )
      ) +
      theme_cowplot(theme_size) + labs(x = paste(dimname, "1"), y = paste(dimname, "2"))

    if (!is.null(text_by)) {
      p1 <- p1 + add_label(sce, dimname, text_by = text_by, ...)
      p2 <- p2 + add_label(sce, dimname, text_by = text_by, ...)
    }

    title_theme <- calc_element("plot.title", theme_cowplot())

    title <- ggdraw() + draw_label(title,
      x = 0.05,
      hjust = title_theme$hjust,
      vjust = title_theme$vjust,
      fontface = title_theme$face,
      color = title_theme$colour,
      size = title_theme$size * 1.3,
      lineheight = title_theme$lineheight,
      angle = title_theme$angle
    )

    # Return compound figure
    plot_grid(title,
      plot_grid(p1, p2, align = "vh", nrow = 1),
      ncol = 1, rel_heights = c(0.1, 1)
    )
  }
}

#' Create a boxplot of expression values
#'
#' This function produces a boxplot to show the gene expression intensity
#' for a grouping of cells.
#'
#' @param sce A `SingleCellExperiment` object.
#' @param features A character (or factor) vector of row names, a logical
#' vector, or integer vector of indices specifying rows of `sce` to visualize.
#' @param columns A character vector of col names, a logical vector, or
#' integer vector of indices specifying the columns (i.e. subset of cells)
#' of `sce` to visualize. By default, all columns (cells) are used.
#' @param group_by A character vector of length no more than 2 indicating the
#' field(s) of `colData(sce)` containing the grouping factor, e.g., cell types
#' or clusters, to group cells by.
#' @param color_by A character string of either `"Group"` or `"Detected"`
#' indicating how to colour the boxplot. For the option of "Detected, the
#' boxplot will be colour according to the proportion of cells with detectable
#' expression values. "Default is "Detected".
#' @param box_colors A character vector of colour codes indicating the colours
#' of the cell groups, or a palette function that creates a vector of colours
#' along a colour map. Default is `NULL`.
#' @param detection_limit A numeric scalar indicating the value above which
#' observations are deemed to be expressed. Default is 0.
#' @param max_detected A numeric value indicating the cap on the proportion
#' of detected expression values. Default is Default is `NULL`.
#' @param exprs_by A string or integer scalar specifying which assay from
#' `sce` to obtain expression values from, for use in boxplot aesthetics.
#' Use `assayNames(sce)` to find all availavle assays in `sce`.
#' Default is "logcounts".
#' @param facet_ncol A numeric scalar indicating the number of columns to show
#' in the facet wrap. Default is `NULL`.
#' @param guides_barheight A numeric or a `grid::unit()` object specifying the
#' width and height of the colourbar. Default is NULL.
#' @param x.text_size A numeric scalar indicating the size of x axis labels.
#' Default is `NULL`.
#' @param x.text_angle A numeric scalar indicating the angle of x axis labels.
#' Possible choices are 0, 45, 90. Default is 0.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#'
#' @return A `ggplot` object
#'
#' @details
#' The function creates a boxplot showing the expression of selected features
#' in groups of cells.
#'
#' @author I-Hsuan Lin
#'
#' @name plotBox
#'
#' @export
#' @import ggplot2
#' @importFrom rlang abort
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom cowplot theme_cowplot
#' @examples
#' # Load demo dataset
#' data(sce)
#'
#' # All cells in 1 group
#' plotBox(sce, features = rownames(sce)[10:13])
#'
#' # Group cells by "label", colour by cell groups, and change theme base size
#' plotBox(sce, features = rownames(sce)[10:13], group_by = "label",
#'   color_by = "Group", theme_size = 14)
#'
#' # Show cells of labels A and B only and change colour palette
#' keep <- sce$label %in% c("A","B")
#' plotBox(sce, features = rownames(sce)[10:13], group_by = "label",
#'   box_colors = rainbow(50), columns = keep)
#'
#' # Group cells by "label" and "CellType", showing 2 features per row
#' plotBox(sce, features = rownames(sce)[10:15],
#'   group_by = c("label","CellType"), box_colors = rainbow(50),
#'   facet_ncol = 2, x.text_angle = 90, x.text_size = 14)
plotBox <- function(sce, features, columns = NULL, group_by = NULL, color_by = "Detected",
                    box_colors = NULL, detection_limit = 0, max_detected = NULL,
                    exprs_by = "logcounts", facet_ncol = NULL, guides_barheight = NULL,
		    x.text_size = NULL, x.text_angle = 0, theme_size = 18) {
  .is.sce(sce)
  .check_assayname(sce, exprs_by)

  # Check validity of group_by
  if (!is.null(group_by)) {
    if (!is.character(group_by) || length(group_by) > 2L) {
      abort("Invalid values for 'group_by'.")
    } else if (!all(as.character(group_by) %in% colnames(colData(sce)))) {
      abort("Some 'group_by' columns not in `colData(sce)`.")
    }
  }

  # Subset features (genes)
  if (is.logical(features)) {
    if (length(features) != nrow(sce)) abort("logical features index should be of length nrow(sce).")
    features <- rownames(sce)[features]
  } else if (is.numeric(features)) {
    if (any(features > nrow(sce)) || any(features <= 0) || any(features != round(features))) abort("All features should be round numbers; > 0 and <= nrow(sce).")
    features <- rownames(sce)[features]
  } else if (is.character(features) || is.factor(features)) {
    if (!all(as.character(features) %in% rownames(sce))) abort("Some features not in input sce.")
  }
  # if it's character, preserve the order
  if (is.character(features)) {
    features <- factor(features, levels = features)
  }

  # Subset columns (cells)
  if (!is.null(columns)) {
    if (is.logical(columns)) {
      if (length(columns) != ncol(sce)) abort("logical columns index should be of length ncol(sce).")
      columns <- colnames(sce)[columns]
    } else if (is.numeric(columns)) {
      if (any(columns > ncol(sce)) || any(columns <= 0) || any(columns != round(columns))) abort("All columns should be round numbers; > 0 and <= ncol(sce).")
      columns <- colnames(sce)[columns]
    } else if (is.character(columns) || is.factor(columns)) {
      if (!all(as.character(columns) %in% colnames(sce))) abort("Some columns not in input sce.")
    }
  } else {
    columns <- colnames(sce)
  }

  # Create new object
  new <- sce[as.character(features), as.character(columns), drop = FALSE]
  colData(new) <- droplevels(colData(new))

  # Create expr data.frame
  expr <- as.data.frame(t(assay(new, exprs_by)))

  # Create coldata data.frame
  if (is.null(group_by)) {
    coldata <- data.frame(Group = rep("all", ncol(new))) # group all cells into one group
  } else {
    if (length(group_by) == 2L) {
      # Concatenate values from 2 columns into a new column
      coldata <- data.frame(Group = paste(colData(new)[, group_by[1]], colData(new)[, group_by[2]], sep = " - "),
                            colData(new)[, group_by])
      coldata$Group <- as.factor(coldata$Group)

      # Use original orders on 'Group' if possible
      if(is.factor(colData(new)[, group_by[1]]) & is.factor(colData(new)[, group_by[2]])) {
        new_order <- paste(rep(levels(colData(new)[, group_by[1]]), each = length(levels(colData(new)[, group_by[2]]))), 
			   levels(colData(new)[, group_by[2]]), sep = " - ")
        new_order <- new_order[new_order %in% levels(coldata$Group)]
        coldata$Group <- factor(coldata$Group, levels = new_order)
      } else {
        if (requireNamespace("gtools")) coldata$Group <- factor(coldata$Group, levels = gtools::mixedsort(levels(coldata$Group)))
      }
    } else {
      coldata <- cbind(data.frame(Group = colData(new)[, group_by[1]]),
		       data.frame(colData(new)[, group_by[1], drop = FALSE]))
      if (!is.factor(coldata$Group)) coldata$Group <- as.factor(coldata$Group)
    }
  }

  summarized <- summarizeAssayByGroup(assay(new, exprs_by), ids = coldata$Group, statistics = c("prop.detected"),
                                      threshold = detection_limit)
  prop <- assay(summarized, "prop.detected")
  num <- data.frame(Symbol = rep(features, ncol(prop)),
                    Group = rep(colnames(summarized), each = nrow(prop)),
                    Detected = as.numeric(prop))

  df <- cbind(expr, coldata)

  long <- data.frame(Group = rep(df$Group, length(features)),
                     Symbol = rep(features, each = nrow(df)),
                     Expression = as.numeric(as.matrix(df[, features, drop = FALSE])))
  long$Symbol <-  factor(long$Symbol, levels = features)

  if (!is.null(group_by)) {
    long <- cbind(long, data.frame(mapply(rep, df[, group_by, drop = FALSE], length(features))))
    if (is.factor(df[, group_by[1]])) long[, group_by[1]] <- factor(long[, group_by[1]], levels = levels(df[, group_by[1]]))

    if (length(group_by) == 2L) {
      if (is.factor(df[, group_by[2]])) long[, group_by[2]] <- factor(long[, group_by[2]], levels = levels(df[, group_by[2]]))
    }
  }

  # Combin long and num DataFrames
  df <- merge(long, num, by = c("Group","Symbol"))
  df$Group <- factor(df$Group, levels = levels(coldata$Group))

  # Add number of cells to group labels
  freq <- data.frame(table(df$Group)/length(features))
  levels(df$Group) <- paste0(freq$Var1, " (", freq$Freq, ")")

  # Create plot
  my.aes <- if (color_by == "Detected") aes_string(x = "Group", y = "Expression", color = "Detected", fill = "Detected")
          else aes_string(x = "Group", y = "Expression", color = "Group", fill = "Group")

  p <- ggplot(df, my.aes) + geom_boxplot(outlier.size = 0.5, alpha = 0.3) +
          ylab(exprs_by) + theme_cowplot(theme_size)

  p <- if (is.null(facet_ncol)) p + facet_wrap(~ Symbol, scales = "free_y")
          else p + facet_wrap(~ Symbol, scales = "free_y", ncol = facet_ncol)

  if (color_by == "Detected") {
    if (is.null(box_colors)) {
      p <- p + scale_colour_gradientn(colours = c("blue", "yellow", "red"), limits = c(0, 1)) +
              scale_fill_gradientn(colours = c("blue", "yellow", "red"), limits = c(0, 1))
    } else {
      p <- p + scale_color_gradientn(colours = box_colors, limits = c(0, 1)) +
	      scale_fill_gradientn(colours = box_colors, limits = c(0, 1))
    }
    p <- p + guides(color = guide_colourbar(title = "Proportion\nDetected", barheight = guides_barheight),
                    fill = guide_colourbar(title = "Proportion\nDetected", barheight = guides_barheight))
  } else {
    if (!is.null(box_colors)) p <- p + scale_color_manual(values = box_colors) + scale_fill_manual(values = box_colors)
  }

  if (!is.null(x.text_size)) p <- p + theme(axis.text.x = element_text(size = x.text_size))

  if (x.text_angle == 45) p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  else if (x.text_angle == 90) p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  else p <- p + theme(axis.text.x = element_text(angle = 0))

  if (!is.null(group_by)) p <- p + xlab(paste(group_by, collapse = " + "))

  return(p)
}
