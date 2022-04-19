#' Discrete palettes
#'
#' Colour palettes designed for discrete, categorical data with no more
#' than 40 categories.
#'
#' @param n Number of colours to return. (default: maximum number of
#' available colours)
#'
#' @return A vector of colours as hex strings.
#'
#' @details
#' The `c30()` palette has 30 unique colours and `c40()` palette has 40
#' unique colours. The `c40()` colour palette is taken from `plotScoreHeatmap()`
#' of the \pkg{SingleR} package (which itself is based on \pkg{DittoSeq} and
#' Okabe-Ito colors).
#'
#' @author I-Hsuan Lin
#'
#' @name discrete_colors
#'
#' @seealso [SingleR::plotScoreHeatmap()]
#'
#' @importFrom rlang warn
#' @examples
#' # Show colours as pie charts
#' pie(rep(1, 30), col = c30(), main = "c30 palette")
#' pie(rep(1, 40), col = c40(), main = "c40 palette")
NULL

#' @export
#' @rdname discrete_colors
c30 <- function(n = 30) {
  if (n > 30) {
    warn("Only 30 colours are available with `c30`, returning 30 colours")
    n <- 30
  }
  hex <- c(
    "#006400", "#ff0000", "#0000ff", "#ff8c00", "#800080",
    "#00ffff", "#ff00ff", "#00008b", "#00ff00", "#1e90ff",
    "#8fbc8f", "#483d8b", "#cd853f", "#adff2f", "#b03060",
    "#4682b4", "#ff1493", "#00cd00", "#808000", "#dc143c",
    "#008b8b", "#da70d6", "#cdcd00", "#f08080", "#00ced1",
    "#bf3eff", "#f0e68c", "#90ee90", "#7f0000", "#696969"
  )
  hex[1:n]
}


#' @export
#' @rdname discrete_colors
c40 <- function(n = 40) {
  if (n > 40) {
    warn("Only 40 colours are available with `c40`, returning 40 colours")
    n <- 40
  }
  hex <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
    "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
    "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
    "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57",
    "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
    "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
    "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
  )
  hex[1:n]
}

#' Change `repr.plot.*` behaviour
#'
#' Change the behaviour of `repr.plot.*`.
#'
#' @param width Plotting area width in inches. Default is 10.
#' @param height Plotting area height in inches. Default is 7.
#' @param pointsize Text height in pt. Default is 12.
#' @param bg Background colour. Default is "white".
#' @param antialias Which kind of antialiasing to use for for lines and text?
#' 'gray', 'subpixel' or 'none'? Default is "gray".
#' @param res PPI for rasterization. Defaults to 120.
#' @param quality Quality of JPEG format in %. Default is 90.
#' @param family Font family. "sans", "serif", "mono" or a specific one.
#' Default is "sans".
#'
#' @details
#' The function uses [options()] to change the behaviour of `repr.plot.*`.
#' It provides a quick and easy way to change plot size and other
#' `repr.plot.*` behaviours when running R in a Jupyter Notebook. If use
#' without indicating any argument, the plot behaviour will be reset to
#' default. The `reset.fig()` is an alias of `fig()`.
#'
#' @author I-Hsuan Lin
#'
#' @name fig
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Change plot area width to 8 inches and height to 5 inches
#' fig(width = 8, height = 5)
#' ggplot(mpg, aes(class)) +
#'   geom_bar()
#'
#' # Reset to default settings
#' fig()
#'
#' # Change plot wider and taller
#' fig(width = 14, height = 10)
#' ggplot(mpg, aes(class)) +
#'   geom_bar()
#'
#' # Alias of fig()
#' reset.fig()
#' }
fig <- function(width = 10, height = 7, pointsize = 12, bg = "white", antialias = "gray",
                res = 120, quality = 90, family = "sans") {
  options(
    repr.plot.width = width, repr.plot.height = height,
    repr.plot.pointsize = pointsize, repr.plot.bg = bg,
    repr.plot.antialias = antialias, repr.plot.res = res,
    repr.plot.quality = quality, repr.plot.family = family
  )
}

#' @rdname fig
#' @export
reset.fig <- fig

#' Choose discrete colours
#'
#' Given a character vector of features and optionally a vector of color codes,
#' this function evaluate if the supplied color codes has sufficient number
#' of colours. It returns a named vector of color codes based on the input
#' features, with the same length as the unique features.
#'
#' @param object A character vector containing the discrete or categorical
#' features.
#' @param color A character vector of color codes. Default is to use `c30()`,
#' a palette with 30 colours.
#' @param quiet Logical scalar indicating whether to disable messaging.
#' Default is `TRUE`.
#'
#' @return A named character vector
#'
#' @details
#' By default, it uses the [c30()] palette when no more than 30 colours are
#' required, then the [c40()] palette, and lastly the [rainbow()] colour palette
#' when requiring more than 40 colours.
#'
#' @author I-Hsuan Lin
#'
#' @name choosePalette
#'
#' @seealso [c30()], [c40()], [grDevices::rainbow()]
#'
#' @export
#' @importFrom rlang inform
#' @importFrom rlang warn
#' @importFrom stats na.omit
#' @importFrom stats setNames
#' @importFrom utils flush.console
#' @examples
#' # When input is a vector
#' feat <- LETTERS[1:5]
#' feat
#' choosePalette(feat) # use c30()
#'
#' # When input is a factor of 3 levels
#' feat <- factor(rep(LETTERS[1:3], 5))
#' feat
#' choosePalette(feat, rainbow(10))
choosePalette <- function(object, color = c30(), quiet = TRUE) {
  if (is.factor(object)) {
    features <- levels(object)
  } else {
    if (requireNamespace("gtools")) {
      features <- gtools::mixedsort(unique(object))
    } else {
      features <- sort(unique(as.character(object)))
    }
  }
  n <- length(features)

  if (is.null(color)) {
    color <- .updatePalette(n, quiet = quiet)
  } else {
    ncol <- length(na.omit(color))

    if (n <= ncol) {
      if (is.null(names(color))) {
        if (!quiet) inform(sprintf("Use first '%d' colours from input palette.", n))
        color <- color[1:n]
      } else {
        name <- names(color)
        if (all(features %in% name)) { # found all matching named colours
          if (!quiet) inform(sprintf("Use '%d' matched colours from input palette.", n))
          return(color[features])
        } else {
          warn(sprintf(
            "Found '%d' colour(s) without matched name(s). Updating palette...",
            sum(features %in% name == FALSE)
          ))
          color <- .updatePalette(n)
        }
      }
    } else {
      warn("There are not enough colours in the input palette. Updating palette...")
      color <- .updatePalette(n)
    }
  }
  if (interactive()) flush.console()
  setNames(color, features)
}

#' Add nudged labels in a parallel sets diagram
#'
#' This is the same function as `geom_parallel_sets_labels()` from the
#' \pkg{ggforce} package but with the ability to nudge labels at a fixed
#' distance. It is especially useful when the labels are too long to fit
#' inside the bars depicting the discrete categories.
#'
#' @inheritParams ggforce::geom_parallel_sets_labels
#' @param stat A string indicating the statistical transformation to use on
#' the data for this layer. Default is "parallel_sets_axes".
#' @param position Position adjustment, either as a string, or the result of
#' a call to a position adjustment function. Cannot be jointy specified with
#' `nudge_x` or `nudge_y`.
#' @param parse If `TRUE`, the labels will be parsed into expressions and
#' displayed as described in [grDevices::plotmath()]. Default is `FALSE`.
#' @param nudge_x,nudge_y Horizontal and vertical adjustment to nudge labels
#' by. Useful for offsetting text from points, particularly on discrete
#' scales. Cannot be jointly specified with `position`.
#' @param check_overlap If `TRUE`, text that overlaps previous text in the
#' same layer will not be plotted.
#' @param ... Other arguments passed on to `layer()`.
#'
#' @return A `layer` object
#'
#' @details A pull request of the nudge enhancement has been
#' submitted to the \pkg{ggforce}'s GitHub repository,
#' \url{https://github.com/thomasp85/ggforce/pull/260}, awaiting approval.
#'
#' @author I-Hsuan Lin
#'
#' @name geom_parallel_sets_labs
#'
#' @seealso [ggforce::geom_parallel_sets_labels()], [ggplot2::layer()]
#'
#' @export
#' @import ggforce
#' @importFrom rlang abort
#' @importFrom ggplot2 position_nudge
#' @importFrom ggplot2 layer
#' @importFrom ggplot2 GeomText
#' @examples
#' library(ggforce)
#' data <- as.data.frame(Titanic)
#' data <- gather_set_data(data, 1:4)
#'
#' # Use nudge_x to offset and hjust = 0 to left-justify label
#' ggplot(data, aes(x, id = id, split = y, value = Freq)) +
#'   geom_parallel_sets(aes(fill = Sex), alpha = 0.3, axis.width = 0.1) +
#'   geom_parallel_sets_axes(axis.width = 0.1) +
#'   geom_parallel_sets_labs(colour = "red", angle = 0, nudge_x = 0.1, hjust = 0)
geom_parallel_sets_labs <- function(mapping = NULL, data = NULL, stat = "parallel_sets_axes",
                                    position = "identity", parse = FALSE,
                                    nudge_x = 0, nudge_y = 0, check_overlap = FALSE,
                                    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) abort("You must specify either `position` or `nudge_x`/`nudge_y`.")
    position <- position_nudge(nudge_x, nudge_y)
  }

  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomText,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(parse = parse, check_overlap = check_overlap, na.rm = na.rm, ...)
  )
}

#' Create a 2-variable parallel sets diagrams
#'
#' This function uses the \pkg{ggforce} package to produce a parallel sets
#' diagram for visualising interaction between 2 variables.
#'
#' @param lab1,lab2 Each is a character vector containing membership
#' information.
#' @param labels A character vector of length 2 containing the labels for
#' `lab1` and `lab2`. Default is `c("label1", "label2")`.
#' @param color A character vector of colour codes indicating the colours of
#' the elements in the parallel sets diagram. Default is `NULL`, and use
#' [choosePalette()] to select a palette.
#' @param add_counts Logical scalar indicating whether to show total counts
#' of each element with the labels. Default is `FALSE`.
#' @param add_breakdown A string to control which side to show the detailed
#' breakdown of counts. Allowable character values are `"left"`, `"right"` and
#' `"both"` (same as `add_breakdown = TRUE`). Use `add_breakdown = FALSE` to
#' disable breakdown. Default is `FALSE`.
#' @param text_size A numeric scalar indicating the size of membership labels.
#' Default is 5.
#' @param xlab_size A numeric scalar indicating the size of main categorical
#' label. Default is 14.
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 18.
#'
#' @return A `ggplot` object
#'
#' @details NULL
#'
#' @author I-Hsuan Lin
#'
#' @name plotParallel
#'
#' @seealso [ggforce::geom_parallel_sets()]
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom dplyr pull
#' @importFrom dplyr tally
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#' @importFrom tidyr drop_na
#' @importFrom rlang abort
#' @importFrom rlang .data
#' @importFrom grid unit
#' @importFrom ggforce gather_set_data
#' @importFrom ggforce geom_parallel_sets
#' @importFrom ggforce geom_parallel_sets_axes
#' @examples
#' library(ggforce)
#' data <- as.data.frame(Titanic)
#' plotParallel(data$Class, data$Age, labels = c("class", "age"))
#'
#' # Load demo dataset
#' data(sce)
#'
#' plotParallel(sce$label, sce$CellType,
#'   labels = c("Cluster", "Cell Type"),
#'   add_counts = TRUE
#' )
plotParallel <- function(lab1, lab2, labels = c("label1", "label2"), color = NULL,
                         add_counts = FALSE, add_breakdown = FALSE, text_size = 5,
                         xlab_size = 14, theme_size = 18) {
  .check_is.null(lab1)
  .check_is.null(lab2)
  if (length(lab1) != length(lab2)) abort("`lab1` and `lab2` have different lengths.")

  # Encode vectors as factors, or drop unused levels
  lab1 <- if (!is.factor(lab1)) as.factor(lab1) else droplevels(lab1)
  lab2 <- if (!is.factor(lab2)) as.factor(lab2) else droplevels(lab2)

  lab1name <- levels(lab1)
  lab2name <- levels(lab2)

  # Calculate and append counts to labels
  if (add_counts == TRUE) {
    lab1name <- .count_label(lab1)
    lab2name <- .count_label(lab2)
  }

  # Calculate and append detailed breakdown to labels
  if (add_breakdown == TRUE | add_breakdown %in% c("both", "left", "right")) {
    tab <- data.frame(table(lab1 = lab1, lab2 = lab2))
    breakdown1 <- group_by(tab, .data$lab1) %>%
      summarise_all(paste, collapse = ",") %>%
      pull(.data$Freq)
    breakdown2 <- group_by(tab, .data$lab2) %>%
      summarise_all(paste, collapse = ",") %>%
      pull(.data$Freq)

    lab1name <- if (add_breakdown != "right") paste0(lab1name, "\n[", breakdown1, "]")
    lab2name <- if (add_breakdown != "left") paste0(lab2name, "\n[", breakdown2, "]")
  }

  levels(lab1) <- lab1name
  levels(lab2) <- lab2name

  # Rename labels if duplicate names used in the 2 sets
  if (length(intersect(levels(lab1), levels(lab2))) > 0) {
    levels(lab1) <- paste0("(", labels[1], ") ", levels(lab1))
    levels(lab2) <- paste0("(", labels[2], ") ", levels(lab2))
  }

  # Build data.frame
  data <- data.frame(lab1 = lab1, lab2 = lab2) %>%
    group_by(.data$lab1, .data$lab2) %>%
    tally() %>%
    ungroup() %>%
    drop_na() %>%
    gather_set_data(1:2) %>%
    mutate(x = factor(.data$x, levels = c("lab1", "lab2")))
  data$y <- droplevels(data$y) # or mutate(data, y = fct_drop(y))

  data_labels <- tibble(group = unique(data[, c("x", "y")])$x) %>%
    mutate(
      hjust = ifelse(.data$group == "lab2", 0, 1),
      nudge_x = ifelse(.data$group == "lab2", 0.1, -0.1)
    )

  # Build ggplot object
  aes <- aes_string(x = "x", id = "id", split = "y", value = "n")
  aes_sets <- aes_string(fill = "lab1")
  aes_axes <- aes_string(fill = "y")
  aes_text <- aes_string(y = "n", split = "y")

  p <- ggplot(data, aes) +
    geom_parallel_sets(aes_sets, alpha = 0.6, axis.width = 0.15) + # edge
    geom_parallel_sets_axes(aes_axes, size = 0.3, axis.width = 0.1) + # annotation
    geom_parallel_sets_labs(aes_text,
      hjust = data_labels$hjust, nudge_x = data_labels$nudge_x,
      fontface = "bold", color = "black", size = text_size
    ) +
    scale_x_discrete(labels = labels) +
    theme_void(base_size = theme_size) +
    theme(
      legend.position = "none", plot.margin = unit(c(1, 1, 1, 1), "lines"),
      axis.text.x = element_text(face = "bold", color = "black", size = xlab_size)
    )

  # Return plot
  p + scale_fill_manual(values = choosePalette(data$y, color))
}
