#' Check and lift cell type compositions between 2 CellChat objects
#'
#' This function check the cell type compositions in 2 input CellChat objects, 
#' use the `liftCellChat` function from CellChat to update the slot related to the 
#' cell-cell communication network, including slots `object@net`, `object@netP` 
#' and `object@idents`, so that both objects will have same cell labels.
#'
#' @param x1 A CellChat object containing one dataset.
#' @param x2 A CellChat object containing one dataset.
#' @return A list containing updated CellChat objects.
#'
#' @details NA
#'
#' @author I-Hsuan Lin
#'
#' @name check_and_liftCC
#'
#' @export
check_and_liftCC <- function(x1, x2) {
    id1 <- levels(x1@idents)
    id2 <- levels(x2@idents)
    if(identical(id1, id2)) {
        message("Identical cell groups, lift not required")
    } else {
        if(!.check_pkg("CellChat")) abort("Requires R package 'CellChat'.")
        print(sprintf("1st object: %s", paste(id1, collapse = ", ")))
        print(sprintf("2nd object: %s", paste(id2, collapse = ", ")))
        if(all(id1 %in% id2)) { # id1 is missing some idents
            message("Different cell groups. Run `liftCellChat` on 1st object")
            x1 <- CellChat::liftCellChat(x1, levels(x2@idents))
        } else if(all(id2 %in% id1)) { # id2 is missing some idents
            message("Different cell groups. Run `liftCellChat` on 2nd object")
            x2 <- CellChat::liftCellChat(x2, levels(x1@idents))
        } else {
            # Both have different idents
            message("Different cell groups. Run `liftCellChat` on both objects")
            x1 <- CellChat::liftCellChat(x1, levels(x2@idents))
            x2 <- CellChat::liftCellChat(x2, levels(x1@idents))
        }
    }
    list(x1, x2)
}

#' Visualization of network using heatmap
#'
#' This function produces heatmap that shows the number of interactions, interaction 
#' strength, or communication probabilities of a chosen pathway in a single dataset 
#' (defined by `pathway =`), as well as showing the differential comparison of the 
#' chosen measures between two datasets.
#'
#' When showing differential comparison:
#' * Using the default `heatmap_colors`, the red (or blue) represents increased 
#' (or decreased) signaling in the second dataset compared to the first one.
#' * The box plots on the top show the original values of the selected `measure` 
#' from cell groups from the two datasets in each column (incoming signaling).
#' * The box plots on the right show the original values of the selected `measure` 
#' from cell groups from the two datasets in each row (outgoing signaling).
#'
#' @param obj A CellChat object containing one or merged datasets.
#' @param comparison A numerical or character vector of length 2, to indicate the 
#' datasets for comparison. Default is `c(1, 2)`.
#' @param measure Define the measure to plot, `count` shows the number of interactions; 
#' `weight` shows the total interaction weights (strength); 
#' `prob` shows thecommunication probabilities of a chosen pathway.
#' @param pathway A string to indicate the name of signaling networks in the CellChat 
#' object. This is required when `measure = "prob"`.
#' @param heatmap_colors A character vector to indicate the colours used to generate a 
#' colour palette. The first colour provided is used to create a sequential palette 
#' for single-data CellChat object. When comparing two datasets, a diverging palette 
#' is created to match the max/min values. Default is `c("#4575b4", "#d73027")`.
#' @param heatmap_colors_min A numeric value to set the minimum value used to generate 
#' the colour palette. Default is `NULL` and is calculated from the data.
#' @param heatmap_colors_max A numeric value to set the maximum value used to generate
#' the colour palette. Default is `NULL` and is calculated from the data.
#' @param cluster_rows A logical value to indicate whether to cluster on rows. 
#' Default is `FALSE`.
#' @param cluster_columns A logical value to indicate whether to cluster on columns. 
#' Default is `FALSE`.
#' @param row_names_rot A numeric value to set the rotation of row names.
#' Default is `0`.
#' @param column_names_rot A numeric value to set the rotation of column names.
#' Default is `90`.
#' @param fontsize A numeric value to set the overall font size. Default is `9`.
#' @param cell_fontsize A numeric value to set the font size of the values in each 
#' cell. Default is `9`.
#' @param show_values A logical value to indicate whether to show the values in each 
#' cell. Default is `TRUE`.
#' @param legend_side A string to indicate the side to put heatmap legend. 
#' Possible options are `right`, `left`, `bottom`, and `top`. Default is `bottom`.
#' @param legend_direction A string to indicate the direction of the heatmap legend. 
#' Possible options are `vertical` and `horizontal`. Default is `horizontal`.
#' @param legend_title A string to indicate the heatmap legend title. Default is `NULL` 
#' and the function will generate a suitable title depending on the selected `measure`.
#' @param legend_title_position A string to indicate the position of title relative 
#' to the heatmap legend. Possible options are `leftcenter` and `lefttop` for horizontal 
#' legend. Vertical legend accepts `topleft`, `topcenter`, `leftcenter-rot` and 
#' `lefttop-rot`. Default is `leftcenter`.
#' @param legend_width A `unit` object to set the width of the the whole heatmap legend 
#' body. It is only used for horizontal continous legend. Default is `unit(5, "cm")`.
#' @param legend_height A `unit` object to set the height of the whole heatmap legend 
#' body. It is only used for vertical continous legend. Default is `NULL`.
#' @param right_anno_width A `unit` object to set the width of the right (row) 
#' annotation on the right. Default is `unit(2, "cm")`.
#' @param top_anno_height A `unit` object to set the height of the top (column) 
#' annotations on the top. Default is `unit(2, "cm")`.
#' @param anno_colors A character vector to indicate the colours used to fill the 
#' box plots showing the two datasets, or a single colour used to fill the bar plots 
#' showing a single dataset. Default is `c("#66c2a5", "#fc8d62")`.
#' @param anno_legend_ncol A integer to set the number of columns in the legend grids. 
#' Default is `1`.
#' @param anno_legend_side A string to indicate the side to put annotation legend.
#' Default is `bottom`.
#' @param anno_legend_direction A string to indicate the direction of the annotation 
#' legend. Default is `horizontal`.
#' @param anno_legend_title A string to indicate the annotation legend title. 
#' Default is "Sample".
#' @param anno_legend_title_position A string to indicate the position of title relative 
#' to the annotation legend. Default is `leftcenter`.
#' @param plot_title A string to indicate the plot (column) title. Default is `NULL`
#' and the function will generate a suitable title depending on the selected `measure`.
#' @param title_prefix A string to append to the plot (column) title. Default is `NULL`.
#' @param row_title A string to indicate the row title. 
#' Default is `Sources (outgoing signaling)`.
#' @param draw A logical value to indicate whether to draw the heatmap. 
#' When `draw = FALSE`, the function retuns a grob (`gTree`) object made with package 
#' `grid`. Default is `TRUE`.
#' @param data_only A logical value to indicate whether to return a list containing the 
#' calculated values. This is useful for users to create their own plots. 
#' When `data_only = TRUE`, the function will not plot or return a `gTree` object.
#' Default is `FALSE`.
#'
#' @return A plot appears on currect plotting device if `draw = FALSE`. 
#' When `draw = FALSE`, the function retuns a grob (`gTree`) object made with package 
#' `grid`. When `data_only = TRUE`, the function retuns `list` object containing the 
#' data values.
#'
#' @details
#' The function offers an alternative method to create heatmap from CellChat object. 
#'
#' Compared to `CellChat::netVisual_heatmap` function, it offers more controls to the 
#' plot aesthetic and offers two output types (defined by `draw =`): to plot on device 
#' or provide a `gTree` object. This can be used with other R package such as 'patchwork' 
#' (might require warpping `gTree` object with `wrap_elements` function), 'cowplot' 
#' (with `plot_grid` function) or 'gridExtra' (with `grid.arrange` function) to combine 
#' multiple heatmaps, or with other type of plot object types such as `ggplot` into a 
#' single graphic.
#'
#' For users who will like to create their own plots, one can use `data_only = TRUE` to
#' obtain a list that contains the matrix with the calculated values, and optionally the
#' original values of the selected `measure` from cell groups from the two datasets when 
#' showing differential comparison.
#'
#' @author I-Hsuan Lin
#'
#' @name netVis_heatmap
#'
#' @export
#' @import grid
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap grid.boxplot
#' @importFrom circlize colorRamp2
#' @importFrom rlang abort
netVis_heatmap <- function(obj, comparison = c(1, 2), measure = c("count", "weight", "prob"),
                           pathway = NULL, heatmap_colors = c("#4575b4", "#d73027"),
                           heatmap_colors_min = NULL, heatmap_colors_max = NULL,
                           cluster_rows = FALSE, cluster_columns = FALSE,
                           row_names_rot = 0, column_names_rot = 90,
                           fontsize = 9, cell_fontsize = 9, show_values = TRUE,
                           legend_side = "bottom", legend_direction = "horizontal",
                           legend_title = NULL, legend_title_position = "leftcenter",
                           legend_width = unit(5, "cm"), legend_height = NULL,
                           right_anno_width = unit(2, "cm"), top_anno_height = unit(2, "cm"),
                           anno_colors = c("#66c2a5", "#fc8d62"), anno_legend_ncol = 1,
                           anno_legend_side = "bottom", anno_legend_direction = "horizontal",
                           anno_legend_title = "Sample", anno_legend_title_position = "leftcenter",
                           plot_title = NULL, title_prefix = NULL,
                           row_title = "Sources (outgoing signaling)",
			   draw = TRUE, data_only = FALSE) {
  if(class(obj) != "CellChat") abort("Input is not a CellChat class object.")

  measure <- if(length(measure) != 1) measure[1] else measure # default "count" from "net"
  type <- if("prob" %in% names(obj@net)) "single" else "multi"

  if(type == "multi") {
    sample.names <- names(obj@net)

    if(is.character(comparison)) {
      if(comparison %in% sample.names) {
        c1 <- comparison[1]
        c2 <- comparison[2]
      } else {
        abort("The given comparison names cannot be found.")
      }
    } else {
      if(max(comparison) > length(sample.names)) {
        abort("The given comparison index is out of bounds.")
      } else {
        c1 <- comparison[1]
        c2 <- comparison[2]
      }
    }
  } else {
    sample.names <- measure
  }

  # Create a named colour vector
  if(!is.null(names(anno_colors))) {
    if(sample.names %in% names(anno_colors)) {
      anno_colors <- anno_colors[sample.names]
    } else {
      names(anno_colors) <- NULL
    }
  }

  if(is.null(names(anno_colors))) {
    if(type == "multi") {
      names(anno_colors) <- sample.names
    } else {
      names(anno_colors) <- measure
    }
  }

  # Setup matrix and title
  missing_pathway <- 0
  if(measure == "prob") {
    if(is.null(pathway)) abort("Need to specify a pathway name.")
    if(type == "multi") {
      pathways <- unique(c(obj@netP[[c1]]$pathways, obj@netP[[c2]]$pathways))
      if(!pathway %in% pathways) {
        abort(paste0("Pathway '", pathway, "' cannot be found."))
      }
      pathways <- intersect(obj@netP[[c1]]$pathways, obj@netP[[c2]]$pathways)
      if(pathway %in% pathways) {
        measure1 <- obj@netP[[c1]]$prob[,,pathway]
        measure2 <- obj@netP[[c2]]$prob[,,pathway]
      } else {
        # Missing pathway in 1 of the datasets
        if(!pathway %in% obj@netP[[c1]]$pathways) { # only in 2nd dataset
          measure2 <- obj@netP[[c2]]$prob[,,pathway]
          measure1 <- matrix(0, nrow = nrow(measure2), ncol = ncol(measure2))
          missing_pathway <- 1
      } else { # only in 1st dataset
          measure1 <- obj@netP[[c1]]$prob[,,pathway]
          measure2 <- matrix(0, nrow = nrow(measure1), ncol = ncol(measure1))
          missing_pathway <- 2
        }
      }
      mat <- measure2 - measure1
      plot_title <- if(is.null(plot_title)) paste("Differential communication prob. of",
                                                  pathway, "signaling network")
    } else {
      if(pathway %in% obj@netP$pathways) {
        mat <- obj@netP$prob[,,pathway]
      } else {
        abort(paste0("Pathway '", pathway, "' cannot be found."))
      }
      plot_title <- if(is.null(plot_title)) paste("Communication prob. of", pathway, "signaling network")
    }
  } else {
    if(type == "multi") {
      measure1 <- obj@net[[c1]][[measure]]
      measure2 <- obj@net[[c2]][[measure]]
      mat <- measure2 - measure1
      if(is.null(plot_title)) {
        plot_title <- paste("Differential", ifelse(measure == "count", "number of interactions",
                                                   "interaction weights/strength"))
      }
    } else {
      mat <- obj@net[[measure]]
      if(is.null(plot_title)) {
        plot_title <- ifelse(measure == "count", "Number of interactions",
                             "Interaction weights/strength")
        }
    }
  }
  plot_title <- if(!is.null(title_prefix)) paste(title_prefix, plot_title) else plot_title

  # Return data (no plot)
  if(data_only) {
    if(measure == "prob") {
        data <- list(title = plot_title, measure = measure, pathway = pathway, mat = mat)
    } else {
        data <- list(title = plot_title, measure = measure, mat = mat)
    }
    if(type == "multi") {
      data[[sample.names[1]]] <- measure1
      data[[sample.names[2]]] <- measure2
    }
    return(data)
  }

  # Set up cell_fun
  cell_fun <- NULL
  if(show_values) {
    cell_fun <- function(j, i, x, y, width, height, fill) {
      grid.text(sprintf(cell_format, mat[i, j]), x, y, gp = gpar(fontsize = cell_fontsize))
    }
  }
  cell_format <- if(isTRUE(all.equal(as.numeric(mat), as.integer(as.numeric(mat))))) "%d" else "%.2f"

  # Set up colour palette
  heatmap_colors_min <- if(is.null(heatmap_colors_min)) min(mat) else heatmap_colors_min
  heatmap_colors_max <- if(is.null(heatmap_colors_max)) max(mat) else heatmap_colors_max

  if(type == "multi" && missing_pathway == 0) {
    # Diverging palette
    heatmap_col <- colorRamp2(c(heatmap_colors_min, 0,heatmap_colors_max),
                              c(heatmap_colors[1], "white", heatmap_colors[2]), space = "RGB")
  } else {
    # Sequential palette
    if(missing_pathway == 1) {
      heatmap_col <- colorRamp2(c(heatmap_colors_min, heatmap_colors_max),
                                c("white", heatmap_colors[2]), space = "RGB")
    } else if(missing_pathway == 2) {
      heatmap_col <- colorRamp2(c(heatmap_colors_min, heatmap_colors_max),
                                c(heatmap_colors[1], "white"), space = "RGB")
    } else {
      # Single dataset
      heatmap_col <- colorRamp2(c(heatmap_colors_min, heatmap_colors_max),
                                c("white", heatmap_colors[1]), space = "RGB")
    }
  }

  font_small <- fontsize*0.75
  font_large <- fontsize*1.25

  # Set up dodge boxplot
  if(type == "multi") {
    rg <- range(c(measure1, measure2))
    rg[1] <- rg[1] - (rg[2] - rg[1])* 0.05
    rg[2] <- rg[2] + (rg[2] - rg[1])* 0.05
    col1 <- anno_colors[sample.names[1]]
    col2 <- anno_colors[sample.names[2]]

    right_boxplots <- function(index) {
      nr = length(index)
      pushViewport(viewport(xscale = rg, yscale = c(0.5, nr + 0.5)))
      for(i in seq_along(index)) {
        grid.rect(y = nr-i+1, height = 1, default.units = "native")
        grid.boxplot(measure1[ index[i], ], pos = nr-i+1 + 0.2, box_width = 0.3,
                     gp = gpar(fill = col1), direction = "horizontal")
        grid.boxplot(measure2[ index[i], ], pos = nr-i+1 - 0.2, box_width = 0.3,
                     gp = gpar(fill = col2), direction = "horizontal")
      }
      grid.xaxis(gp = gpar(fontsize = font_small))
      popViewport()
    }

    top_boxplots <- function(index) {
      nr = length(index)
      pushViewport(viewport(yscale = rg, xscale = c(0.5, nr + 0.5)))
      for(i in seq_along(index)) {
        grid.rect(x = i, width = 1, default.units = "native")
        grid.boxplot(measure1[ , index[i] ], pos = i - 0.2, box_width = 0.3,
                     gp = gpar(fill = col1), direction = "vertical")
        grid.boxplot(measure2[ , index[i] ], pos = i + 0.2, box_width = 0.3,
                     gp = gpar(fill = col2), direction = "vertical")
        }
        grid.yaxis(gp = gpar(fontsize = font_small))
        popViewport()
    }

    # Set up heatmap annotations
    ha <- HeatmapAnnotation(measure = top_boxplots, height = top_anno_height,
                            show_annotation_name = FALSE)
    ra <- rowAnnotation(measure = right_boxplots, width = right_anno_width,
                        show_annotation_name = TRUE, annotation_name_offset = top_anno_height/2,
                        annotation_name_side = "top", annotation_name_rot = 0,
                        annotation_name_gp = gpar(fontsize = fontsize))
    ra@anno_list$measure@label <- measure

    legend_title <- if(is.null(legend_title)) bquote(Delta ~ .(measure)) else legend_title
  } else {
    # Set up heatmap annotations (barplots showing total)
    fill <- anno_colors[1]
    ha <- HeatmapAnnotation(measure = anno_barplot(colSums(mat), baseline = 0, gp = gpar(fill = fill)),
                            height = top_anno_height, show_annotation_name = FALSE)
    ra <- rowAnnotation(measure = anno_barplot(rowSums(mat), baseline = 0, gp = gpar(fill = fill)),
                        width = right_anno_width, show_annotation_name = TRUE,
                        annotation_name_offset = top_anno_height/2,
                        annotation_name_side = "top", annotation_name_rot = 0,
                        annotation_name_gp = gpar(fontsize = fontsize))
    ra@anno_list$measure@label <- paste("total\n", measure)

    legend_title <- if(is.null(legend_title)) measure else legend_title
  }

  column_names_centered <- if(column_names_rot == 0) TRUE else FALSE
  row_names_centered <- if(row_names_rot == 90) TRUE else FALSE

  # Set up main heatmap
  ht <- Heatmap(mat, top_annotation = ha, right_annotation = ra,
                cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                row_title = row_title, column_title = plot_title,
                row_title_side = "left", column_title_side = "top",
                row_title_gp = gpar(fontsize = font_large),
                column_title_gp = gpar(fontsize = font_large),
                row_names_rot = row_names_rot, column_names_rot = column_names_rot,
                row_names_side = "left", column_names_side = "bottom",
                row_names_gp = gpar(fontsize = fontsize),
                column_names_gp = gpar(fontsize = fontsize),
                row_names_centered = row_names_centered,
                column_names_centered = column_names_centered,
                col = heatmap_col, na_col = "white",
                rect_gp = gpar(col = "black", lwd = 0.5),
                heatmap_legend_param = list(title = legend_title,
                                            title_gp = gpar(fontsize = fontsize),
                                            title_position = legend_title_position,
                                            labels_gp = gpar(fontsize = fontsize),
                                            direction = legend_direction,
                                            legend_width = legend_width,
                                            legend_height = legend_height),
                cell_fun = cell_fun)

  if(type == "multi") {
    # Set up annotation legend
    annotation_legend_gp_fill <- if(!is.null(anno_colors)) anno_colors else NULL
    lgd = list(Legend(title = anno_legend_title,
                      title_gp = gpar(fontsize = fontsize),
                      title_position = anno_legend_title_position,
                      labels = sample.names,
                      labels_gp = gpar(fontsize = fontsize),
                      legend_gp = gpar(fill = annotation_legend_gp_fill),
                      direction = anno_legend_direction,
                      type = "grid", ncol = anno_legend_ncol))

    # Output
    if(draw) {
      return({
        draw(ht, heatmap_legend_side = legend_side, merge_legend = TRUE,
             annotation_legend_list = lgd, annotation_legend_side = anno_legend_side)
        })
    } else {
        return({
          grid.grabExpr(draw(ht, heatmap_legend_side = legend_side, merge_legend = TRUE,
                             annotation_legend_list = lgd, annotation_legend_side = anno_legend_side))
        })
    }
  } else {
    # Output
    if(draw) {
      return(draw(ht, heatmap_legend_side = legend_side))
    } else {
      return(grid.grabExpr(draw(ht, heatmap_legend_side = legend_side)))
    }
  }
}
