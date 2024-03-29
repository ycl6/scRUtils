% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_sc.R
\name{plotBox}
\alias{plotBox}
\title{Create a boxplot of expression values}
\usage{
plotBox(
  sce,
  features,
  columns = NULL,
  group_by = NULL,
  color_by = "Detected",
  box_colors = NULL,
  detection_limit = 0,
  max_detected = NULL,
  exprs_by = "logcounts",
  facet_ncol = NULL,
  guides_barheight = NULL,
  x.text_size = NULL,
  x.text_angle = 0,
  theme_size = 18
)
}
\arguments{
\item{sce}{A \code{SingleCellExperiment} object.}

\item{features}{A character (or factor) vector of row names, a logical
vector, or integer vector of indices specifying rows of \code{sce} to visualize.}

\item{columns}{A character vector of col names, a logical vector, or
integer vector of indices specifying the columns (i.e. subset of cells)
of \code{sce} to visualize. By default, all columns (cells) are used.}

\item{group_by}{A character vector of length no more than 2 indicating the
field(s) of \code{colData(sce)} containing the grouping factor, e.g., cell types
or clusters, to group cells by.}

\item{color_by}{A character string of either \code{"Group"} or \code{"Detected"}
indicating how to colour the boxplot. For the option of "Detected, the
boxplot will be colour according to the proportion of cells with detectable
expression values. "Default is "Detected".}

\item{box_colors}{A character vector of colour codes indicating the colours
of the cell groups, or a palette function that creates a vector of colours
along a colour map. Default is \code{NULL}.}

\item{detection_limit}{A numeric scalar indicating the value above which
observations are deemed to be expressed. Default is 0.}

\item{max_detected}{A numeric value indicating the cap on the proportion
of detected expression values. Default is Default is \code{NULL}.}

\item{exprs_by}{A string or integer scalar specifying which assay from
\code{sce} to obtain expression values from, for use in boxplot aesthetics.
Use \code{assayNames(sce)} to find all availavle assays in \code{sce}.
Default is "logcounts".}

\item{facet_ncol}{A numeric scalar indicating the number of columns to show
in the facet wrap. Default is \code{NULL}.}

\item{guides_barheight}{A numeric or a \code{grid::unit()} object specifying the
width and height of the colourbar. Default is NULL.}

\item{x.text_size}{A numeric scalar indicating the size of x axis labels.
Default is \code{NULL}.}

\item{x.text_angle}{A numeric scalar indicating the angle of x axis labels.
Possible choices are 0, 45, 90. Default is 0.}

\item{theme_size}{A numeric scalar indicating the base font size.
Default is 18.}
}
\value{
A \code{ggplot} object
}
\description{
This function produces a boxplot to show the gene expression intensity
for a grouping of cells.
}
\details{
The function creates a boxplot showing the expression of selected features
in groups of cells.
}
\examples{
# Load demo dataset
data(sce)

# All cells in 1 group
plotBox(sce, features = rownames(sce)[10:13])

# Group cells by "label", colour by cell groups, and change theme base size
plotBox(sce, features = rownames(sce)[10:13], group_by = "label",
  color_by = "Group", theme_size = 14)

# Show cells of labels A and B only and change colour palette
keep <- sce$label \%in\% c("A","B")
plotBox(sce, features = rownames(sce)[10:13], group_by = "label",
  box_colors = rainbow(50), columns = keep)

# Group cells by "label" and "CellType", showing 2 features per row
plotBox(sce, features = rownames(sce)[10:15],
  group_by = c("label","CellType"), box_colors = rainbow(50),
  facet_ncol = 2, x.text_angle = 90, x.text_size = 14)
}
\author{
I-Hsuan Lin
}
