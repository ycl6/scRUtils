% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_basics.R
\name{plotParallel}
\alias{plotParallel}
\title{Create a 2-variable parallel sets diagrams}
\usage{
plotParallel(
  lab1,
  lab2,
  labels = c("label1", "label2"),
  color = NULL,
  add_counts = FALSE,
  add_breakdown = FALSE,
  text_size = 5,
  xlab_size = 14,
  theme_size = 18
)
}
\arguments{
\item{lab1, lab2}{Each is a character vector containing membership
information.}

\item{labels}{A character vector of length 2 containing the labels for
\code{lab1} and \code{lab2}. Default is \code{c("label1", "label2")}.}

\item{color}{A character vector of colour codes indicating the colours of
the elements in the parallel sets diagram. Default is \code{NULL}, and use
\code{\link[=choosePalette]{choosePalette()}} to select a palette.}

\item{add_counts}{Logical scalar indicating whether to show total counts
of each element with the labels. Default is \code{FALSE}.}

\item{add_breakdown}{A string to control which side to show the detailed
breakdown of counts. Allowable character values are \code{"left"}, \code{"right"} and
\code{"both"} (same as \code{add_breakdown = TRUE}). Use \code{add_breakdown = FALSE} to
disable breakdown. Default is \code{FALSE}.}

\item{text_size}{A numeric scalar indicating the size of membership labels.
Default is 5.}

\item{xlab_size}{A numeric scalar indicating the size of main categorical
label. Default is 14.}

\item{theme_size}{A numeric scalar indicating the base font size.
Default is 18.}
}
\value{
A \code{ggplot} object
}
\description{
This function uses the \pkg{ggforce} package to produce a parallel sets
diagram for visualising interaction between 2 variables.
}
\examples{
library(ggforce)
data <- as.data.frame(Titanic)
plotParallel(data$Class, data$Age, labels = c("class", "age"))

# Load demo dataset
data(sce)

plotParallel(sce$label, sce$CellType,
  labels = c("Cluster", "Cell Type"),
  add_counts = TRUE
)
}
\seealso{
\code{\link[ggforce:geom_parallel_sets]{ggforce::geom_parallel_sets()}}
}
\author{
I-Hsuan Lin
}
