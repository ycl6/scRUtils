% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_sc.R
\name{plotVarianceVsMean}
\alias{plotVarianceVsMean}
\title{Plot logcounts variance against logcounts mean}
\usage{
plotVarianceVsMean(
  sce,
  top_n = 5,
  point_size = 2,
  point_alpha = 0.8,
  text_size = 4,
  text_color = "black",
  box.padding = 0.5,
  max.overlaps = Inf,
  seed = 12321,
  xlab = "Mean log-counts",
  ylab = "Variance of log-counts",
  title = NULL,
  theme_size = 18,
  ...
)
}
\arguments{
\item{sce}{A \code{SingleCellExperiment} object containing the required QC
metrics returned by \code{addPerCellQC()}/\code{perFeatureQCMetrics()}.}

\item{top_n}{An integer scalar indicating the first \code{top_n} genes with
the highest calculated mean to label. Default is 5.}

\item{point_size}{A numeric scalar indicating the size of the points.
Default is 2.}

\item{point_alpha}{A numeric scalar (between 0 and 1) indicating the
transparency. Default is 0.8.}

\item{text_size}{A numeric scalar indicating the size of the label. This
is passed to \code{geom_text_repel()}. Default is 4.}

\item{text_color}{A string indicating the colour of the label. This is
passed to \code{geom_text_repel()}. Default is "black".}

\item{box.padding}{A scalar indicating the amount of padding around
bounding box, as unit or number. This is passed to \code{geom_text_repel()}.
Default is 0.5.}

\item{max.overlaps}{Exclude text labels that overlap too many things.
This is passed to \code{geom_text_repel()}. Default is Inf.}

\item{seed}{Random seed passed to \code{\link[=set.seed]{set.seed()}}. This is passed to
\code{geom_text_repel()}. Default is 12321.}

\item{xlab}{The title of the x-axis. Default is "Mean log-counts".}

\item{ylab}{The title of the y-axis. Default is "Variance of log-counts".}

\item{title}{Plot title. Default is \code{NULL}.}

\item{theme_size}{A numeric scalar indicating the base font size.
Default is 18.}

\item{...}{Other arguments passed on to \code{geom_text_repel()}.}
}
\value{
A \code{ggplot} object
}
\description{
This function calculates the per-feature mean and variance using normalised
counts (accessible via \code{logcounts()}) from a \code{SingleCellExperiment} object,
and returns a logcounts mean-variance scatter plot.
}
\details{
When the \code{rowData} slot of the \code{SingleCellExperiment} input
contains QC metrics for features that were returned by \code{addPerCellQC()} or
\code{perFeatureQCMetrics()}, this function will use \code{detected} (percentage of
expressed features above the detection limit) to calculate \code{pct_dropout}
(percentage of dropouts) and colour the points accordingly.
}
\examples{
library(SingleCellExperiment)

# Load demo dataset
data(sce)

plotVarianceVsMean(sce, title = "logcounts mean-variance plot")
}
\seealso{
\code{\link[scater:plotRowData]{scater::plotRowData()}}, \code{\link[scuttle:logNormCounts]{scuttle::logNormCounts()}},
\code{\link[scuttle:calculateCPM]{scuttle::calculateCPM()}}, \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_text_repel()}}
}
\author{
I-Hsuan Lin
}
