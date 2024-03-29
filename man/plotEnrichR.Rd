% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_markers.R
\name{plotEnrichR}
\alias{plotEnrichR}
\title{Visualise Enrichr results as barplots}
\usage{
plotEnrichR(
  object,
  db,
  showTerms = 20,
  numChar = 50,
  y = "Count",
  order_by = "P.value",
  theme_size = 16,
  prefix = NULL,
  dir_path = NULL,
  ...
)
}
\arguments{
\item{object}{A named list of list of DataFrames storing outputs returned
by \code{runEnrichR()}, each of which contains results for the corresponding
group of comparison.}

\item{db}{A string indicating the name of one gene-set library to plot.}

\item{showTerms}{An integer scalar indicating the number of terms to show.
Default is 20.}

\item{numChar}{An integer scalar indicating the number characters to keep
in the terms' descriptions. Default is 50.}

\item{y}{A string indicating the variable that should be mapped to the
y-axis. It can be \code{"Count"} and \code{"Ratio"}. Default is "Count".}

\item{order_by}{A string indicating how to order the Enrichr terms before
selecting the first \code{showTerms} terms to plot. It can be \code{"P.value"} or
\code{"Combined.Score"}. Default is "P.value".}

\item{theme_size}{A numeric scalar indicating the base font size.
Default is 16.}

\item{prefix}{A string indicating the prefix of output file. When
\code{prefix = NULL}, the plots are shown in the current graphics device using
\code{\link[grid:grid.draw]{grid::grid.draw()}}. The plots are saved to PDF files using \code{ggsave()}
when \code{prefix} is not \code{NULL}. Default is \code{NULL}.}

\item{dir_path}{The directory path for exported PDF. The files are exported
to the current working directory when \code{file_path = NULL}. Default is \code{NULL}.}

\item{...}{Other arguments passed on to \code{ggsave()}.}
}
\value{
None
}
\description{
This function takes a list of list of DataFrames containing Enrichr results
returned by \code{runEnrichR()} and create barplots from a selected gene-set
library.
}
\details{
The function prints barplots from each group of comparison using
\code{\link[grid:grid.draw]{grid::grid.draw()}} one after another, therefore the plotting feature is
suitable only in Jupyter Notebook or R Markdown. Otherwise, disable plotting
by specifying the \code{prefix} argument to save barplots to PDF files.
}
\examples{
\dontrun{
# Load demo dataset
data(res_deseq2)

# We construct a list of 3 DataFrames using the same DESeq2 output (`res_deseq2`)
# to mimic having a list containing results from multiple sets of comparisons.
res.de <- list(A_B = res_deseq2, A_C = res_deseq2, B_C = res_deseq2)

# Select gene-set libraries
dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018")

# Run enrichR using the D. melanogaster specific modEnrichr site, and
# specify the gene symbols are stored in the 'external_gene_name' column
res.ora <- runEnrichR(res.de,
  dbs = dbs, site = "FlyEnrichr",
  column_by = "external_gene_name"
)

# Print plots on to standard output
plotEnrichR(res.ora, db = "GO_Biological_Process_2018", theme_size = 14)

# Save plots to PDF files in the per-session temporary directory
plotEnrichR(res.ora,
  db = "GO_Biological_Process_2018", prefix = "Enrichr",
  dir_path = tempdir(), width = 12, height = 5
)
}
}
\seealso{
\code{\link[=runEnrichR]{runEnrichR()}}, \code{\link[ggplot2:ggsave]{ggplot2::ggsave()}}
}
\author{
I-Hsuan Lin
}
