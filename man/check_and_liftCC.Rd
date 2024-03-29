% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_net.R
\name{check_and_liftCC}
\alias{check_and_liftCC}
\title{Check and lift cell type compositions between 2 CellChat objects}
\usage{
check_and_liftCC(x1, x2)
}
\arguments{
\item{x1}{A CellChat object containing one dataset.}

\item{x2}{A CellChat object containing one dataset.}
}
\value{
A list containing updated CellChat objects.
}
\description{
This function check the cell type compositions in 2 input CellChat objects,
use the \code{liftCellChat} function from CellChat to update the slot related to the
cell-cell communication network, including slots \code{object@net}, \code{object@netP}
and \code{object@idents}, so that both objects will have same cell labels.
}
\details{
NA
}
\author{
I-Hsuan Lin
}
