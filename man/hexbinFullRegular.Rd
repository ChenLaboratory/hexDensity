% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hexbinFullRegular.R
\name{hexbinFullRegular}
\alias{hexbinFullRegular}
\title{Binning with regular hexagon}
\usage{
hexbinFullRegular(
  x,
  y = NULL,
  xbins = 128,
  xbnds = range(x),
  ybnds = range(y),
  xlab = NULL,
  ylab = NULL,
  IDs = FALSE,
  weight = NULL
)
}
\arguments{
\item{x, y}{Coords of the points or a single plotting structure to be used in binning. See xy.coords.}

\item{xbins}{Number of bins in a row.}

\item{xbnds, ybnds}{Horizontal and vertical limits of the binning region in x or y units respectively, must encompass range(x) or range(y) respectively; Vector of length 2}

\item{xlab, ylab}{Optional character strings used as labels for x and y. If NULL, sensible defaults are used.}

\item{IDs}{Logical indicating if the individual cell “IDs” should be returned, see hexbin.}
}
\value{
hexbin object
}
\description{
Adapted from hexbin to give a full grid (include hexagon with no count) of regular hexagon. See hexbin for more details.
}
