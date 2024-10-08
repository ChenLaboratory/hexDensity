% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xy2hcell.R
\name{xy2hcell}
\alias{xy2hcell}
\title{Find the hexagon cells from xy coordinates given a hexbin object.}
\usage{
xy2hcell(
  hexbin = NULL,
  x,
  y = NULL,
  xbins = NULL,
  xbnds = NULL,
  ybnds = NULL,
  shape = NULL
)
}
\arguments{
\item{hexbin}{hexbin object to be referenced to.}

\item{x, y}{coordinates or vectors of coordinates of the points.}

\item{xbins}{number of bins partitioning the range of xbnds}

\item{xbnds, ybnds}{horizontal and vertical limit of the binning region. Must be numeric vector of length 2.}

\item{shape}{shape = yheight/xwidth of the plotting regions.}
}
\value{
a vector the same length as x with the hexagonal cell ID for each point
}
\description{
Find the hexagon cells IDs from xy coordinates given a hexbin object. Useful 
if you want to get the KDE value at a certain coordinate.
}
\details{
If a hexbin object is not provided, parameters of the binning region (xbins, xbnds, ybnds, shape) can be used instead.
For finding the xy coordinates of the hexagons for a hexbin object, see \link[hexbin]{hcell2xy}.
}
\examples{
library(hexbin)
set.seed(133)
d=hexDensity(x=rnorm(20000),y=rnorm(20000),xbins=50)
#Get KDE value at the coordinate x=0,y=0
loc = xy2hcell(d,x=0,y=0)
d@count[loc]
}
