% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHexDensity.R
\name{plotHexDensity}
\alias{plotHexDensity}
\title{Plot method for hexDensity}
\usage{
plotHexDensity(
  hexDensity,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xaxt = TRUE,
  yaxt = TRUE,
  lcex = 1,
  colramp = colorRampPalette(viridis(11)),
  colorcut = 1024,
  legend = T,
  legendWidth = 0.05,
  legendDistance = 0.15,
  aspectRatio = diff(hexDensity@xbnds)/diff(hexDensity@ybnds),
  margin = 0.18,
  newpage = T
)
}
\arguments{
\item{hexDensity}{hexbin object returned by hexDensity}

\item{main}{Main title}

\item{xlab, ylab}{x-axis and y-axis label}

\item{xaxt, yaxt}{logical for whether to plot x-axis and y-axis values}

\item{lcex}{Expansion factor for all letters.}

\item{colramp}{Color function that accept an integer n and return n colors.}

\item{colorcut}{An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1].}

\item{legend}{Legend is currently non-functional and should be ignored.}

\item{legendWidth}{Expansion factor for legend width.}

\item{legendDistance}{Expansion factor for the space between the plot and the legend.s}

\item{aspectRatio}{width to height ratio of the plot. Default is the (inverse of) shape value of hexDensity.}

\item{margin}{Minimum guaranteed margin for the plot. Different aspect ratio between the screen and the plot means that margin can be larger on certain sides.}

\item{newpage}{logical for whether to plot on a new page.}
}
\description{
Adapted plotting function for hexbin object.
}
\examples{
set.seed(133)
d = hexDensity(x=rnorm(200),y=rnorm(200),bandwidth=0.15)
plotHexDensity(d)
}
