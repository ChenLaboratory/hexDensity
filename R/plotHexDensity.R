#' Plot method for hexDensity
#'
#' Adapted plotting function for hexbin object. 
#' 
#' @param hexDensity hexbin object returned by hexDensity
#' @param colramp Color function that accept an integer n and return n colors.
#' @param main Main title
#' @param legend Legend is currently non-functional and should be ignored.
#' @param colorcut Vector of values covering [0, 1] that determine hexagon color class boundaries and hexagon legend size boundaries. Alternatively, an integer (<= maxcnt) specifying the number of equispaced colorcut values in [0,1]. See hexbin.
#'
#' @return
#' @export
#'
#' @examples
plotHexDensity = function(hexDensity, 
                      colramp =colorRampPalette(viridis::viridis(11)),
                      main=deparse(substitute(hexDensity)),
                      legend=F,
                      colorcut=seq(0,1,length = 1024)) {
  hexbin::plot(hexDensity,
               mincnt=-0.0000001, 
               main=main,
               colramp=colramp, 
               legend=legend,
               colorcut=colorcut) 
}
