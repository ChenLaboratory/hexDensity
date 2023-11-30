# library(hexbin)
# library(RColorBrewer)

#' Title
#'
#' @param kernel 
#' @param colramp 
#' @param main 
#' @param legend 
#' @param colorcut 
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
