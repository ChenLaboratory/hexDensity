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
plotHexDensity = function(kernel, 
                      colramp =colorRampPalette(viridis::viridis(11)),
                      main=deparse(substitute(kernel)),
                      legend=F,
                      colorcut=seq(0,1,length = 1024)) {
  hexbin::plot(kernel,
               mincnt=-0.0000001, 
               main=main,
               colramp=colramp, 
               legend=legend,
               colorcut=colorcut) 
}
