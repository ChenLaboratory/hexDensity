library(hexbin)
# library(RColorBrewer)

plotKernel = function(kernel, 
                      colramp =colorRampPalette(viridis::viridis(11)),
                      main=deparse(substitute(kernel)),
                      legend=F,
                      colorcut=seq(0,1,length = 2*kernel@ncells)) {
  hexbin::plot(kernel,
               mincnt=-0.0000001, 
               main=main,
               colramp=colramp, 
               legend=legend,
               colorcut=colorcut) 
}
