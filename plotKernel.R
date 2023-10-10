library(hexbin)
library(RColorBrewer)

plotKernel = function(kernel, 
                      colramp =colorRampPalette(rev(brewer.pal(11,'Spectral'))),
                      legend=F,
                      colorcut=seq(0,1,length = 2*kernel@ncells)) {
  hexbin::plot(kernel,mincnt=-0.0000001, colramp=colramp, legend=legend,colorcut=colorcut) 
}