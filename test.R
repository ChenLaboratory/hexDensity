source("hexDensity.R")
source("plotKernel.R")
library(RColorBrewer)
library(spatstat.data)
library(plotly)

###preparing data
library(hexbin)
my_colors=colorRampPalette(viridis::viridis(11))
data = bei
nbins = 128

plot(bei)
###Testing hexDensity
#normal hexbin
bins = hexbin(data$x, data$y,xbins=nbins,shape = 1)
plot(bins)
smoothbin = smooth.hexbin(bins,wts=c(48,4,1))
smoothbin2=smooth.hexbin(smoothbin,wts=c(48,4,1))
smoothbin3=smooth.hexbin(smoothbin2,wts=c(48,4,1))
smoothbin4=smooth.hexbin(smoothbin3,wts=c(48,4,1))
smoothbin5=smooth.hexbin(smoothbin4,wts=c(48,4,1))
plot(smoothbin,main="",colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F)
plot(smoothbin2,main="",colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F)
plot(smoothbin3,main="",colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F)
plot(smoothbin4,main="",colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F)
plot(smoothbin5,main="",colramp=my_colors,colorcut=seq(0,1,length=3000),legend=F)
plot(bins, main="" , colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F) 
hexbinplot(bins, main="" , colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F) 
gplot.hexbin(bins, main="1" , colramp=my_colors,colorcut=seq(0,1,length = 3000),legend=F) 

#hexDensity w increasing sigma
hex = hexDensity(data,xbins=nbins,sigma=40)
print(hex)
plotKernel(hexDensity(data,xbins=nbins,sigma=2)) 
plotKernel(hexDensity(data,xbins=nbins,sigma=10)) 
plotKernel(hexDensity(data,xbins=nbins,sigma=40)) 
plotKernel(hexDensity(data,xbins=nbins,sigma=70,edge=TRUE,diggle=FALSE)) 

raw=second.moment.calc(data, sigma=40, what="smooth",
                   kernel="gaussian",
                   scalekernel=is.character(kernel),
                   weights=NULL, varcov=NULL)
both=second.moment.calc(data, sigma=40, what="smoothedge",
                        kernel="gaussian",
                        scalekernel=is.character(kernel),
                        weights=NULL, varcov=NULL)

#density.ppp comparison
library(spatstat.explore)
square=density.ppp(data,sigma=40)
plot(square)
plot(density.ppp(data,diggle = TRUE))
plot(density.ppp(data,sigma=10, eps=diff(range(data$x))/128) ,col=colorRampPalette(viridis::viridis(11)))
plot(density.ppp(fractional=FALSE, data,sigma=30, eps=diff(range(data$x))/128),col = colorRampPalette(viridis::viridis(11)))

###Plot_ly volcano plot
fig = plot_ly(z=both$edge$v) %>% add_surface()
###Testing with bio data
library(MerfishData)
spe = MouseHypothalamusMoffitt2018()
cdat = data.frame(colData(spe),spatialCoords(spe))
cdat = subset(cdat, cell_class!= "Ambiguous",select = -c(cell_id,sample_id,sex,behavior,neuron_cluster_id))
cdat = subset(cdat,z == -0.14)
cdat.inhibitory = subset(cdat, cell_class == ("Inhibitory"))
cdat.ependymal = subset(cdat, cell_class == ("Ependymal"))
#spatstat compare
cdat.inhibitory.ppp = ppp(cdat.inhibitory$x,cdat.inhibitory$y,window = owin(range(cdat.inhibitory$x),range(cdat.inhibitory$y)))
cdat.ependymal.ppp = ppp(cdat.ependymal$x,cdat.ependymal$y,window = owin(range(cdat.ependymal$x),range(cdat.ependymal$y)))
plot.im(density.ppp(cdat.inhibitory.ppp,sigma=70, eps=diff(range(cdat.inhibitory.ppp$x))/128),col=colorRampPalette(viridis::viridis(11)))
plot.im(density.ppp(fractional=TRUE,cdat.ependymal.ppp,sigma=20, eps=diff(range(cdat.ependymal.ppp$x))/128),col=colorRampPalette(viridis::viridis(11)))
a=density.ppp(fractional=TRUE,cdat.inhibitory.ppp,sigma=20, eps=diff(range(cdat.inhibitory.ppp$x))/128)
#hexDensity
plotKernel(hexDensity(cdat.inhibitory.ppp,sigma=70))
plotKernel(hexDensity(cdat.ependymal.ppp,sigma=30))




#Miscellaneous
nr=128
nc=128
sigma = 100
xcol.ker = c(0:(nc-1),-(nc:1))
yrow.ker = c(0:(nr-1),-(nr:1))
densX.ker <- dnorm(xcol.ker, sd=sigma)
densY.ker <- dnorm(yrow.ker, sd=sigma)
Kern <- outer(densY.ker, densX.ker, "*")
print(sum(Kern))
Kern = Kern/sum(Kern)
dyn.load("hbin.so")
ans <- .Fortran("hbin",
                x = 10,
                y = 10,
                cell = 10,
                cnt = 10,
                xcm = 10,
                ycm = 10,
                xbins = 10,
                shape = 10,
                xbnds = 10,
                ybnds = 10,
                dim = c(10,10),
                n = 10,
                cID = -1)[-(1:2)]

library(SpatialKDE)
library(dplyr)
library(sp)
library(sf)
library(tmap)
data(meuse)

meuse <- data.frame(bei) %>%
  st_as_sf(coords = c("x", "y"), dim = "XY") %>%
  st_set_crs(28992) %>%
  select()

cell_size <- 8
band_width <- 35

grid_meuse <- meuse %>%
  create_grid_hexagonal(cell_size = cell_size, side_offset = band_width)
kde <- meuse %>%
  kde(band_width = band_width, kernel = "quartic", grid = grid_meuse)

tm_shape(kde) +
  tm_polygons(col = "kde_value",style="cont", palette = "viridis", title = "KDE Estimate",legend.show=FALSE) #+
  #tm_shape(meuse) +
  #tm_bubbles(size = 0.1, col = "red")

#check value
library(spatstat.data)
library(spatstat.explore)
data=bei
sigma=70
square=density.ppp(data,sigma=sigma,edge=FALSE)

plot(square)
x=0.1#square$xcol[1] #+ square$xstep
y=0.1#square$yrow[1]
sum=0

for (i in 1:3604) {
  sum = sum + (1/(2*pi*(sigma**2)))*exp(-(((data$x[i]-x)**2) +(data$y[i]-y)**2)/(2*(sigma**2)))
}

hex=hexDensity(data,xbins=128,sigma=sigma,edge=FALSE)
plotKernel(hex)

