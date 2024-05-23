# Kernel Density with Hexagon
Features:

* Fast Kernel Density calculation using hexagonal grid and plotting of result.

* Edge correction including Jones-Diggllle algorithm as described in Jones, M.C. (1993) Simple boundary corrections for kernel density estimation. Statistics and Computing 3, 135--146.

* Default bandwidth is the diagonal normal scale bandwidth  
## Demonstration
### Data Preparation
Using bei spatial dataset from spatstat.explore
```
library(spatstat.data)
data=bei
```
### Kernel density
Calculating kernel density using hexagonal grid
```
#specify the x, y vectors
density = hexDensity(x=data$x, y=data$y)
#or just let hexDensity figure it out.
density = hexDensity(data)
```
### Plot result
```
plotHexDensity(density)
```
![Rplot10](https://github.com/ChenLaboratory/hexDensity/assets/99466326/887bdf16-40c1-4970-8cde-44ee589ea5b0)


Comparing to density.ppp by spatstat which use square-grid. (eps is to ensure square instead of rectangular grid)
```
library(ks)
bandwidth = sqrt(diag(Hns.diag(cbind(data$x,data$y)))) #Set bandwidth to be the same default plug-in bandwidth in hexDensity.

library(spatstat.explore)
#eps variable is used to turn the grid square instead of rectangle 
density = density.ppp(data,sigma=bandwidth, eps=diff(range(data$x))/128)
plot.im(density,col=colorRampPalette(viridis::viridis(11)))
```
![Rplot08](https://github.com/ChenLaboratory/hexDensity/assets/99466326/55dbd3dd-058d-4c7a-8b78-a29658c039c1)

Comparison to SpatialKDE package, which can also do hexagonal kernel density but really slow to compute and plot. Selected "bandwidth" and "cell size" values are chosen to best fit with the above examples but may not match perfectly. Note that SpatialKDE does not have option for different bandwidth values in different directions and  does not have edge correction.

```
library(SpatialKDE)
library(dplyr)
library(sp)
library(sf)
library(tmap)
#Prepare data
data <- data.frame(bei) %>%
  st_as_sf(coords = c("x", "y"), dim = "XY") %>%
  st_set_crs(28992) %>%
  select()
cell_size <- 12
band_width <- 160
#Create grid
grid <- data %>%
  create_grid_hexagonal(cell_size = cell_size, side_offset = band_width)
#Calculate KDE
kde <- data %>%
  kde(band_width = band_width, kernel = "quartic", grid = grid)
#Plot
tm_shape(kde) +
  tm_polygons(col = "kde_value",style="cont", palette = "viridis", title = "KDE Estimate",legend.show=FALSE)
```
![Rplot14](https://github.com/ChenLaboratory/hexDensity/assets/99466326/1f3577ad-f7bf-46d1-b66f-a6cfcb18a57f)


### MERFISH dataset
Using spatial transcriptomic dataset. This example use MERFISH mouse hypothalamic preoptic region dataset from [MerfishData package](https://bioconductor.org/packages/release/data/experiment/html/MerfishData.html) on bioconductor. Finding density for cells of 'inhibitory' cell class.

```
library(MerfishData)
spe = MouseHypothalamusMoffitt2018()
```

```
#filter for just Inhibitory cells on z-layer -0.14.
cdat = data.frame(colData(spe),spatialCoords(spe))
cdat = subset(cdat, cell_class!= "Ambiguous",select = -c(cell_id,sample_id,sex,behavior,neuron_cluster_id))
cdat = subset(cdat,z == -0.14)
cdat.inhibitory = subset(cdat, cell_class == ("Inhibitory"))


#hexDensity
density.inhibitory = hexDensity(cdat.inhibitory$x,cdat.inhibitory$y)
plotHexDensity(density.inhibitory)
```
![Rplot09](https://github.com/ChenLaboratory/hexDensity/assets/99466326/e3027494-5d0e-4e6a-940c-5a8994106525)

```
#comparison with density.ppp by spatstat

#need to convert into ppp object
cdat.inhibitory = ppp(cdat.inhibitory$x,cdat.inhibitory$y,window = owin(range(cdat.inhibitory$x),range(cdat.inhibitory$y)))

bandwidth = sqrt(diag(Hns.diag(cbind(cdat.inhibitory$x,cdat.inhibitory$y))))
density.inhibitory = density.ppp(cdat.inhibitory,sigma=bandwidth, eps=diff(range(cdat.inhibitory$x))/128)
plot.im(density.inhibitory, col=colorRampPalette(viridis::viridis(11)))
```
![Rplot11](https://github.com/ChenLaboratory/hexDensity/assets/99466326/281873c4-9a05-426a-a60f-01d8ceac38dd)
