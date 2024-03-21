# Kernel Density with Hexagon
Features:

* Fast Kernel Density calculation using hexagonal grid and plotting of result.

* Work with SpatialExperiment class from Bioconductor for spatial transcriptomic data.

* Edge correction including Jones-Diggllle algorithm as described in Jones, M.C. (1993) Simple boundary corrections for kernel density estimation. Statistics and Computing 3, 135--146.

* Bandwidth scales based on the coordinate data the same as the density calculated by the spatstat.explore package.
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
density = hexDensity(x=data$x, y=data$y,sigma=25)
#or just let hexDensity figure it out.
density = hexDensity(data,sigma=25)
```
### Plot result
```
plotHexDensity(density)
```
![Rplot02](https://github.com/ChenLaboratory/hexDensity/assets/99466326/5b705995-c6c6-408f-b85d-88430609c851)

Comparing to density.ppp by spatstat which use square-grid. (eps is to ensure square instead of rectangular grid)
```
library(spatstat.explore)
#eps variable is used to turn the grid square instead of rectangle 
density = density.ppp(data,sigma=25, eps=diff(range(data$x))/128)
plot.im(density,col=colorRampPalette(viridis::viridis(11)))
```

![Rplot04](https://github.com/ChenLaboratory/hexDensity/assets/99466326/de5ef328-1239-4b54-ae9f-72e656164ca0)

Comparison to SpatialKDE package which can also do hexagonal kernel density but really slow to compute and plot. Selected "bandwidth" and "cell size" values are chosen to best fit with the above examples but may not match perfectly. Note that SpatialKDE does not have option for Gaussian kernel or edge correction.

```
library(SpatialKDE)
library(dplyr)
library(sp)
library(sf)
library(tmap)
#Prepare data
bei <- data.frame(bei) %>%
  st_as_sf(coords = c("x", "y"), dim = "XY") %>%
  st_set_crs(28992) %>%
  select()
cell_size <- 8
band_width <- 35
#Create grid
grid_bei <- bei %>%
  create_grid_hexagonal(cell_size = cell_size, side_offset = band_width)
#Calculate KDE
kde <- bei %>%
  kde(band_width = band_width, kernel = "quartic", grid = grid_bei)
#Plot
tm_shape(kde) +
  tm_polygons(col = "kde_value",style="cont", palette = "viridis", title = "KDE Estimate",legend.show=FALSE)
```
![SpatialKDEViridis](https://github.com/ChenLaboratory/Hoang/assets/99466326/380d7e48-9529-4fbf-81a3-067b4415d695)


### MERFISH dataset
Using spatial transcriptomic dataset. This example use MERFISH mouse hypothalamic preoptic region dataset from [MerfishData package](https://bioconductor.org/packages/release/data/experiment/html/MerfishData.html) on bioconductor. Finding density for cells of 'inhibitory' cell class.

```
library(MerfishData)
spe = MouseHypothalamusMoffitt2018()
```

```
#hexDensity
plotHexDensity(hexDensity(spe,assay='exprs',sigma=20,
               weight='cell_class', #look for cell_class in either rownames or colData
               weightTransform='Inhibitory') #get only "Inhibitory")
```
![Rplot05](https://github.com/ChenLaboratory/hexDensity/assets/99466326/4952443f-e4ea-489d-b4ca-1706aeeae540)

```
#comparison with density.ppp by spatstat
cdat = data.frame(colData(spe),spatialCoords(spe))
cdat.inhibitory = subset(cdat, cell_class == ("Inhibitory"))
#need to convert into ppp object
cdat.inhibitory = ppp(cdat.inhibitory$x,cdat.inhibitory$y,window = owin(range(cdat.inhibitory$x),range(cdat.inhibitory$y)))

density.inhibitory = density.ppp(cdat.inhibitory,sigma=20, eps=diff(range(cdat.inhibitory$x))/128)
plot.im(density.inhibitory, col=colorRampPalette(viridis::viridis(11)))
```
![Rplot06](https://github.com/ChenLaboratory/hexDensity/assets/99466326/d9d04d84-8a0f-41cc-a649-a3bb61e9d601)
