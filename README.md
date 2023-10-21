# Kernel Density with Hexagon
## Functions
### hexKernel: calculating kernel density with regular hexagonal grid
Edge correction has been implemented, including Jones-Diggle algorithm as described in Jones, M.C. (1993) Simple boundary corrections for kernel density estimation. Statistics and Computing 3, 135--146.
```
hexKernel(
  x,
  y=NULL,
  xbins=128,
  sigma=1,
  edge = TRUE,
  diggle = FALSE
)
```
#### Arguments:
x, y: Coords of the points or a single plotting structure to be passed into hexbinFullRegular. See [xy.coords](https://www.rdocumentation.org/packages/grDevices/versions/3.6.2/topics/xy.coords)

xbins: number of bins in row

sigma: bandwidth for kernel density calculation

edge: logical value for whether to apply edge correction

diggle: logical value for apply edge correction with the Jones-Diggle improved edge correction which is more accurate. (need 'edge' to be TRUE to take effect).

### plotKernel: plot result from hexDensity
Adapted from plotting function of hexbin. Hacky workaround since hexbin plot things into discrete bins whereas results from hexKernel is more like a gradient and need a continuous color spectrum. Will be changed in the future.
```
plotKernel(
  kernel,
  main = deparse(substitute(kernel)),
  legend = F,
  colramp = colorRampPalette(rev(brewer.pal(11,'Spectral'))),
  colcut = seq(0,1,length = 100*kernel@ncells))
```

#### Arguments:
kernel: resulting hexbin object from hexKernel

main: header of the plot

legend: legend. Currently non-functional

colramp: color pallete. function accepting an integer n and return n colors

colcut: ignore for now

### hexbinFullRegular: hexagonal binning with regular hexagon
Adapted from hexbin to use regular hexagon and also returns hexagon with no counts
```
hexbinFullRegular(
  x,
  y = NULL,
  xbins = 128,
  xbnds = range(x),
  ybnds = range(y),
  xlab = NULL,
  ylab = NULL,
  IDs = FALSE
)
```
#### Arguments
See [hexbin](https://www.rdocumentation.org/packages/hexbin/versions/1.29.0/topics/hexbin) function from hexbin package. 'shape' argument is omitted to ensure hexagon is regular.
## Demonstration
### Data Preparation
Using bei spatial dataset from spatstat.explore
```
library(spatstat.data)
data=bei
```

Alternatively, use a real biological spatial transcriptomic dataset (larger & slower to load). This example use MERFISH mouse hypothalamic preoptic region dataset from [MerfishData package](https://bioconductor.org/packages/release/data/experiment/html/MerfishData.html) on bioconductor. Data is subsetted for cells of 'inhibitory' cell class only.
```
library(MerfishData)
spe = MouseHypothalamusMoffitt2018()
#subsetting data 
cdat = data.frame(colData(spe),spatialCoords(spe))
cdat = subset(cdat, cell_class!= "Ambiguous",select = -c(cell_id,sample_id,sex,behavior,neuron_cluster_id))
cdat = subset(cdat,z == -0.14)
cdat.inhibitory = subset(cdat, cell_class == ("Inhibitory"))
#convert to spatstat 'ppp' object for later comparison using density.ppp 
cdat.inhibitory = ppp(cdat.inhibitory$x,cdat.inhibitory$y,window = owin(range(cdat.inhibitory$x),range(cdat.inhibitory$y)))
```

### Kernel density
Calculating kernel density using hexagonal grid
```
density = hexKernel(data,sigma=48)
```

### Plot result
As mentioned in the Functions section, plotKernel adapted the plotting function from hexbin which is used to plot discrete data instead of the continuous range of KDE so the image's color may not accurately reflect the underlying value. 
```
plotKernel(density)
```
![hexKernelWEdgeCorrection](https://github.com/ChenLaboratory/Hoang/assets/99466326/1d76789f-045d-47c2-bbae-861e8fd3b271)

Comparing to density.ppp by spatstat which use square-grid (sigma and color may need to be tweaked for better visual match)
```
library(spatstat.explore)
#eps variable is used to turn the grid square instead of rectangle 
density = density.ppp(data,sigma=25, eps=diff(range(data$x))/128)
plot.im(density,col=colorRampPalette(viridis::viridis(11)))
```
![densitypppWEdgeCorrection](https://github.com/ChenLaboratory/Hoang/assets/99466326/88bb3de7-52f1-4e7a-8879-95cb49df4a40)

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

Using MERFISH dataset
```
plotKernel(hexKernel(cdat.inhibitory,sigma=20))
```

![MerfishHexWEdgeCorrection](https://github.com/ChenLaboratory/Hoang/assets/99466326/2d637d87-58a4-4ba6-a2ff-86a3121d0014)

and with density.ppp by spatstat
```
plot.im(density.ppp(cdat.inhibitory,sigma=20, eps=diff(range(data$x))/128),col=colorRampPalette(viridis::viridis(11)))
```
![MerfishDensityWEdgeCorrection](https://github.com/ChenLaboratory/Hoang/assets/99466326/1a0db07c-2339-438c-ba98-653480c84930)


