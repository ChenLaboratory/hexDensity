# Kernel Density with Hexagon
## Functions
### hexKernel: calculating kernel density with regular hexagonal grid
```
hexKernel(
  x,
  y=NULL,
  xbins=128,
  sigma=1
)
```
#### Arguments:
x, y: Coords of the points or a single plotting structure to be passed into hexbinFullRegular. See [xy.coords](https://www.rdocumentation.org/packages/grDevices/versions/3.6.2/topics/xy.coords)

xbins: number of bins in row

sigma: bandwidth for kernel density calculation

### plotKernel: plot result from hexDensity
Adapted from plotting function of hexbin. Hacky workaround since hexbin plot things into discrete bins whereas results from hexKernel is more like a gradient. Will be changed in the future.
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
density = hexKernel(data,sigma=20)
```

### Plot result
```
plotKernel(density)
```
![hexkernel](https://github.com/ChenLaboratory/Hoang/assets/99466326/5db9d6f4-64c1-4c29-bff7-6bbb3fe46832)

Comparing to density.ppp by spatstat which use square-grid (sigma and color may need to be tweaked for better visual match)
```
library(spatstat.explore)
density = density.ppp(cdat.inhibitory.ppp,sigma=20)
plot.im(density)
```
![spatstatkernel](https://github.com/ChenLaboratory/Hoang/assets/99466326/a8c77cb5-b566-4061-b144-e69cdbcdd8ed)

Using MERFISH dataset
```
plotKernel(hexKernel(cdat.inhibitory,sigma=20))
```
![b](https://github.com/ChenLaboratory/Hoang/assets/99466326/56a462c1-2ba7-46df-8157-0e7d4e615d4c)

and with density.ppp by spatstat
```
plot.im(density.ppp(cdat.inhibitory.ppp,sigma=20))
```
![a](https://github.com/ChenLaboratory/Hoang/assets/99466326/cf05ed52-7100-4a6d-940a-f9bf488a98c5)

