#' Plot method for hexDensity
#'
#' Adapted plotting function for hexbin object. 
#' 
#' @param hexDensity hexbin object returned by hexDensity
#' @param main Main title
#' @param xlab,ylab x-axis and y-axis label
#' @param xaxt,yaxt logical for whether to plot x-axis and y-axis values
#' @param lcex Expansion factor for all letters.
#' @param colramp Color function that accept an integer n and return n colors.
#' @param colorcut An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1].
#' @param legend Legend is currently non-functional and should be ignored.
#' @param legendWidth Expansion factor for legend width.  
#' @param legendDistance Expansion factor for the space between the plot and the legend.s
#' @param aspectRatio width to height ratio of the plot. Default is the (inverse of) shape value of hexDensity.
#' @param margin Minimum guaranteed margin for the plot. Different aspect ratio between the screen and the plot means that margin can be larger on certain sides.
#' @param newpage logical for whether to plot on a new page.
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(133)
#' d = hexDensity(x=rnorm(200),y=rnorm(200),bandwidth=0.15)
#' plotHexDensity(d)
#' @importFrom grid grid.newpage viewport pushViewport upViewport grid.xaxis grid.yaxis grid.text grid.rect gpar unit grid.pretty
#' @importFrom viridis viridis
plotHexDensity = function(hexDensity, 
                      main=NULL, xlab=NULL, ylab=NULL,
                      xaxt=TRUE, yaxt=TRUE,
                      lcex=1,
                      colramp = colorRampPalette(viridis(11)), colorcut=1024,
                      legend=T, legendWidth=0.05, legendDistance=0.15,
                      aspectRatio=1/hexDensity@shape,
                      margin=0.18,
                      newpage=T) {
  if(!is(hexDensity,"hexbin"))
    stop("first argument must be a hexbin object")
  if (length(colorcut) > 1) { # a sequence 0,...,1
    if(colorcut[1] != 0)
      stop("Colorcut lower boundary must be 0")
    if(colorcut[length(colorcut)] != 1)
      stop("Colorcut upper boundary must be 1")
  }
  else {
    colorcut <-
      if(colorcut > 1) seq(0, 1, length = colorcut)
    else 1
  }
  
  ## -----Prepare viewports ----------------------
  screenRatio = dev.size()[1]/dev.size()[2]
  plotsize = 1-margin*2
  if (aspectRatio > screenRatio) {
    w = unit(plotsize,'npc')
    h = unit(plotsize*screenRatio/aspectRatio,'npc')
  }
  else {
    w = unit(plotsize*aspectRatio/screenRatio,'npc')
    h = unit(plotsize,'npc')
  }
  
  if (legend) {
    legendWidth = unit(legendWidth*w,'npc')
    legendDistance = unit(legendDistance*w,'npc')
    #legend size is based on plot height, which is bigger, instead of width
    if (aspectRatio<1) {
      legendWidth = legendWidth/aspectRatio
      legendDistance = legendDistance/aspectRatio
    }
    
    #need rescale to maintain margin
    if (c(w + legendWidth + legendDistance) > plotsize) {
      rescaleFactor = plotsize/c(w + legendWidth + legendDistance)
      w = w*rescaleFactor
      h = h*rescaleFactor
      legendWidth = legendWidth*rescaleFactor
      legendDistance = legendDistance*rescaleFactor
    }
    legendViewport = viewport(x=unit(0.5,'npc') + unit(w/2 + legendDistance/2,'npc'),
                              width = legendWidth,
                              height = h,
                              yscale=range(hexDensity@count),
                              clip=F
      )
  }
  hexViewport = viewport(x=unit(0.5,'npc')-unit(as.numeric(legend)*(legendWidth+legendDistance)/2,'npc'),
                         width=w,
                         height=h,
                         xscale=hexDensity@xbnds, yscale=hexDensity@ybnds,
                         default.units = 'native',
                         clip=F
                         )
  
  ## ----- plotting starts ------------------------
  if (newpage) grid.newpage()
  pushViewport(hexViewport)
  #labels
  if(xaxt) grid.xaxis(gp=gpar(cex=lcex))
  if(yaxt) grid.yaxis(gp=gpar(cex=lcex))
  ## xlab, ylab, main :
  if(is.null(xlab)) xlab <- hexDensity@xlab
  if(is.null(ylab)) ylab <- hexDensity@ylab
  if(nchar(xlab) > 0)
    grid.text(xlab, 
              y = unit(-1.5 - xaxt*1, "lines"),
              gp = gpar(fontsize = 16,cex=lcex))
  if(nchar(ylab) > 0)
    #Ensured no overlap with the yaxis ticks with 'strwidth'.
    #Not perfect since ticks can be like c(0,0.05,0.1) where the last value is 
    #not the longest. -1.5 lines should be enough wiggle-room.
    yTicks = grid.pretty(hexViewport$yscale)
    grid.text(ylab, x = unit(-1.5, "lines") - 
                unit(yaxt*1,'strwidth',as.character(yTicks[length(yTicks)]))
              , gp = gpar(fontsize = 16,cex=lcex), rot = 90)
  if(!is.null(main) && nchar(main) > 0)
    grid.text(main, y = unit(1,'npc') + unit(2, "lines"),
              gp = gpar(fontsize = 18,cex=lcex))
  
  #hexagons
  upViewport()
  hexViewport$clip = T
  pushViewport(hexViewport)
  grid.hexagontile(hexDensity,
                colorcut = colorcut,
                colramp = colramp)
  grid.rect(gp=gpar(fill=NA))

  upViewport()
  # ----- Legend ------------------------
  if (legend) {
    pushViewport(legendViewport)
    #Make ribbon legend
    grid.rect(y = unit(seq(0,1-1/256,length=256), "npc"),
              height = unit(1/256, "npc"), 
              just = "bottom",
              gp = gpar(col = NA, fill = colramp(256)))
    grid.yaxis(gp=gpar(cex=lcex),main=F)
    grid.rect(gp=gpar(col="black",fill=NA))
    upViewport()
  }
}
