#' Plot method for hexDensity
#'
#' Adapted plotting function for hexbin object. 
#' 
#' @param hexDensity hexbin object returned by hexDensity
#' @param colramp Color function that accept an integer n and return n colors.
#' @param main Main title
#' @param legend Legend is currently non-functional and should be ignored.
#' @param colorcut An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1] (e.g. c(0,0.5,0.75,1) would assign different colors, as determined by colramp, for the bottom 50% values, the next 25%, and the top 25% values).
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom grid grid.newpage viewport pushViewport upViewport grid.xaxis grid.yaxis grid.text grid.rect gpar unit
#' @importFrom hexbin hexcoords hexpolygon hcell2xy

plotHexDensity = function(hexDensity, 
                      main=deparse(substitute(hexDensity)), xlab=NULL, ylab=NULL,
                      xaxt=TRUE, yaxt=TRUE,
                      lcex=1,
                      colramp = colorRampPalette(viridis::viridis(11)), colorcut=1024,
                      legend=T, legendWidth=0.05, legendDistance=0.15,
                      aspectRatio=1/hexDensity@shape,
                      margin=0.2,
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
    grid.text(xlab, y = unit(-2, "lines"), gp = gpar(fontsize = 16,cex=lcex))
  if(nchar(ylab) > 0)
    grid.text(ylab, x = unit(-2, "lines"), gp = gpar(fontsize = 16,cex=lcex), rot = 90)
  if(nchar(main) > 0)
    grid.text(main, y = unit(1, "npc") + unit(1.5, "lines"),
              gp = gpar(fontsize = 18,cex=lcex))
  
  #hexagons
  upViewport()
  hexViewport$clip = T
  pushViewport(hexViewport)
  grid.hexagons(hexDensity,
                check.erosion = FALSE,
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

grid.hexagons <-
  function(dat,
           use.count=TRUE, cell.at=NULL,
           check.erosion = TRUE,
           trans = NULL,
           colorcut = seq(0, 1, length = 17),
           colramp = function(n){ LinGray(n,beg = 90, end = 15) },
           def.unit = "native"
  {
    ## Warning:	 presumes the plot has the right shape and scales
    ##		 See plot.hexbin()
    ## Arguments:
    ##	dat   = hexbin object
    ##	style = type of plotting
    ##		'centroids' =  symbol area is a function of the count,
    ##			approximate location near cell center of
    ##			mass without overplotting
    ##		'lattice' = symbol area is a function of the count,
    ##			plot at lattice points
    ##		'colorscale'   = gray scale plot,
    ##			color number determined by
    ##			transformation and colorcut,
    ##			area = full hexagons.
    ##		'nested.lattice'= plots two hexagons
    ##		   background hexagon
    ##			area=full size
    ##			color number by count in powers of 10 starting at pen 2
    ##		   foreground hexagon
    ##			area by log10(cnt)-floor(log10(cnt))
    ##			color number by count in powers of 10 starting at pen 12
    ##		'nested.centroids' = like nested.lattice
    ##			but counts < 10 are plotted
    ##
    ##	minarea =  minimum symbol area as fraction of the binning cell
    ##	maxarea =  maximum symbol area as fraction of the binning cell
    ##	mincnt	=  minimum count accepted in plot
    ##	maxcnt	=  maximum count accepted in plot
    ##	trans	=  a transformation scaling counts into [0,1] to be applied
    ##		   to the counts for options 'centroids','lattice','colorscale':
    ##			 default=(cnt-mincnt)/(maxcnt-mincnt)
    ##	colorcut=  breaks for translating values between 0 and 1 into
    ##		   color classes.  Default= seq(0,1,17),
    ##	density =  for hexagon graph paper
    ##	border	   plot the border of the hexagon, use TRUE for
    ##		   hexagon graph paper
    ## Symbol size encoding:
    ##	 Area= minarea + scaled.count*(maxarea-minarea)
    ##	 When maxarea==1 and scaled.count==1, the hexagon cell
    ##	 is completely filled.
    ##
    ##	 If small hexagons are hard to see increase minarea.
    ## For gray scale encoding
    ##	  Uses the counts scaled into [0,1]
    ##	  Default gray cutpoints seq(0,1,17) yields 16 color classes
    ##	  The color number for the first class starts at 2.
    ##	       motif coding: black 15 white puts the first of the
    ##			     color class above the background black
    ##	  The function subtracts 1.e-6 from the lower cutpoint to include
    ##	  the boundary
    ## For nested scaling see the code
    ## Count scaling alternatives
    ##
    ##	log 10 and Poisson transformations
    ##	trans <- function(cnt) log10(cnt)
    ##	   min	inv   <- function(y) 10^y
    ##
    ##	     trans <- function(cnt) sqrt(4*cnt+2)
    ##	     inv   <- function(y) (y^2-2)/4
    ## Perceptual considerations.
    ##	  Visual response to relative symbol area is not linear and varies from
    ##	  person to person.  A fractional power transformation
    ##	  to make the interpretation nearly linear for more people
    ##	  might be considered.	With areas bounded between minarea
    ##	  and 1 the situation is complicated.
    ##
    ##	  The local background influences color interpretation.
    ##	  Having defined color breaks to focus attention on
    ##	  specific countours can help.
    ##
    ## Plotting the symbols near the center of mass is not only more accurate,
    ##	  it helps to reduce the visual dominance of the lattice structure.  Of
    ##	  course higher resolution binning reduces the possible distance between
    ##	  the center of mass for a bin and the bin center.  When symbols
    ##	  nearly fill their bin, the plot appears to vibrate.  This can be
    ##	  partially controlled by reducing maxarea or by reducing
    ##	  contrast.
    
    
    ##____________________Initial checks_______________________
    if(!is(dat,"hexbin"))
      stop("first argument must be a hexbin object")
    ##_______________ Collect computing constants______________
    
    if(use.count){
      cnt <- dat@count
    }
    else{
      cnt <- cell.at
      if(is.null(cnt)){
        if(is.null(dat@cAtt)) stop("Cell attribute cAtt is null")
        else cnt <- dat@cAtt
      }
    }

    
    ##___________Transform Counts to range [0,1]_____________________
   if(!is.null(trans)) {
     cnt = trans(cnt) 
     if(any(is.na(rcnt)))
        stop("bad count transformation")
   }
   range = range(cnt)
   rcnt <- {
       if(range[1] == range[2]) rep.int(1, length(cnt))
       else (cnt - range[1])/(range[2]-range[1])
   }

    ##______________Set Colors_____________________________
             ## MM: Following is quite different from bin2d's
   nc <- length(colorcut)
   if(colorcut[1] > colorcut[nc]){
     colorcut[1] <- colorcut[1] + 1e-06
     colorcut[nc] <- colorcut[nc] - 1e-06
   } else {
     colorcut[1] <- colorcut[1] - 1e-06
     colorcut[nc] <- colorcut[nc] + 1e-06
   }
   colgrp <- cut(rcnt, colorcut,labels = FALSE)
   if(any(is.na(colgrp))) colgrp <- ifelse(is.na(colgrp),0,colgrp)
   ##NL: colramp must be a function accepting an integer n
   ##    and returning n colors
   clrs <- colramp(length(colorcut) - 1)
   pen <- clrs[colgrp]
   
   # Speed up plotting a bit by setting the most frequent color as background  
   # so don't have to plot those hexagons. 
   # Only worth if safe ~>1000 hexagons when tested.
   mostFreqPen = which.max(table(pen))
   if (mostFreqPen > 1000) { 
     mostFreqPen = names(mostFreqPen)
     grid.rect(gp=gpar(col=F,fill=mostFreqPen))
     notMostFreq=(pen!=mostFreqPen)
     pen = pen[notMostFreq]
     dat@cell=dat@cell[notMostFreq] #safe since R is pass-by-value
     # if (verbose) cat(length(pen),"hexagons drawn out of",length(dat@count),"hexagons")
   }
    ##__________________ Construct a hexagon___________________
   
   xbins <- dat@xbins
   shape <- dat@shape
   tmp <- hcell2xy(dat, check.erosion = check.erosion)
   xnew <- tmp$x
   ynew <- tmp$y
   sx <- xbins/diff(dat@xbnds)
   sy <- (xbins * shape)/diff(dat@ybnds)
   
    ## The inner and outer radius for hexagon in the scaled plot
    inner <- 0.5
    outer <- (2 * inner)/sqrt(3)
    ## Now construct a point up hexagon symbol in data units
    dx <- inner/sx
    dy <- outer/(2 * sy)
    # rad <- sqrt(dx^2 + dy^2)
    hexC <- hexcoords(dx, dy, sep=NULL)
    ##_______________ Full Cell	 Plotting_____________________
  

   hexpolygon(xnew, ynew, hexC,
              fill = pen)
   ## and that's been all for these styles
   return(invisible(paste("done")))
  }