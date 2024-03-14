#' Draw hexagon tiles with grid package
#'
#' Adapted of grid.hexagons from hexbin with some speedup specific for hexagonal tiling.
#'
#' @param hexDensity \link[hexbin]{hexbin} object returned by hexDensity.
#' @param use.count logical specifying if counts from hexbin object should be used.
#' @param cell.at numeric vector to be plotted instead of counts, must be same length as the number of cells.
#' @param trans a transformation function (or NULL) for the counts, e.g., sqrt.
#' @param colorcut An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1].
#' @param colramp Color function that accept an integer n and return n colors.
#' @param def.unit Default \link[grid]{unit} to be used.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom viridis viridis 
#' @importFrom hexbin hexcoords hexpolygon hcell2xy
grid.hexagontile <-
  function(hexDensity,
           use.count=TRUE, cell.at=NULL,
           trans = NULL,
           colorcut = seq(0, 1, length = 1024),
           colramp = colorRampPalette(viridis(11)),
           def.unit = "native")
  {
    ##____________________Initial checks_______________________
    if(!is(hexDensity,"hexbin"))
      stop("first argument must be a hexbin object")
    ##_______________ Collect computing constants______________
    
    if(use.count){
      cnt <- hexDensity@count
    }
    else{
      cnt <- cell.at
      if(is.null(cnt)){
        if(is.null(hexDensity@cAtt)) stop("Cell attribute cAtt is null")
        else cnt <- hexDensity@cAtt
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
    # Only worth if can save ~>1000 hexagons when tested.
    mostFreqPen = which.max(table(pen))
    if (mostFreqPen > 1000) { 
      mostFreqPen = names(mostFreqPen)
      grid.rect(gp=gpar(col=F,fill=mostFreqPen))
      notMostFreq=(pen!=mostFreqPen)
      pen = pen[notMostFreq]
      hexDensity@cell=hexDensity@cell[notMostFreq] #safe since R is pass-by-value
      # if (verbose) cat(length(pen),"hexagons drawn out of",length(hexDensity@count),"hexagons")
    }
    ##__________________ Construct a hexagon___________________
    
    xbins <- hexDensity@xbins
    shape <- hexDensity@shape
    tmp <- hcell2xy(hexDensity, check.erosion = F)
    xnew <- tmp$x
    ynew <- tmp$y
    sx <- xbins/diff(hexDensity@xbnds)
    sy <- (xbins * shape)/diff(hexDensity@ybnds)
    
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