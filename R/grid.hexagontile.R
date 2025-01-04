#' Draw hexagon tiles with grid package
#'
#' Adapted from \link[hexbin]{grid.hexagons} by hexbin with speedup specific 
#' for hexagonal tiling (avoid plotting the most abundance hexagons by setting 
#' its color as background).
#'
#' @param hexDensity \link[hexbin]{hexbin} object returned by hexDensity.
#' @param use.count logical specifying if counts from hexbin object should be used.
#' @param cell.at numeric vector to be plotted instead of counts, must be same length as the number of cells.
#' @param trans a transformation function (or NULL) for the counts, e.g., sqrt.
#' @param colorcut An integer for the number of equi-spaced colorcut in [0,1] to assign colors to values. Alternatively, a vector of custom colorcut spacing between [0, 1].
#' @param colramp Color function that accept an integer n and return n colors.
#' @param def.unit Default \link[grid]{unit} to be used.
#'
#' @export
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom grid grid.polygon
#' @section SIDE EFFECTS: Adds hexagons to plot 
#' 
#' @returns No return value
#' @author Dan Carr <dcarr@voxel.galaxy.gmu.edu>; ported and extended by 
#' Nicholas Lewin-Koh nikko@hailmail.net. Modified by Quoc Hoang Nguyen 
#' <nguyen.q@wehi.edu.au> for hexDensity.
grid.hexagontile <- function(hexDensity,
    use.count=TRUE, cell.at=NULL,
    trans = NULL,
    colorcut = seq(0, 1, length = 1024),
    colramp = colorRampPalette(col.viridis),
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
        grid.rect(gp=gpar(col=FALSE,fill=mostFreqPen))
        notMostFreq=(pen!=mostFreqPen)
        pen = pen[notMostFreq]
        hexDensity@cell=hexDensity@cell[notMostFreq] #safe since R is pass-by-value
    }
    ##__________________ Construct a hexagon___________________

    xbins <- hexDensity@xbins
    shape <- hexDensity@shape
    tmp <- hcell2xyInt2(hexDensity)
    col <- tmp$j
    row <- tmp$i
    n6=rep.int(6,hexDensity@ncells)
    sx <- xbins/diff(hexDensity@xbnds)
    sy <- (xbins * shape)/diff(hexDensity@ybnds)
    ## The inner and outer radius for hexagon in the scaled plot
    inner <- 0.5
    outer <- (2 * inner)/sqrt(3)
    ## Now construct a point up hexagon symbol in data units
    dx <- inner/sx
    dy <- outer/(2 * sy)
    
    ## This prevents any gaps between hexagons
    offset = c(hexDensity@xbnds[1] - dx,
               hexDensity@ybnds[1] - dy*2)
    xc = offset[1] + (0:(hexDensity@dimen[2]*2+1)) * dx
    yc = offset[2] + (0:(hexDensity@dimen[1]*3+1)) * dy
    
    xhex = xc[rep.int(2*col+!(row%%2),n6)+c(0,1,1,0,-1,-1)]
    yhex = yc[rep.int(3*row,n6)+c(2,1,-1,-2,-1,1)]
    
    grid.polygon(x = unit(xhex,"native"),
                 y = unit(yhex,"native"),
                 id.lengths = n6,
                 gp = gpar(col=0, fill= pen))

    return(invisible())
}

#Turn hbin@cell into row&col representation
hcell2xyInt2 <- function(hbin) {
    cell = hbin@cell-1
    return(
        list(
            i = cell %/% hbin@dimen[2] + 1,
            j = cell %% hbin@dimen[2] + 1
        )
    )
}