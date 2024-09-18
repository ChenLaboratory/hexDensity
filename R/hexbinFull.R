#' Hexagonal binning with whole grid output. 
#'
#' Adapted from hexbin to output the full grid (include hexagon with 0 count) and also to include weights. Default is to use regular hexagon. See \link[hexbin]{hexbin} for extra details.
#' @param x,y Coords of the points or a single plotting structure to be used in binning. See xy.coords.
#' @param xbins Number of bins in a row.
#' @param shape shape = yheight/xwidth of the plotting regions
#' @param xbnds,ybnds Horizontal and vertical limits of the binning region in x or y units respectively, must encompass range(x) or range(y) respectively; Vector of length 2
#' @param xlab,ylab Optional character strings used as labels for x and y. If NULL, sensible defaults are used.
#' @param IDs Logical indicating if the hexagonal cell ID for each point should be returned, see hexbin.
#' @param weight Numeric weight vector to be assigned to points.
#'
#' @return hexbin object
#' @export
#'
#' @examples
#' set.seed(133)
#' d=hexbinFull(x=rnorm(20000),y=rnorm(20000),xbins=50)
#' plotHexDensity(d)
#' 
#' @importClassesFrom hexbin hexbin
#' @importFrom grDevices xy.coords
#' @importFrom methods new
hexbinFull <-
    function(x, y = NULL, xbins = 128, shape = NULL,
	     xbnds = range(x), ybnds = range(y),
	     xlab = NULL, ylab = NULL, IDs = FALSE,
	     weight = NULL)
{
    call <- match.call()
    ## (x,y, xlab, ylab) dealing
    xl <- if (!missing(x)) deparse(substitute(x))
    yl <- if (!missing(y)) deparse(substitute(y))
    xy <- xy.coords(x, y, xl, yl)
    ch0 <- function(u) if(is.null(u)) "" else u
    xlab <- if (is.null(xlab)) ch0(xy$xlab) else xlab
    ylab <- if (is.null(ylab)) ch0(xy$ylab) else ylab
    if(! (is.character(xlab) || is.expression(xlab)))
        stop("xlab must be a character or expression")
    if(! (is.character(ylab) || is.expression(ylab)))
        stop("ylab must be a character or expression")
    
    x <- xy$x
    y <- xy$y
    na <- is.na(x) | is.na(y)
    has.na <- any(na)
    if (has.na) {
    	ok <- !na
    	x <- x[ok]
    	y <- y[ok]
    }
    n <- length(x)
    
    if(diff(xbnds) <= 0)
	stop("xbnds[1] < xbnds[2] is not fulfilled")
    if(any(sign(xbnds - range(x)) == c(1,-1)))
	stop("'xbnds' must encompass range(x)")
    if(diff(ybnds) <= 0)
	stop("ybnds[1] < ybnds[2] is not fulfilled")
    if(any(sign(ybnds - range(y)) == c(1,-1)))
	stop("'ybnds' must encompass range(y)")
    jmax <- floor(xbins + 1.5001)
    #default shape to make regular hexagon
    if (is.null(shape)) shape = diff(ybnds)/diff(xbnds)
    c1 <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
    imax <- trunc((jmax*c1 -1)/jmax + 1)
    lmax <- jmax * imax

    if(is.null(weight)) {
      weight = rep.int(1,times=n)
    }
    if(length(weight) != n) {
      stop("weight must be a vector with same length as x")
    }
    
    ans <- .Fortran(`hbin`,
	      x = as.double(x),
	      y = as.double(y),
	      cell = integer(lmax),
	      cnt = double(lmax),
	      xcm = double(lmax),
	      ycm = double(lmax),
	      xbins = as.double(xbins),
	      shape = as.double(shape),
	      xbnds = as.double(xbnds),
	      ybnds = as.double(ybnds),
	      dim = as.integer(c(imax, jmax)),
	      n = as.integer(n),
	      cID = if(IDs) integer(n) else as.integer(-1),
	      weight = as.double(weight))[-(1:2)]
    ## cut off extraneous stuff
    if(!IDs) ans$cID <- NULL
    if(IDs && has.na) {
      ok <- as.integer(ok)
      ok[!na] <- ans$cID
      ok[na] <- NA
      ans$cID <- ok
    }
    nc <- ans$n
    length(ans$cell) <- nc
    length(ans$cnt) <- nc
    length(ans$xcm) <- nc
    length(ans$ycm) <- nc
    if(sum(ans$cnt) != sum(weight)) warning("Lost counts in binning")
    new("hexbin",
    cell = ans$cell, count = ans$cnt,
    xcm = ans$xcm, ycm = ans$ycm, xbins = ans$xbins,
    shape = ans$shape, xbnds = ans$xbnds , ybnds = ans$ybnds,
    dimen = c(imax, jmax), n = n, ncells = ans$n,
    call = call, xlab = xlab, ylab = ylab, cID = ans$cID, cAtt = integer(0))
    }