library(hexbin)
dyn.load("hbin.so")
#modified hexbin that give full grid (include hex with no count) of regular hex
hexbinFullRegular <-
    function(x, y = NULL, xbins = 30,
	     xbnds = range(x), ybnds = range(y),
	     xlab = NULL, ylab = NULL, IDs = FALSE)
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
    n <- length(x)
    na <- is.na(x) | is.na(y)
    has.na <- any(na)
    if (has.na) {
	ok <- !na
	x <- x[ok]
	y <- y[ok]
        n0 <- n
        na.pos <- which(na)
	n <- length(x)
    }
    if(diff(xbnds) <= 0)
	stop("xbnds[1] < xbnds[2] is not fulfilled")
    if(!missing(xbnds) && any(sign(xbnds - range(x)) == c(1,-1)))
	stop("'xbnds' must encompass range(x)")
    if(diff(ybnds) <= 0)
	stop("ybnds[1] < ybnds[2] is not fulfilled")
    if(!missing(ybnds) && any(sign(ybnds - range(y)) == c(1,-1)))
	stop("'ybnds' must encompass range(y)")
    jmax <- floor(xbins + 1.5001)
    #imax <- 2 * floor((xbins * shape)/sqrt(3) + 1.5001)
    shape = diff(range(y))/diff(range(x))
    c1 <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
    imax <- trunc((jmax*c1 -1)/jmax + 1)
    lmax <- jmax * imax
    ans <- .Fortran("hbin",
	      x = as.double(x),
	      y = as.double(y),
	      cell = integer(lmax),
	      cnt = integer(lmax),
	      xcm = double(lmax),
	      ycm = double(lmax),
	      xbins = as.double(xbins),
	      shape = as.double(shape),
	      xbnds = as.double(xbnds),
	      ybnds = as.double(ybnds),
	      dim = as.integer(c(imax, jmax)),
	      n = as.integer(n),
	      cID = if(IDs) integer(n) else as.integer(-1))[-(1:2)]
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
    if(sum(ans$cnt) != n) warning("Lost counts in binning")
    new("hexbin",
    cell = ans$cell, count = ans$cnt,
    xcm = ans$xcm, ycm = ans$ycm, xbins = ans$xbins,
    shape = ans$shape, xbnds = ans$xbnds , ybnds = ans$ybnds,
    dimen = c(imax, jmax), n = n, ncells = ans$n,
    call = call, xlab = xlab, ylab = ylab, cID = ans$cID, cAtt = integer(0))
    #dimen = ans$dim
    }