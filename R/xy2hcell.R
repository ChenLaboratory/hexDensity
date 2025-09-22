#' Find the hexagon cells from xy coordinates given a hexbin object.
#' 
#' Find the hexagon cells IDs from xy coordinates given a hexbin object. Useful 
#' if you want to get the KDE value at a certain coordinate. 
#' @param hexbin hexbin object to be referenced to.
#' @param x,y coordinates or vectors of coordinates of the points. 
#' @param xbins number of bins partitioning the range of xbnds
#' @param xbnds,ybnds horizontal and vertical limit of the binning region. Must be numeric vector of length 2. 
#' @param shape shape = yheight/xwidth of the plotting regions.
#' @return a vector the same length as x with the hexagonal cell ID for each point
#' @details
#' If a hexbin object is not provided, parameters of the binning region (xbins, xbnds, ybnds, shape) can be used instead.
#' For finding the xy coordinates of the hexagons for a hexbin object, see \link[hexbin]{hcell2xy}.
#' @export
#' @importFrom grDevices xy.coords
#'
#' @examples
#' library(hexbin)
#' set.seed(133)
#' d=hexDensity(x=rnorm(20000),y=rnorm(20000),xbins=50)
#' #Get KDE value at the coordinate x=0,y=0
#' loc = xy2hcell(d,x=0,y=0)
#' d@count[loc]
xy2hcell <- function(hexbin = NULL, x, y=NULL, xbins = NULL, xbnds=NULL, ybnds=NULL, shape=NULL) 
{
  if(!is.null(hexbin)) {
    xbins=hexbin@xbins
    xbnds=hexbin@xbnds
    ybnds=hexbin@ybnds
    shape=hexbin@shape
  }
  
  xl <- if (!missing(x)) deparse(substitute(x))
  yl <- if (!missing(y)) deparse(substitute(y))
  xy <- xy.coords(x, y, xl, yl)
  
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
  
  jmax <- ceiling(xbins + 0.5)
  #default shape to make regular hexagon
  c1 <- 2 * floor((xbins*shape)/sqrt(3) + 1.5001)
  imax <- trunc((jmax*c1 -1)/jmax + 1)
  lmax <- jmax * imax
  weight = rep.int(1,times=n)
  
  #get cell IDs for all x,y at the same hexbin's paras
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
                  cID = integer(n),
                  weight = as.double(weight))$cID
  if(has.na) {
    ok <- as.integer(ok)
    ok[!na] <- ans
    ok[na] <- NA
    ans <- ok
  }
  return(ans)
}