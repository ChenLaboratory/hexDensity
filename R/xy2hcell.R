#' Get the hexagonal cell ID from a hexbin object at the specified xy coordinates.
#' 
#' @param hexbin hexbin object to be referenced to.
#' @param x,y coordinates or vectors of coordinates of the points. 
#'
#' @return a vector the same length as x with the hexagonal cell ID for each point
#' @export
#'
#' @examples
#' set.seed(133)
#' d=hexbin(x=rnorm(20000),y=rnorm(20000),xbins=50)
#' xy2hcell(d,x=0.5,y=0.2)
xy2hcell <- function(hexbin, x, y=NULL) 
{
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
  
  jmax <- floor(hexbin@xbins + 1.5001)
  #default shape to make regular hexagon
  # if (is.null(shape)) shape = diff(ybnds)/diff(xbnds)
  c1 <- 2 * floor((hexbin@xbins*hexbin@shape)/sqrt(3) + 1.5001)
  imax <- trunc((jmax*c1 -1)/jmax + 1)
  lmax <- jmax * imax
  weight = rep(1,length=n)
  
  #get cell IDs for all x,y at the same hexbin's specs
  ans <- .Fortran(`hbin`,
                  x = as.double(x),
                  y = as.double(y),
                  cell = integer(lmax),
                  cnt = double(lmax),
                  xcm = double(lmax),
                  ycm = double(lmax),
                  xbins = as.double(hexbin@xbins),
                  shape = as.double(hexbin@shape),
                  xbnds = as.double(hexbin@xbnds),
                  ybnds = as.double(hexbin@ybnds),
                  dim = as.integer(c(imax, jmax)),
                  n = as.integer(n),
                  cID = integer(n),
                  weight = as.double(weight))$cID
  # print(ans)
  if(has.na) {
    ok <- as.integer(ok)
    ok[!na] <- ans
    ok[na] <- NA
    ans <- ok
  }
  return(ans)
}