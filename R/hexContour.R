#' Generate contour for a hexagonal grid.
#'
#' Algorithm is a modification of the meandering triangles as described in 
#' https://blog.bruce-hill.com/meandering-triangles. See 
#' \link[isoband]{isolines} for details about the output.
#' 
#' @param hexDensity hexDensity object to be contoured.
#' @param levels Numeric vector for which contour lines should be generated
#' @return list of x, y, and ID, for the contour line at each levels. 
#' @importFrom hexbin hcell2xy
#' 
#' @details
#' This function is made to follow the same behaviour as
#' \link[isoband]{isolines}. A dedicated plotting function is in the work. 
#' Meanwhile, see example of how to plot the output with ggplot2::geom_path.
#' 
#' @export 
#' @examples
#' 
#' set.seed(133)
#' x=rnorm(200)
#' y=rnorm(200)
#' d = hexDensity(x=x,y=y,bandwidth=0.4)
#' cutoff=quantile(d@count,0.9)
#' lines = hexContour(d,cutoff)
#' 
#' library(ggplot2)
#' library(hexbin)
#' #plot against density
#' ggplot()+
#'   geom_point(
#'     aes(x=hcell2xy(d)$x,
#'         y=hcell2xy(d)$y,
#'         col=d@count)
#'  ) +
#'   scale_color_viridis_c()+
#'   geom_path(
#'     aes(
#'      x = lines[[1]]$x, y = lines[[1]]$y, group = lines[[1]]$id
#'     )
#'   )
#' 
#' #plot against data points
#' ggplot() +
#'   geom_point(
#'     aes(x=x,y=y)) +
#'   geom_path(
#'     aes(
#'       x = lines[[1]]$x, y = lines[[1]]$y, group = lines[[1]]$id
#'     )
#'   )
hexContour = function(hexDensity,levels) {
  # Prepping values
  coords = hcell2xy(hexDensity)
  x.coords.left = coords$x[1:hexDensity@dimen[2]]
  x.coords.right = coords$x[(hexDensity@dimen[2]+1):(2*hexDensity@dimen[2])]
  y.coords = coords$y[((1:length(coords$y))-1)%%hexDensity@dimen[2]==0]
  z=matrix(hexDensity@count,ncol=hexDensity@dimen[2],byrow=TRUE)
  
  return(meanderingTriangles(x.coords.left,x.coords.right,y.coords,z,levels))
}

#' Meandering triangles for hexagonal grid in C++
#' @param x.coords.left Vector for x coords of left-aligned rows (row 1,3,5,...)
#' @param x.coords.right Vector for x coords of right-aligned rows (row 2,4,6,...)
#' @param y.coords Vector for y coords of all rows.
#' @param z Matrix for elevation values for the grid point
#' @param levels Vector of z value cutoffs for contouring.
#' @details
#' This function is not meant to be used as is, unless you are very familiar 
#' with how hexContour works.
#' @export 
meanderingTriangles = function(x.coords.left,x.coords.right,y.coords,z,levels) {
  res = .Call(`meanderingTrianglesC`,x.coords.left,x.coords.right,y.coords,z,levels)
  names(res) = as.character(levels)
  return(res)
}