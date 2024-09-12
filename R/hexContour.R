#' Generate contour using a hexagonal grid.
#'
#' This is made to follow the same behaviour as the isoband package. A dedicated plotting function is in the work. In the mean time, see example of how to plot the output with ggplot2::geom_path
#' Algorithm is a modification of the meandering triangles as described here https://blog.bruce-hill.com/meandering-triangles
#' @param hexDensity hexDensity object to be contoured.
#' @param levels Numeric vector for which contour lines should be generated
#' @return list of x, y, and ID, for the lines contouring each levels. 
#' @importFrom hexbin hcell2xy
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

# Essentially treat the hexagonal grid as a rectangular grid for calculation then
# shift every other rows by half bin width afterwards to get back the hexagonal grid.
hexContour = function(hexDensity,levels) {
  # Prepping values
  coords = hcell2xy(hexDensity)
  x.coords.left = coords$x[1:hexDensity@dimen[2]]
  x.coords.right = coords$x[(hexDensity@dimen[2]+1):(2*hexDensity@dimen[2])]
  y.coords = coords$y[((1:length(coords$y))-1)%%hexDensity@dimen[2]==0]
  z=matrix(hexDensity@count,ncol=hexDensity@dimen[2],byrow=T)
  
  isolines=meanderingTriangles(x.coords.left,x.coords.right,y.coords,z,levels)
  return(isolines)
}

#Meandering triangles c++
meanderingTriangles = function(x.coords.left,x.coords.right,y.coords,z,levels) {
  res = .Call(`meanderingTrianglesC`,x.coords.left,x.coords.right,y.coords,z,levels)
  names(res) = as.character(levels)
  return(res)
}

#Keep for posterity
# meanderingTriangles2 = function(x.coords.left,x.coords.right,y.coords,z,levels) {
#   n = (length(x.coords.left)-1)*(length(y.coords)-1)*2
#   triangles = vector(mode='list', length=n)
# 
#   i=1
#   for (x in 1:(length(x.coords.left)-1)) {
#     for (y in 1:(length(y.coords)-1)) {
#       if(y%%2==1) {
#         triangles[[i]]=list(v1=c(x,y),v2=c(x+1,y),v3=c(x,y+1))
#         triangles[[i+1]]=list(v1=c(x+1,y),v2=c(x,y+1),v3=c(x+1,y+1))
#       } else {
#         triangles[[i]]=list(v1=c(x,y),v2=c(x+1,y+1),v3=c(x,y+1))
#         triangles[[i+1]]=list(v1=c(x,y),v2=c(x+1,y+1),v3=c(x+1,y))
#       }
#       i=i+2
#     }
#   }
# 
# 
#   ## Find triangles that intersect contour lines
#   res = list()
#   for(level in levels) {
#     interpolatedPos = dict()
#     contour_segments = list()
#     for (triangle in triangles) {
# 
#       below = Filter(function(v) z[v[2],v[1]]<level,triangle)
#       above = Filter(function(v) z[v[2],v[1]]>=level,triangle)
# 
#       # No contour line for all above or below
#       if (length(below)==0 || length(above)==0){
#         next
#       }
# 
#       minority = `if`(length(above) < length(below),above,below)
#       majority = `if`(length(above) > length(below),above,below)
# 
#       contour_points = list()
#       crossed_edges = list(list(e1=minority[1],e2=majority[1]),
#                            list(e1=minority[1],e2=majority[2]))
# 
#       for (triangle_edge in crossed_edges) {
#         e1 = triangle_edge$e1[[1]]
#         e2 = triangle_edge$e2[[1]]
# 
#         crossing_point=c(
#           (e1[1]-e2[1])*0.5+e2[1],
#           (e1[2]-e2[2])*0.5+e2[2]
#         )
# 
#         #interpolation. Doing this roundabout way to avoid calculation error
#         how_far = ((level - z[e2[2],e2[1]]) /
#                      (z[e1[2],e1[1]] - z[e2[2],e2[1]]))
#         e1.x = `if`(e1[2]%%2==1,x.coords.left[e1[1]],x.coords.right[e1[1]])
#         e2.x = `if`(e2[2]%%2==1,x.coords.left[e2[1]],x.coords.right[e2[1]])
#         interpolated=c(
#           # (x.coords[e1[1]]-x.coords[e2[1]])*how_far+x.coords[e2[1]],
#           (e1.x-e2.x)*how_far+e2.x,
#           (y.coords[e1[2]]-y.coords[e2[2]])*how_far+y.coords[e2[2]]
#         )
#         interpolatedPos$set(crossing_point,interpolated)
# 
#         contour_points[[length(contour_points)+1]] = crossing_point
#       }
#       contour_segments[[length(contour_segments)+1]]=list(e1=contour_points[[1]],e2=contour_points[[2]])
# 
#     }
# 
#     #no contour. Can probably save time by checking this earlier
#     if (length(contour_segments)==0){
#       res[[as.character(level)]]=list(x=numeric(0),y=numeric(0),id=integer(0))
#       next
#     }
# 
#     ## Joining up
#     unused_segments = ordered_dict(rep(NA,length(contour_segments)),contour_segments)
#     segments_by_point = dict()
# 
#     for (segment in contour_segments) {
#       if (segments_by_point$has(segment$e1)){
#         temp = segments_by_point$get(segment$e1)
#         temp[length(temp)+1] = list(segment)
#         segments_by_point$set(segment$e1,temp)
#       } else {
#         segments_by_point$set(segment$e1,list(segment))
#       }
# 
#       if (segments_by_point$has(segment$e2)){
#         temp = segments_by_point$get(segment$e2)
#         temp[length(temp)+1] = list(segment)
#         segments_by_point$set(segment$e2,temp)
#       } else {
#         segments_by_point$set(segment$e2,list(segment))
#       }
#     }
#     n_unused_segments = length(contour_segments)
#     contour_lines=list()
#     while(n_unused_segments) {
#       line = deque(unused_segments$popitem()$key)
#       n_unused_segments = n_unused_segments-1
#       while (TRUE) {
# 
#         tail_candidates = segments_by_point$get(line$peek())
# 
#         tail_candidates = Filter(function(x) unused_segments$has(x),
#                                  tail_candidates)
# 
#         if (length(tail_candidates)>0) {
#           tail=tail_candidates[[1]]
#           line$push(`if`(all(tail$e2==line$peek()),tail$e1,tail$e2))
#           unused_segments$remove(tail)
#           n_unused_segments=n_unused_segments-1
#         }
#         head_candidates = segments_by_point$get(line$peekleft())
#         head_candidates = Filter(function(x) unused_segments$has(x),
#                                  head_candidates)
#         if(length(head_candidates)>0) {
#           head=head_candidates[[1]]
#           line$pushleft(`if`(all(head$e2==line$peekleft()),head$e1,head$e2))
#           unused_segments$remove(head)
#           n_unused_segments=n_unused_segments-1
#         }
#         if (length(tail_candidates)==0 && length(head_candidates)==0) {
#           contour_lines[length(contour_lines)+1] = list(line$as_list())
#           break
#         }
#       }
#     }
# 
#     isoline_format=list(x=c(),y=c(),id=c())
#     n=1
#     for(i in 1:length(contour_lines)) {
#       for (j in 1:length(contour_lines[[i]])) {
#         temp = interpolatedPos$get(contour_lines[[i]][[j]])
#         isoline_format$x[n]=temp[1]
#         isoline_format$y[n]=temp[2]
# 
#         isoline_format$id[n]=i
#         n=n+1
#       }
#     }
# 
#   res[[as.character(level)]]=isoline_format
#   }
#   return(res)
# }
