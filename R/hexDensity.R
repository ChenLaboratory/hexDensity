#' Kernel Density Estimation with Hexagonal grid.
#' 
#' @param x,y Coords of the points or a single plotting structure to be used in binning. See \link[grDevices]{xy.coords}.
#' @param xbins Number of bins in a row.
#' @param bandwidth Bandwidth for the smoothing kernel. Either a scalar, a vector of length 2, or a 2x2 variance-covariance matrix for the bandwidths in the x and y directions.
#' @param edge Logical value for whether to apply edge correction. Default is TRUE.
#' @param diggle Logical value for apply edge correction with the more accurate Jones-Diggle method (need 'edge' to be TRUE).
#' @param weight numeric weight vector to be assigned to points.
#' @param ... arguments for \link[hexDensity]{hexbinFull}
#' @return an S4 object of class \link[hexbin]{hexbin}.
#' @importFrom spatstat.geom fft2D
#' @importFrom grDevices xy.coords
#' @importFrom stats dnorm
#' 
#' @details Default bandwidth is 1/8 of the range of the smaller dimensions.
#' 
#' @references Diggle, P. J. (2010) Nonparametric methods. Chapter 18, pp. 
#' 299--316 in A.E. Gelfand, P.J. Diggle, M. Fuentes and P. Guttorp (eds.) 
#' Handbook of Spatial Statistics, CRC Press, Boca Raton, FL.
#' @references Jones, M. C. (1993) Simple boundary corrections for kernel 
#' density estimation. Statistics and Computing 3, 135--146.
#' 
#' @export
#' @examples
#' set.seed(133)
#' d = hexDensity(x=rnorm(200),y=rnorm(200),bandwidth=0.15)
hexDensity = function(x,y=NULL,
                      xbins = 128, #128 is the magic number in spatstat
                      bandwidth = NULL,
                      edge = TRUE,
                      diggle = FALSE,
                      weight = NULL,...) {
  hbin = hexbinFull(x,y,xbins=xbins, weight=weight,...) 
  row = hbin@dimen[1]
  col = hbin@dimen[2]
  
  xy <- xy.coords(x, y)
  
  # bandwidth
  if (is.null(bandwidth)) {
    bandwidth = min(diff(range(xy$x)),diff(range(xy$y)))/8
  }
  if (length(bandwidth)==1) {
    bandwidth=c(bandwidth,bandwidth)
  }
  
  xhex = diff(hbin@xbnds)/xbins
  yhex = xhex*diff(hbin@ybnds)/(diff(hbin@xbnds)*hbin@shape)
  
  # Calculate KDE for hex within the boundaries only since hexbin use more 
  # hexagons than needed. Important for edge correction.
  # topRow = floor(((hbin@ybnds[2]-hbin@ybnds[1]+yhex/sqrt(3))*2/sqrt(3))/yhex+1)
  topRow = floor((2*sqrt(3)*diff(hbin@ybnds)/yhex+5)/3)
  
  staggeredBin = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  for (i in seq(row-topRow+1,row)) {
    staggeredBin[i,(i-i%/%2):(i-i%/%2+col-1)] = hbin@count[((row-i)*col+1):((row-i)*col+col)]
  }
  
  #Make kernel
  xcol.ker.left = xhex*seq(col,1-col)
  yrow.ker.left = yhex*sqrt(3)*seq(row/2-1,-row/2)
    
  xcol.ker.right = xhex*(seq(col,1-col)-0.5)
  yrow.ker.right = yhex*sqrt(3)*(seq(row/2-1,-row/2)+0.5)
  if (length(bandwidth)==2) {
    kernel.left.hori = dnorm(xcol.ker.left,sd=bandwidth[1])
    kernel.left.verti = dnorm(yrow.ker.left,sd=bandwidth[2])
    kernel.left=outer(kernel.left.verti,kernel.left.hori)
        
    kernel.right.hori = dnorm(xcol.ker.right,sd=bandwidth[1])
    kernel.right.verti = dnorm(yrow.ker.right,sd=bandwidth[2])
    kernel.right=outer(kernel.right.verti,kernel.right.hori)
        
  } else if (is.matrix(bandwidth) && all.equal(dim(bandwidth),c(2,2))) {
    # anisotropic kernel
    Sinv = solve(bandwidth)
    halfSinv = Sinv/2
        
    Ypad = matrix(0, ncol=2*col, nrow=row)
    xsq.left <- matrix((xcol.ker.left^2)[col(Ypad)], ncol=2*col, nrow=row)
    ysq.left <- matrix((yrow.ker.left^2)[row(Ypad)], ncol=2*col, nrow=row)
    xy.left <- outer(yrow.ker.left, xcol.ker.left, "*")
    kernel.left <- exp(-(xsq.left * halfSinv[1,1]
                         - xy.left * (halfSinv[1,2]+halfSinv[2,1])
                         + ysq.left * halfSinv[2,2]))
    
    xsq.right <- matrix((xcol.ker.right^2)[col(Ypad)], ncol=2*col, nrow=row)
    ysq.right <- matrix((yrow.ker.right^2)[row(Ypad)], ncol=2*col, nrow=row)
    xy.right <- outer(yrow.ker.right, xcol.ker.right, "*")
    kernel.right <- exp(-(xsq.right * halfSinv[1,1]
                          - xy.right * (halfSinv[1,2]+halfSinv[2,1])
                          + ysq.right * halfSinv[2,2]))
  }
    
  #staggered bin
  kernel = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  for (i in seq(1, row)) {
    kernel[i*2-1,(i:(i+2*col-1))] = kernel.right[i,]
    kernel[i*2,(i:(i+2*col-1))] = kernel.left[i,]
  }
  kernel = kernel/sum(kernel)
  
  #inverse the kernel so that the center is now in top left for convolution.
  kernel.inv = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  #going clock-wise from the top-left section of inverse kernel
  kernel.inv[1:(row+1),1:(col+row/2)] = kernel[row:(2*row),(col+row/2):(2*col+row-1)]
  kernel.inv[1:(row+1),(col+row/2+1):(2*col+row-1)] = kernel[row:(2*row),1:(col+row/2-1)]
  kernel.inv[(row+2):(2*row),(col+row/2+1):(2*col+row-1)] = kernel[1:(row-1),1:(col+row/2-1)]
  kernel.inv[(row+2):(2*row),1:(col+row/2)] = kernel[1:(row-1),(col+row/2):(2*col+row-1)]
  
  fK = fft2D(kernel.inv)
  
  #edge correction
  if (edge) {
    mask = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
    for (i in seq(row-topRow+1, row, by=2)) {
      mask[i:(i+1),(i-i%/%2):(i-i%/%2+col-1)] = matrix(1,2,col)
    }
    fM = fft2D(mask)
    con = fft2D(fM * fK, inverse=TRUE)
    con = con/(2*row*(2*col+(row-1)))
    edg <- Mod(con[1:row, 1:(col+(row-1)/2)])
    
    if(diggle) {
      staggeredBin[1:row,1:(col+(row-1)/2)] = staggeredBin[1:row,1:(col+(row-1)/2)]/edg
      # Remove NaN (from 0/0 in the boundary of the staggered bin) for fft2D.
      staggeredBin[is.nan(staggeredBin)] = 0
    }
  }
  
  #KDE calculation by convolution
  fY = fft2D(staggeredBin)
  sm = fft2D(fY*fK,inverse = TRUE)/(2*row*(2*col+(row-1)))
  
  #edge correction
  if(edge && !diggle){
    #No need to clean up the NaN since will discard them anyway
    sm[1:row,1:(col+(row-1)/2)] = Re(sm[1:row,1:(col+(row-1)/2)])/edg
  }
  
  #extract back to hexbin class
  count=vector(mode="numeric", length = topRow*col)
  for (i in seq(row,row-topRow+1)) {
    count[(1+col*(row-i)):(col*(row-i)+col)] = Re(sm[i,(i-i%/%2):(i-i%/%2+col-1)])
  }
  
  hbin@count = count/(hexAreaFromWidth(xhex)*yhex/xhex)
  hbin@ncells=length(hbin@count)
  hbin@cell=hbin@cell[1:hbin@ncells]
  hbin@xcm=hbin@xcm[1:hbin@ncells]
  hbin@ycm=hbin@ycm[1:hbin@ncells]
  return(hbin)
}