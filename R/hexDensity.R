#' Kernel Density Estimation with Hexagonal binning.
#'
#' @param x,y Coords of the points or a single plotting structure to be used in binning. See xy.coords.
#' @param xbins Number of bins in a row.
#' @param sigma Bandwidth for kernel density calculation.
#' @param edge Logical value for whether to apply edge correction. Default is TRUE.
#' @param diggle Logical value for apply edge correction with the more accurate Jones-Diggle methor (need 'edge' to be TRUE).
#' @param weight numeric weight vector to be assigned to points. For SpatialExperiment, this can be a name in rownames or colData (transformation may be needed to turn data into numerical values)
#' @return hexbin object.


#' @examples
#' @importFrom spatstat.geom fft2D
#' @importClassesFrom SpatialExperiment SpatialExperiment
#' @export
hexDensity = function(x,...) UseMethod('hexDensity')

#' @rdname hexDensity
#' @param assay assay to be used of the SpatialExperiment object 
#' @param weightTransform a function that can take in the weight vector to return a vector of same length. For SpatialExperiment, one may want to use this to transform string into numerical values (e.g.: transform cell_class in colData into binary boolean for selection of a specific class)
#' @export
hexDensity.SpatialExperiment = function(x,
                                        assay = assayNames(x)[1],
                                        xbins = 128, #128 is the magic number in pixellate of spatstat
                                        sigma = 1,
                                        edge = TRUE,
                                        diggle = FALSE,
                                        weight = NULL, # value will become weight
                                        weightTransform = NULL) { 
  xy = spatialCoords(x)
  
  # if weight is a strings
  if (is.character(weight)) {
    #weight is a row 
    if (weight %in% rownames(x)) {
      weight = (assays(x)[[assay]])[weight,]
    }
    #weight is a colData
    else if (weight %in% colnames(colData(x))) {
      weight = colData(x)[[weight]]
    }
    else {
      stop(paste("Cannot find",weight,'in neither rownames nor colData'))
    }
  }
  
  if (!is.null(weightTransform)) {
    if (is.vector(weightTransform)) {
      weight = weight  %in% weightTransform
    } 
    else if (is.function(weightTransform)) {
      tryCatch(
        {weight = weightTransform(weight)},
        error=function(e) {
          message("weightTransform must be a function that can accept the weight variable")
          print(e)
        }
      )
    }
  }
  
  #TODO: see if I can replace this with the ...
  hexDensity.default(x=xy,
                     xbins=xbins,
                     sigma=sigma,
                     edge=edge,
                     diggle = diggle,
                     weight=weight)
}

#' @rdname hexDensity
#' @export
hexDensity.default = function(x,y=NULL, 
                     xbins = 128, #128 is the magic number in spatstat
                     sigma = 1,
                     edge = TRUE,
                     diggle = FALSE,
                     weight = NULL) {
  hbin = hexbinFullRegular(x,y,xbins=xbins, weight=weight) 
  row = hbin@dimen[1]
  col = hbin@dimen[2]

  hexSize = diff(hbin@xbnds)/xbins
  
  #convert hexbin representation to staggered bin
  staggeredBin = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  for (i in seq(1,row,by=2)) {
    staggeredBin[i,(i-i%/%2):(i-i%/%2+col-1)] = hbin@count[((row-i)*col+1):((row-i)*col+col)]
    staggeredBin[i+1,(i-i%/%2):(i-i%/%2+col-1)] = hbin@count[((row-i-1)*col+1):((row-i-1)*col+col)]
  }
  #Make kernel
  
  # Older,slower kernel. Keeping just in case
  # kernel = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  # center.r = row
  # center.q = col+row/2
  # for (i in seq(1, row)) {
  #   for (j in seq(1,2*col)) {
  #     q = j+row-i
  #     kernel[i*2,q] = dnorm(hexSize*euclidDistance(center.q,center.r,q,2*i),sd=sigma) * dnorm(0,sd=sigma)
  #     kernel[i*2-1,q] = dnorm(hexSize*euclidDistance(center.q,center.r,q,2*i-1),sd=sigma) * dnorm(0,sd=sigma)
  #   }
  # }
  kernel.left.hori = dnorm(hexSize*c(seq(col,1),seq(0,col-1)),sd=sigma)
  kernel.left.verti = dnorm(hexSize*sqrt(3)*c(seq(row/2-1,0),seq(1,row/2)),sd=sigma)
  kernel.left=outer(kernel.left.verti,kernel.left.hori)

  kernel.right.hori = dnorm(hexSize*(c(seq(col,1)-0.5,seq(0,col-1)+0.5)),sd=sigma)
  kernel.right.verti = dnorm(hexSize*sqrt(3)*c(seq(row/2-1,0)+0.5,seq(1,row/2)-0.5),sd=sigma)
  kernel.right=outer(kernel.right.verti,kernel.right.hori)

  kernel = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  for (i in seq(1, row)) {
    kernel[i*2-1,(1+row-i):(1+row-i+2*col-1)] = rev(kernel.right[i,])
    kernel[i*2,(1+row-i):(1+row-i+2*col-1)] = rev(kernel.left[i,])
  }
  #This reverse the kernel
  kernel = kernel[,ncol(kernel):1]
  
  kernel = kernel/sum(kernel)

  #inverse the kernel
  kernel.inv = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  #going clock-wise from top-left of inverse kernel
  kernel.inv[1:(row+1),1:(col+row/2)] = kernel[row:(2*row),(col+row/2):(2*col+row-1)]
  kernel.inv[1:(row+1),((col+row/2)+1):(2*col+row-1)] = kernel[row:(2*row),1:(col+row/2)-1]
  kernel.inv[(row+2):(2*row),((col+row/2)+1):(2*col+row-1)] = kernel[1:(row-1),1:(col+row/2)-1]
  kernel.inv[(row+2):(2*row),1:(col+row/2)] = kernel[1:(row-1),(col+row/2):(2*col+row-1)]
  
  fK = fft2D(kernel.inv)

  #edge correction
  if (edge) {
    mask = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
    for (i in seq(1, row, by=2)) {
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

  #KDE calculation
  fY = fft2D(staggeredBin)
  sm = fft2D(fY*fK,inverse = TRUE)/(2*row*(2*col+(row-1)))

  #extract back to hexbin class
  count = c()
  if(edge && !diggle){
    #No need to clean up the NaN since will discard them anyway
    sm[1:row,1:(col+(row-1)/2)] = Re(sm[1:row,1:(col+(row-1)/2)])/edg
  }

  for (i in seq(row,1,by=-2)) {
    count = append(count,Re(sm[i,(i-i%/%2):(i-i%/%2+col-1)]))
    count = append(count,Re(sm[i-1,(i-i%/%2):(i-i%/%2+col-1)]))
  }

  hbin@count = count/hexAreaFromWidth(hexSize)
  return(hbin)
}




distance = function(q1,r1,q2,r2) {
  return((abs(q1-q2)
         + abs(q1+r1-q2-r2)
         +abs(r1-r2))/2)
}

euclidDistance = function(q1,r1,q2,r2) {
  x1 = sqrt(3)*q1 +(sqrt(3)/2) *r1
  y1 = (3/2)*r1
  
  x2 = sqrt(3)*q2 +(sqrt(3)/2) *r2
  y2 = (3/2)*r2
  return (sqrt(((x1-x2)**2 + (y1-y2)**2)/3))
}

hexAreaFromWidth = function(w) {
  return((w)**2*sqrt(3)/2)
}

