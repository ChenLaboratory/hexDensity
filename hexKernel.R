library(spatstat.data)
library(spatstat.explore)
source("hexbinFullRegular.R")

hexKernel = function(x,y=NULL, 
                     xbins = 128, 
                     sigma = 1,
                     edge = TRUE) {
  hbin = hexbinFullRegular(x,y=NULL,xbins=xbins) #128 is the magic number in pixellate of spatstat
  row = hbin@dimen[1]
  col = hbin@dimen[2]
  print(row)
  print(col)
  hexSize = diff(hbin@xbnds)/xbins
  
  #convert hexbin representation to staggered bin
  staggeredBin = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  for (i in seq(1, row, by=2)) {
    staggeredBin[i:(i+1),((row-i+1)/2):((row-i+1)/2+col-1)] = apply(matrix(rev(rev(hbin@count)[((i-1)*col+1):(col*(i+1))]),2,col,byrow=TRUE),2,rev)
  }
  
  #Make kernel
  kernel = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  center.r = row
  center.q = col+row/2
  for (i in seq(1, row)) {
    for (j in seq(1,2*col)) {
      q = j+row-i
      kernel[i*2,q] = dnorm(hexSize*euclidDistance(center.q,center.r,q,2*i),sd=sigma) * dnorm(0,sd=sigma)
      kernel[i*2-1,q] = dnorm(hexSize*euclidDistance(center.q,center.r,q,2*i-1),sd=sigma) * dnorm(0,sd=sigma)
    }
  } 
  #print(sum(kernel))
  kernel = kernel/sum(kernel)
  
  #inverse the kernel
  kernel.inv = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
  #going clock-wise from top-left of inverse kernel
  kernel.inv[1:(row+1),1:(col+row/2)] = kernel[row:(2*row),(col+row/2):(2*col+row-1)]
  kernel.inv[1:(row+1),((col+row/2)+1):(2*col+row-1)] = kernel[row:(2*row),1:(col+row/2)-1]
  kernel.inv[(row+2):(2*row),((col+row/2)+1):(2*col+row-1)] = kernel[1:(row-1),1:(col+row/2)-1]
  kernel.inv[(row+2):(2*row),1:(col+row/2)] = kernel[1:(row-1),(col+row/2):(2*col+row-1)]
  
  fK = fft2D(kernel.inv)
  fY = fft2D(staggeredBin)
  sm = fft2D(fY*fK,inverse = TRUE)/(2*row*(2*col+row-1))

  if (edge) {
    mask = matrix(0,nrow = 2*row, ncol = 2*col+row-1)
    for (i in seq(1, row, by=2)) {
      mask[i:(i+1),((row-i+1)/2):((row-i+1)/2+col-1)] = matrix(1,2,col)
    }
    fM = fft2D(mask)
    con = fft2D(fM * fK, inverse=TRUE)
    con = con/(2*row*(2*col+row-1))
    edg <- Mod(con[1:row, 1:(col+(row-1)/2)])
    print(edg)
  }
  
  #extract back to hexbin class
  count = c()
  if(edge){
    for (i in seq(row,1,by=-2)) {
      count = append(count,Re(sm[i,(1+(row-i)/2):(col+(row-i)/2)]/edg[i,(1+(row-i)/2):(col+(row-i)/2)]))
      count = append(count,Re(sm[i-1,(1+(row-i)/2):(col+(row-i)/2)]/edg[i-1,(1+(row-i)/2):(col+(row-i)/2)]))
    }
  } else{
    for (i in seq(row,1,by=-2)) {
      count = append(count,Re(sm[i,(1+(row-i)/2):(col+(row-i)/2)]))
      count = append(count,Re(sm[i-1,(1+(row-i)/2):(col+(row-i)/2)]))
    }
  }

  hbin@count = count
  #print("done")
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
  
  return (sqrt((x1-x2)**2 + (y1-y2)**2))
}
