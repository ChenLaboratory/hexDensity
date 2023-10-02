library(hexbin)
library(spatstat.data)
library(spatstat.explore)

hexKernel = function(x, sigma = 1) {
  hbin = hexbin(x,xbins=128) #128 is the magic number in pixellate of spatstat
  row = hbin@dimen[1]
  col = hbin@dimen[2]
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
      kernel[i*2,q] = dnorm(distance(center.q,center.r,q,2*i),sd=sigma) * dnorm(0,sd=sigma)
      kernel[i*2-1,q] = dnorm(distance(center.q,center.r,q,2*i-1),sd=sigma) * dnorm(0,sd=sigma)
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
  #sm = convolve(staggeredBin,kernel.inv)#/(2*row*(2*col+row-1))

  #extract back to hexbin class
  count = c()
  for (i in seq(row,1,by=-2)) {
    count = Re(append(count,sm[i,(1+(row-i)/2):(col+(row-i)/2)]))
    count = Re(append(count,sm[i-1,(1+(row-i)/2):(col+(row-i)/2)]))
  }
  #print(count)
  hbin@count = count
  
  return(hbin)
  
}
  
distance = function(q1,r1,q2,r2) {
  return((abs(q1-q2)
         + abs(q1+r1-q2-r2)
         +abs(r1-r2))/2)
}
data = bei
bins = hexbin(data$x, data$y,xbins=13)
b=density.ppp(data,sigma=1)
b2=density.ppp(data,sigma=10)
plot(b)
plot(b2)
plot(density.ppp(data,sigma=100))
hexbin::plot(bins)
a = pixellate(bei,padzero = TRUE)
plot(data)
plot(bins)

test = hexKernel(data)
test2 = hexKernel(data,sigma=3)
plot(test)
plot(test2)


nr=6
nc=5
xcol.ker = c(0:(nc-1),-(nc:1))
yrow.ker = c(0:(nr-1),-(nr:1))
densX.ker <- dnorm(xcol.ker, sd=1)
densY.ker <- dnorm(yrow.ker, sd=1)
Kern <- outer(densY.ker, densX.ker, "*")
Kern = Kern/sum(Kern)

