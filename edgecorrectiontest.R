library(spatstat.data)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.utils)
library(plotly)

mydensity.ppp <- function(x, sigma=NULL, ...,
                        weights=NULL, 
                        edge=TRUE, varcov=NULL,
                        at="pixels", leaveoneout=TRUE,
                        adjust=1, diggle=FALSE,
                        se=FALSE, wtype=c("value", "multiplicity"),
                        kernel="gaussian",
                        scalekernel=is.character(kernel),
                        positive=FALSE, verbose=TRUE) {
  verifyclass(x, "ppp")
  
  output <- pickoption("output location type", at,
                       c(pixels="pixels",
                         points="points"))
  
  if(any(sidelengths(Frame(x)) == 0)) { ## pixels will have zero area
    val <- npoints(x)/0 # Inf or NaN
    return(as.im(val, W=Frame(x), ...)) 
  }
  
  if(!identical(kernel, "gaussian")) {
    validate2Dkernel(kernel)
    if(verbose && scalekernel &&
       (is.function(sigma) || (is.null(sigma) && is.null(varcov))))
      warning("Bandwidth selection will be based on Gaussian kernel")
  }
  
  ker <- resolve.2D.kernel(..., sigma=sigma, varcov=varcov, x=x, adjust=adjust)
  sigma <- ker$sigma
  varcov <- ker$varcov
  ## sigma.is.infinite <- ker$infinite
  
  if(is.im(weights)) {
    weights <- safelookup(weights, x) # includes warning if NA
  } else if(is.expression(weights)) 
    weights <- eval(weights, envir=as.data.frame(x), enclos=parent.frame())
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL 
  
  if(se) {
    ## compute standard error
    wtype <- match.arg(wtype)
    SE <- denspppSEcalc(x, sigma=sigma, varcov=varcov,
                        kernel=kernel,
                        ...,
                        weights=weights, wtype=wtype, edge=edge, at=output,
                        leaveoneout=leaveoneout, adjust=adjust,
                        diggle=diggle)
    if(positive) SE <- posify(SE)
  }
  
  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nx <- npoints(x)
    single <- is.null(dim(weights))
    totwt <- if(is.null(weights)) nx else
      if(single) sum(weights) else colSums(weights)
    if(!edge) totwt <- 0 * totwt
    W <- Window(x)
    A <- area.owin(W)
    switch(output,
           pixels = {
             E <- solapply(totwt/A, as.im, W=W, ...)
             names(E) <- colnames(weights)
             if(single) E <- E[[1L]]
           },
           points = {
             numerator <- rep(totwt, each=nx)
             if(!single) numerator <- matrix(numerator, nrow=nx)
             if(leaveoneout && edge) 
               numerator <- numerator - (weights %orifnull% 1)
             E <- numerator/A
             if(!single)
               colnames(E) <- colnames(weights)
           })
    result <- if(se) list(estimate=E, SE=SE) else E
    return(result)
  }
  
  if(output == "points") {
    # VALUES AT DATA POINTS ONLY
    result <- densitypointsEngine(x, sigma,
                                  varcov=varcov,
                                  kernel=kernel,
                                  scalekernel=scalekernel,
                                  weights=weights, edge=edge,
                                  leaveoneout=leaveoneout,
                                  diggle=diggle, ...)
    if(verbose && !is.null(uhoh <- attr(result, "warnings"))) {
      switch(uhoh,
             underflow=warning("underflow due to very small bandwidth"),
             warning(uhoh))
    }
    ## constrain values to be positive
    if(positive) 
      result <- posify(result)
    if(se) 
      result <- list(estimate=result, SE=SE)
    return(result)
  }
  
  # VALUES AT PIXELS
  if(!edge) {
    # no edge correction
    edg <- NULL
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              kernel=kernel,
                              scalekernel=scalekernel,
                              weights=weights, varcov=varcov)
    raw <- divide.by.pixelarea(raw) 
    smo <- raw
  } else if(!diggle) {
    # edge correction e(u)
    both <- second.moment.calc(x, sigma, what="smoothedge", ...,
                               kernel=kernel,
                               scalekernel=scalekernel,
                               weights=weights, varcov=varcov)
    raw <- divide.by.pixelarea(both$smooth)
    edg <- both$edge
    ## Math.im / Math.imlist not yet working
    smo <- imagelistOp(raw, edg, "/")
  } else {
    # edge correction e(x_i)
    edg <- second.moment.calc(x, sigma, what="edge", ...,
                              scalekernel=scalekernel,
                              kernel=kernel, varcov=varcov)
    return(edg)
    wi <- 1/safelookup(edg, x, warn=FALSE)
    wi[!is.finite(wi)] <- 0
    # edge correction becomes weight attached to points
    if(is.null(weights)) {
      newweights <- wi
    } else if(is.matrix(weights) || is.data.frame(weights)) {
      stopifnot(nrow(weights) == npoints(x))
      newweights <- weights * wi
    } else {
      stopifnot(length(weights) == npoints(x))
      newweights <- weights * wi
    }
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              kernel=kernel,
                              scalekernel=scalekernel,
                              weights=newweights, varcov=varcov)
    raw <- divide.by.pixelarea(raw)
    smo <- raw
  }
  
  result <- if(is.im(smo)) smo[x$window, drop=FALSE]
  else solapply(smo, "[", i=x$window, drop=FALSE)
  
  # internal use only
  spill <- resolve.1.default(list(spill=FALSE), list(...))
  if(spill)
    return(list(result=result, sigma=sigma, varcov=varcov, raw = raw, edg=edg))
  
  # constrain values to be positive
  if(positive) 
    result <- posify(result)
  
  # normal return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  attr(result, "kernel") <- kernel
  if(se)
    result <- list(estimate=result, SE=SE)
  return(result)
}


divide.by.pixelarea <- function(x) {
  if(is.im(x)) {
    x$v <- x$v/(x$xstep * x$ystep)
  } else {
    for(i in seq_along(x))
      x[[i]]$v <- with(x[[i]], v/(xstep * ystep))
  }
  return(x)
}

x = bei
sigma=40
edg = mydensity.ppp(x,sigma = 100000,diggle=TRUE,eps=diff(range(data$x))/128)
plot_ly(z=edg$v) %>% add_surface()
safelookup = safelookup(edg, x, warn=FALSE)

plot(test)
test_no_correct = mydensity.ppp(a,sigma=sigma,edge=FALSE)
plot(test_no_correct)

obswin = as.owin(x)
X= pixellate(x,padzero = TRUE)
xstep <- X$xstep
ystep <- X$ystep
Y <- X$v
nr <- nrow(Y)
nc <- ncol(Y)
lengthYpad <- 4 * nc * nr
xcol.ker <- xstep *  c(0:(nc-1),-(nc:1))
yrow.ker <- ystep *  c(0:(nr-1),-(nr:1))
densX.ker <- dnorm(xcol.ker, sd=sigma)
densY.ker <- dnorm(yrow.ker, sd=sigma)
Kern <- outer(densY.ker, densX.ker, "*")
Kern <- Kern/sum(Kern)
plot_ly(z=Kern) %>% add_surface()

fK <- fft2D(Kern)
M <- as.mask(obswin, xy=list(x=X$xcol, y=X$yrow))$m
Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
Mpad[1:nr, 1:nc] <- M
plot_ly(z=Mpad) %>% add_surface()

xcol.pad <- X$xcol[1] + xstep * (0:(2*nc-1))
yrow.pad <- X$yrow[1] + ystep * (0:(2*nr-1))

lengthMpad <- 4 * nc * nr
fM <- fft2D(Mpad)
con <- fft2D(fM * fK, inverse=TRUE)/lengthMpad
edg <- Mod(con)
plot_ly(z=edg) %>% add_surface()

edg = edg[1:nr,1:nc]
edg <- im(edg, xcol.pad[1:nc], yrow.pad[1:nr])
(edg$v)[32,32]
