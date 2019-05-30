splinebins <- function(bEdges, bCounts, m=NULL, numIterations=16, monoMethod=c("hyman", "monoH.FC")) {
  monoMethod <- match.arg(monoMethod)
  L <- length(bCounts)
  tot <- sum(bCounts)
  if(is.null(m))  # no mean provided, so make one up
    m <- sum(0.5*(c(bEdges[1:(L-1)],2.0*bEdges[L-1])+c(0, bEdges[1:(L-1)]))*bCounts/tot)
  sumbAreas <- vapply(1:L, function(i) {sum(bCounts[1:i])/tot}, numeric(1))
  tailEnd <- 1.05*bEdges[L-1] # start with a really short tail
  x <- c(0, bEdges[1:(L-1)], tailEnd, tailEnd*1.01)
  y <- c(0, sumbAreas, 1) # The last x,y pair forces smoothness at x=tailEnd
  f <- splinefun(x,y, method=monoMethod)
  est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
  shrinkFactor <- 1
  shrinkMultiplier <- 0.995
  while(est_mean > m) { # binning must be skewed, because supplied mean is too small
    x <- x*shrinkMultiplier # shrink the bins
    tailEnd <- tailEnd*shrinkMultiplier
    f <- splinefun(x,y, method=monoMethod)
    est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
    shrinkFactor <- shrinkFactor*shrinkMultiplier
  }
  if(bCounts[L] > 0) {
    # search for tail length
    while(est_mean < m) {
      tailEnd <- tailEnd * 2
      x[L+1] <- tailEnd
      x[L+2] <- tailEnd*1.01
      f <- splinefun(x,y, method=monoMethod)
      est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
    }
    # now the correct tail end will be between bEdges[L-1] and tailEnd
    l <- bEdges[L-1]
    r <- tailEnd
    for(i in 1:numIterations) {
      tailEnd <- (l+r)/2
      x[L+1] <- tailEnd
      x[L+2] <- tailEnd*1.01
      f <- splinefun(x,y, method=monoMethod)
      est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
      if(est_mean > m)
        r <- tailEnd
      else
        l <- tailEnd
    }
  }
  xfix <- seq(0, tailEnd, length.out = 200) # ADDED in v0.2.0: sample CDF to invert
  yfix <- f(xfix)
  finv <- splinefun(yfix, xfix, method=monoMethod) # ADDED in v0.2.0: approximate inverse CDF
  splineCDF <- function(x){
    ifelse(x<0, 0, ifelse(x>tailEnd, 1, f(x)))
  }
  splinePDF <- function(x){
    ifelse(x<0 | x>tailEnd, 0, f(x,deriv=1))
  }
  splineInvCDF <- function(x){ # ADDED in v0.2.0: approximate inverse CDF
    ifelse(x<0, 0, ifelse(x>1, tailEnd, finv(x)))
  }
  return(list(splinePDF=splinePDF, splineCDF=splineCDF, E=tailEnd, est_mean=est_mean, shrinkFactor=shrinkFactor, splineInvCDF=splineInvCDF))
}
