rsubbins <- function(bEdges, bCounts, m=NULL, eps1 = 0.25, eps2 = 0.75, depth = 3,
                     tailShape = c("onebin", "pareto", "exponential"),
                     nTail=16, numIterations=20, pIndex=1.160964, tbRatio=0.8) {
  if(is.null(m)) { # no mean provided, so make one up
    L <- length(bCounts)
    m <- sum(0.5*(c(bEdges[1:(L-1)],2.0*bEdges[L-1])+c(0, bEdges[1:(L-1)]))*bCounts/sum(bCounts))
  }
  tailShape <- match.arg(tailShape)
  if(tailShape == "onebin")
    rsubbinsNotail(bEdges, bCounts, m, eps1, eps2, depth)
  else
    rsubbinsTail(bEdges, bCounts, m, eps1, eps2, depth, tailShape, nTail, numIterations, pIndex, tbRatio)
}

rsubbinsTail <- function(bEdges, bCounts, m, eps1, eps2, depth,
                         tailShape = c("pareto", "exponential"),
                         nTail, numIterations, pIndex, tbRatio) {
  tailShape <- match.arg(tailShape)
  eps3 <- (1 - eps2) / 2
  L <- length(bCounts)
  tot <- sum(bCounts)
  e <- c(0,bEdges[1:(L-1)],numeric(nTail)) # preallocate tail edges
  shrinkFactor <- 1
  shrinkMultiplier <- 0.995
  tailCount <- bCounts[L]
  bbtot <- tot-tailCount
  bbMean <- sum((e[2:(L)]+e[1:(L-1)])*bCounts[1:(L-1)])/(2*bbtot) # mean of bounded bins as a step function
  while(m<bbMean) { # binning must be skewed, because supplied mean is less than mean of bounded bins
    e <- e*shrinkMultiplier # shrink the bins
    bbMean <- sum((e[2:(L)]+e[1:(L-1)])*bCounts[1:(L-1)])/(2*bbtot)
    shrinkFactor <- shrinkFactor*shrinkMultiplier
  }
  if(tailCount>0) {
    # there is remaining area in final bin, so make a tail of bins
    L <- L+nTail-1
    tailArea <- tailCount/tot
    bAreas <- c(bCounts/tot, numeric(nTail-1))
    tbWidth <- e[L-nTail+1]-e[L-nTail] # initial tail bin width = width of last bounded bin
    h <- bAreas[L-nTail]/tbWidth # height of last bounded bin
    if(tailShape=="pareto")
      tailUnscaled <- (1:nTail)^(-1-pIndex)
    else
      tailUnscaled <- tbRatio^(1:nTail)
    bAreas[(L-nTail+1):(L)] <- tailUnscaled/sum(tailUnscaled)*tailArea
    repeat {
      e[(L-nTail+2):(L+1)] <- e[L-nTail+1]+(1:(nTail))*tbWidth
      bMean <- sum((e[2:(L+1)]+e[1:L])*bAreas)/2
      if(bMean>m) break # now the correct tbWidth is between 0 and tbWidth
      tbWidth <- tbWidth*2
    }
    l <- 0
    r <- tbWidth
    for(i in 1:numIterations) {
      tbWidth <- (l+r)/2
      e[(L-nTail+2):(L+1)] <- e[L-nTail+1]+(1:(nTail))*tbWidth
      bMean <- sum((e[2:(L+1)]+e[1:L])*bAreas)/2
      if (bMean<m) # bin mean is too small; lengthen bins
        l <- tbWidth
      else
        r <- tbWidth
    }
  }
  else { # final bin is empty, so get rid of it
    L <- L-1
    bAreas <- bCounts[1:L]/tot
    e <- e[1:(L+1)]
  }
  # subdivide bins
  binAreas <- bAreas
  binHeights <- c(0, binAreas/(e[2:(L+1)]-e[1:L]), 0)
  binEdges <- e
  for(i in 1:depth){
    L <- length(binEdges)
    binDiffs <- binHeights[2:(L + 1)] - binHeights[1:L]
    newbinHeights <- 1:(3 * (L - 1))
    binWidths <- binEdges[2:L] - binEdges[1:(L-1)]
    binEdges2 <- binEdges[1:(L - 1)] + binWidths * eps3
    binEdges3 <- binEdges[2:L] - binWidths * eps3
    binEdges <- sort(c(binEdges, binEdges2, binEdges3), decreasing = FALSE)
    for(j in 1:(L-1)) {
      new_middlebin_height <- (((binDiffs[j] - binDiffs[(j + 1)]) * eps1 * eps3) / eps2) + binHeights[j + 1]
      if(new_middlebin_height>0)
      {
        newbinHeights[(3 * j - 2)] <- binHeights[j + 1] - binDiffs[j] * eps1
        newbinHeights[(3 * j)] <- binHeights[(j + 1)] + binDiffs[(j + 1)] * eps1
        newbinHeights[(3 * j - 1)] <- new_middlebin_height
      }
      else # don't shift when middle bin would go negative
      {
        newbinHeights[(3*j-2)] <- binHeights[(j+1)]
        newbinHeights[(3*j)] <- binHeights[(j+1)]
        newbinHeights[(3*j-1)] <- binHeights[(j+1)]
      }
    }
    binHeights <- c(0, newbinHeights, 0)
  }
  L <- length(binEdges)
  binWidths <- binEdges[2:L] - binEdges[1:(L-1)]
  binAreas <- binWidths * binHeights[2:L]
  cAreas <- vapply(1:length(binAreas), function(x){sum(binAreas[1:x])}, numeric(1))
  rsubCDF <- approxfun(binEdges, c(0, cAreas), yleft=0, yright=1, rule=2)
  rsubPDF <- stepfun(binEdges, binHeights)
  return(list(rsubPDF=rsubPDF, rsubCDF=rsubCDF, E=binEdges[L], shrinkFactor=shrinkFactor))
}

rsubbinsNotail <- function(bEdges, bCounts, m, eps1, eps2, depth) {
  eps3 <- (1 - eps2) / 2
  L <- length(bCounts)
  tot <- sum(bCounts)
  p <- c(1,1-vapply(1:(L-1), function(x){sum(bCounts[1:x])}, numeric(1))/tot)
  e <- c(0,bEdges[1:(L-1)],0) # preallocate last entry; replace with E later
  A <- 0.5*sum((e[2:L]-e[1:(L-1)])*(p[2:L]+p[1:(L-1)]))
  shrinkFactor <- 1
  shrinkMultiplier <- 0.995
  while(m<A) { # binning must be skewed, because supplied mean is too small
    e <- e*shrinkMultiplier # shrink the bins
    A <- 0.5*sum((e[2:L]-e[1:(L-1)])*(p[2:L]+p[1:(L-1)]))
    shrinkFactor <- shrinkFactor*shrinkMultiplier
  }
  E <- ifelse(p[L]>0, bEdges[L-1]+2*(m-A)/p[L], bEdges[L-1]*1.001) # make width of final bin small if empty
  e[L+1] <- E # this is the right border of the last bin so the mean is preserved
  # subdivide bins
  binAreas <- bCounts/tot
  binHeights <- c(0, binAreas/(e[2:(L+1)]-e[1:L]), 0)
  binEdges <- e
  for(i in 1:depth){
    L <- length(binEdges)
    binDiffs <- binHeights[2:(L + 1)] - binHeights[1:L]
    newbinHeights <- 1:(3 * (L - 1))
    binWidths <- binEdges[2:L] - binEdges[1:(L-1)]
    binEdges2 <- binEdges[1:(L - 1)] + binWidths * eps3
    binEdges3 <- binEdges[2:L] - binWidths * eps3
    binEdges <- sort(c(binEdges, binEdges2, binEdges3), decreasing = FALSE)
    for(j in 1:(L-1)) {
      new_middlebin_height <- (((binDiffs[j] - binDiffs[(j + 1)]) * eps1 * eps3) / eps2) + binHeights[j + 1]
      if(new_middlebin_height>0)
      {
        newbinHeights[(3 * j - 2)] <- binHeights[j + 1] - binDiffs[j] * eps1
        newbinHeights[(3 * j)] <- binHeights[(j + 1)] + binDiffs[(j + 1)] * eps1
        newbinHeights[(3 * j - 1)] <- new_middlebin_height
      }
      else # don't shift when middle bin would go negative
      {
        newbinHeights[(3*j-2)] <- binHeights[(j+1)]
        newbinHeights[(3*j)] <- binHeights[(j+1)]
        newbinHeights[(3*j-1)] <- binHeights[(j+1)]
      }
    }
    binHeights <- c(0, newbinHeights, 0)
  }
  L <- length(binEdges)
  binWidths <- binEdges[2:L] - binEdges[1:(L-1)]
  binAreas <- binWidths * binHeights[2:L]
  cAreas <- vapply(1:length(binAreas), function(x){sum(binAreas[1:x])}, numeric(1))
  rsubCDF <- approxfun(binEdges, c(0, cAreas), yleft=0, yright=1, rule=2)
  rsubPDF <- stepfun(binEdges, binHeights)
  return(list(rsubPDF=rsubPDF, rsubCDF=rsubCDF, E=E, shrinkFactor=shrinkFactor))
}
