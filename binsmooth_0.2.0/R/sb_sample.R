sb_sample <- function(splinebinFit, n = 1) {
  iCDF <- splinebinFit$splineInvCDF
  return(iCDF(stats::runif(n)))
}
