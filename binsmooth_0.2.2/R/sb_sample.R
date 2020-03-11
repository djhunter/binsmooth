sb_sample <- function(splinebinFit, n = 1) {
  if(splinebinFit$fitWarn) {
    warning("Inaccurate fit detected. Returning NA.\n")
    return(NA)
  }
  iCDF <- splinebinFit$splineInvCDF
  return(iCDF(stats::runif(n)))
}
