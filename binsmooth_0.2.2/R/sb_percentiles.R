sb_percentiles <- function(splinebinFit, p = seq(0,100,25)) {
  if(splinebinFit$fitWarn) {
    warning("Inaccurate fit detected. Returning NA.\n")
    return(NA)
  }
  iCDF <- splinebinFit$splineInvCDF
  percentiles <- iCDF(p/100)
  names(percentiles) <- paste0(p, "%")
  return(percentiles)
}