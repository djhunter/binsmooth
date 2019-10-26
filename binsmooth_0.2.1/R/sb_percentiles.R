sb_percentiles <- function(splinebinFit, p = seq(0,100,25)) {
  iCDF <- splinebinFit$splineInvCDF
  percentiles <- iCDF(p/100)
  names(percentiles) <- paste0(p, "%")
  return(percentiles)
}