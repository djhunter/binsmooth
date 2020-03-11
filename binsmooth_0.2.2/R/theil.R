theil <- function(binFit) {
  PDF <- binFit[[1]]
  CDF <- binFit[[2]]
  E <- binFit[[3]]
  cdf_mean <- E - pracma::integral(CDF, 0, E)
  return(pracma::integral(function(x){PDF(x)*x/cdf_mean*log(x/cdf_mean)}, 0, E))
}