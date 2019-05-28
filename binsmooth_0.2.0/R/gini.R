gini <- function(binFit) {
  CDF <- binFit[[2]]
  E <- binFit[[3]]
  cdf_mean <- E - pracma::integral(CDF, 0, E)
  return(1-pracma::integral(function(x){(1-CDF(x))^2}, 0, E)/cdf_mean)
}