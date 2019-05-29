stats_from_distribution <- function(binFit) {
  PDF <- binFit[[1]]
  CDF <- binFit[[2]]
  E <- binFit[[3]]
  cdf_mean <- E - pracma::integral(CDF, 0, E)
  v <- pracma::integral(function(x){2*x-2*x*CDF(x)}, 0 ,E) - cdf_mean^2
  g <- 1-pracma::integral(function(x){(1-CDF(x))^2}, 0, E)/cdf_mean
  t <- pracma::integral(function(x){PDF(x)*x/cdf_mean*log(x/cdf_mean)}, 0, E)
  statistics <- c(cdf_mean, v, sqrt(v), g, t)
  names(statistics) <- c("mean", "variance", "SD", "Gini", "Theil")
  return(statistics)
}