simcounty <- function(numCounties, minPop=1000, maxPop=100000,
                      bin_minimums=c(0,10000,15000,20000,25000,30000,35000,40000,45000,
                                     50000,60000,75000,100000,125000,150000,200000)) {
  # Preallocate space for the whole data frames
  numBins <- length(bin_minimums)
  # variables for county_true data frame
  fips <- 100+1:numCounties
  mean_true <- numeric(numCounties)
  median_true <- numeric(numCounties)
  gini_true <- numeric(numCounties)
  theil_true <- numeric(numCounties)
  numRows <- numBins*numCounties
  #variables for county_bins data frame
  # later: fips <- rep(fips, each=numBins)
  households <- numeric(numRows)
  bin_min <- rep(bin_minimums, numCounties)
  bin_max <- rep(c(bin_minimums[2:numBins]-1,NA), numCounties)
  county <- character(numRows)
  state <- rep(" Simulated", numRows)
  minPop <- 1000
  maxPop <- 100000
  dSamps <- list(numeric(maxPop), numeric(maxPop), numeric(maxPop)) # 3 samples for each county
  distributions <- c("lognormal", "triangle", "weibull", "gamma")
  for(countyI in 1:numCounties) {
    distsToUse <- sample(distributions, 3, replace=TRUE)
    i <- 1
    for(d in distsToUse)
    {
      sSize <- sample(maxPop,1)+minPop
      switch(d,
             lognormal = {
               lnMu <- sample(c(10.5,10.6,10.7,10.8,10.9),1)
               lnSigma <- sample(c(0.6,0.7,0.8),1)
               dSamps[[i]] <- rlnorm(sSize, lnMu, lnSigma)
             },
             triangle = {
               triB <- sample(c(300000, 400000, 500000),1)
               triC <- sample(c(10000, 15000, 20000, 25000),1)
               dSamps[[i]] <- triangle::rtriangle(sSize,a=0,b=triB,c=triC)
             },
             weibull = {
               weibA <- sample(c(1.4, 1.5, 1.6, 1.7),1)
               weibB <- sample(c(50000, 60000, 70000, 80000),1)
               dSamps[[i]] <-  rweibull(sSize,weibA, weibB)
             },
             gamma = {
               gammaA <- sample(seq(1,3,by=0.2),1)
               dSamps[[i]] <- rgamma(sSize,gammaA,gammaA/50000)
             }
      )
      i <- i+1
    }
    simPopSamp <- unlist(dSamps)
    households[(1+(countyI-1)*numBins):(numBins*countyI)] <- c(vapply(1:(numBins-1),
                                                                      function(x) {sum(simPopSamp >= bin_minimums[x] & simPopSamp < bin_minimums[x+1])},
                                                                      numeric(1)),
                                                               sum(simPopSamp>bin_minimums[numBins]))
    county[(1+(countyI-1)*numBins):(numBins*countyI)] <- rep(paste(distsToUse, collapse=" "), numBins)
    mean_true[countyI] <- mean(simPopSamp)
    median_true[countyI] <- median(simPopSamp)
    gini_true[countyI] <- ineq::Gini(simPopSamp, corr=TRUE)
    theil_true[countyI] <- ineq::Theil(simPopSamp)
    countyI <- countyI + 1
  } #end of for countyI loop
  county_true <- data.frame(fips, mean_true, median_true, gini_true, theil_true)
  fips <- rep(fips, each=numBins)
  county_bins <- data.frame(fips,households,bin_min,bin_max,county,state)
  return(list(county_bins=county_bins, county_true=county_true))
}
