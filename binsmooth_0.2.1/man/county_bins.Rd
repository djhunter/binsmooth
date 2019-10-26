\name{county_bins}
\alias{county_bins}
\docType{data}
\title{
ACS County Income Data, 2006-2010
}
\description{
Binned income data from 3,221 counties in the U.S. and Puerto Rico.
}
\usage{data("county_bins")}
\format{
  A data frame with 51536 observations on the following 6 variables.
  \describe{
    \item{\code{fips}}{Number identifying the county}
    \item{\code{households}}{Bin counts}
    \item{\code{bin_min}}{Left endpoints of bins (US Dollars)}
    \item{\code{bin_max}}{Right endpoints of bins}
    \item{\code{county}}{County name}
    \item{\code{state}}{State name}
  }
}
\source{
U.S. Census Bureau, American Community Survey: \url{https://www.census.gov/programs-surveys/acs/}
}

\seealso{\code{\link{county_true}}}

\examples{
data(county_bins)
data(county_true)
binedges <- county_bins$bin_max[county_bins$fips=="6083"]+0.5 # continuity correction
bincounts <- county_bins$households[county_bins$fips=="6083"]
smean <- county_true$mean_true[county_true$fips=="6083"]
plot(splinebins(binedges, bincounts, smean)$splinePDF, 0, 300000,
     n=500, main="Santa Barbara County")
plot(stepbins(binedges, bincounts, smean)$stepPDF, do.points=FALSE, col="red", add=TRUE)
}
\keyword{datasets}
