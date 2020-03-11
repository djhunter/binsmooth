\name{county_true}
\alias{county_true}
\docType{data}
\title{
ACS County Income Statistics, 2006-2010
}
\description{
Statistics computed from raw data on 3,221 counties in the U.S. and Puerto Rico.
}
\usage{data("county_true")}
\format{
  A data frame with 3221 observations on the following 4 variables.
  \describe{
    \item{\code{fips}}{Number identifying the county}
    \item{\code{mean_true}}{Sample mean}
    \item{\code{median_true}}{Sample median}
    \item{\code{gini_true}}{Gini coefficient}
  }
}
\source{
U.S. Census Bureau, American Community Survey: \url{https://www.census.gov/programs-surveys/acs/}
}

\seealso{\code{\link{county_bins}}}

\examples{
data(county_bins)
data(county_true)
binedges <- county_bins$bin_max[county_bins$fips=="6083"]+0.5 # continuity correction
bincounts <- county_bins$households[county_bins$fips=="6083"]
smean <- county_true$mean_true[county_true$fips=="6083"]
plot(stepbins(binedges, bincounts, smean)$stepPDF, do.points=FALSE,
     main="Santa Barbara County")
}
\keyword{datasets}
