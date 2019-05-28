\name{simcounty}
\alias{simcounty}
\title{
Simulate data to mimic \code{\link{county_bins}} and \code{\link{county_true}}
}
\description{
Samples from a selection of distributions (Gamma, Lognormal, Weibull, Triangle) to simulate income data in the
format used in the American Community Survey data (\code{\link{county_bins}} and \code{\link{county_true}}).
}
\usage{
simcounty(numCounties, minPop = 1000, maxPop = 100000,
          bin_minimums = c(0, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000,
                           50000, 60000, 75000, 100000, 125000, 150000, 200000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numCounties}{
The number of counties to simulate data for
}
  \item{minPop}{
Minimum population to sample (default = 1000)
}
  \item{maxPop}{
Maximum population to sample (default = 100000)
}
  \item{bin_minimums}{
Bin edges. Defaults to the edges used in the Census data.
}
}
\details{
The county names will tell which distributions were sampled to simulate each county.
}
\value{
Returns a list of two data frames:
\item{county_bins}{Simulated binned income data}
\item{county_true}{Statistics computed from the raw data}
}
\references{
Paul T. von Hippel, David J. Hunter, McKalie Drown. \emph{Better Estimates from Binned Income Data: Interpolated CDFs and Mean-Matching}, Sociological Science, November 15, 2017. \url{https://www.sociologicalscience.com/articles-v4-26-641/}
}
\author{
David J. Hunter and McKalie Drown
}
\seealso{
\code{\link{county_bins}}, \code{\link{county_true}}
}
\examples{
l <- simcounty(5)
cb <- l$county_bins
ct <- l$county_true
sbl <- splinebins(cb$bin_max[cb$fips==103], cb$households[cb$fips==103],
                  ct$mean_true[ct$fips==103])
stl <- stepbins(cb$bin_max[cb$fips==105], cb$households[cb$fips==105],
                ct$mean_true[ct$fips==105])
plot(sbl$splinePDF, 0, 300000, n=500)
plot(stl$stepPDF, do.points=FALSE, main=cb$county[cb$fips==105][1])
}
