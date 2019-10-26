\name{sb_percentiles}
\alias{sb_percentiles}
\title{
Estimate percentiles from splinebins
}
\description{
Estimates percentiles of a smoothed distribution obtained using \code{\link{splinebins}}.
}
\usage{
sb_percentiles(splinebinFit, p = seq(0,100,25))
}
\arguments{
  \item{splinebinFit}{
A list as returned by \code{\link{splinebins}}.
}
  \item{p}{
A vector of percentages in the range \eqn{0 \le p \le 100}.
}
}
\details{
The approximate inverse of the CDF calculated by \code{\link{splinebins}} is used to approximate percentiles of the smoothed distribution.
}
\value{
A vector of percentiles.
}
\references{
Paul T. von Hippel, David J. Hunter, McKalie Drown. \emph{Better Estimates from Binned Income Data: Interpolated CDFs and Mean-Matching}, Sociological Science, November 15, 2017. \url{https://www.sociologicalscience.com/articles-v4-26-641/}
}
\author{
David J. Hunter and McKalie Drown
}

\examples{
# 2005 ACS data from Cook County, Illinois
binedges <- c(10000,15000,20000,25000,30000,35000,40000,45000,
              50000,60000,75000,100000,125000,150000,200000,NA)
bincounts <- c(157532,97369,102673,100888,90835,94191,87688,90481,
               79816,153581,195430,240948,155139,94527,92166,103217)
splinefit <- splinebins(binedges, bincounts, 76091)
sb_percentiles(splinefit)
sb_percentiles(splinefit, c(27, 32, 93))
}

