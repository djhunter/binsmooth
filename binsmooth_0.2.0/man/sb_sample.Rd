\name{sb_sample}
\alias{sb_sample}
\title{
Random sample from splinebins distribution
}
\description{
Draw a random sample of points from a smoothed distribution obtained using \code{\link{splinebins}}.
}
\usage{
sb_sample(splinebinFit, n = 1)
}
\arguments{
  \item{splinebinFit}{
A list as returned by \code{\link{splinebins}}.
}
  \item{n}{
A positive integer giving the sample size.
}
}
\details{
The approximate inverse of the CDF calculated by \code{\link{splinebins}} is used to generate random values of the smoothed distribution.
}
\value{
A vector of random deviates.
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
sb_sample(splinefit, 5)
hist(sb_sample(splinefit, 3000))
}

