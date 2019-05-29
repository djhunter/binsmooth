\name{stats_from_distribution}
\alias{stats_from_distribution}
\title{
Estimate various statistics
}
\description{
Estimates the mean, variance, standard deviation, Gini coefficient, and Theil index from a smoothed distribution.
}
\usage{
stats_from_distribution(binFit)
}
\arguments{
  \item{binFit}{
A list as returned by \code{\link{splinebins}}, \code{\link{stepbins}}, or \code{\link{rsubbins}}. (Alternatively, a list containing a PDF of non-negative support, its CDF, and an upper bound for the support of the PDF.)
}
}
\details{
The mean and variance are calculated from the CDF. For details on the other statistics, see \code{\link{gini}} and \code{\link{theil}}.
}
\value{
A vector of five statistics.
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
stepfit <- stepbins(binedges, bincounts, 76091)
splinefit <- splinebins(binedges, bincounts, 76091)
stats_from_distribution(stepfit)
stats_from_distribution(splinefit) # More accurate
}

