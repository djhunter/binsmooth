\name{splinebins}
\alias{splinebins}
\title{
Optimized spline PDF and CDF fitted to binned data
}
\description{
Creates a smooth cubic spline CDF and piecewise-quadratic PDF based on a set of binned data (edges and counts).
}
\usage{
splinebins(bEdges, bCounts, m = NULL,
           numIterations = 16, monoMethod = c("hyman", "monoH.FC"))
}
\arguments{
  \item{bEdges}{
A vector \eqn{e_1, e_2, \ldots, e_n} giving the right endpoints of each bin. The value in \eqn{e_n} is ignored and assumed to be \code{Inf} or \code{NA}, indicating that the top bin is unbounded. The edges determine \eqn{n} bins on the intervals \eqn{e_{i-1} \le x \le e_i}, where \eqn{e_0} is assumed to be 0.
}
  \item{bCounts}{
A vector \eqn{c_1, c_2, \ldots, c_n} giving the counts for each bin (i.e., the number of data elements in each bin). Assumed to be nonnegative.}
  \item{m}{
An estimate for the mean of the distribution. If no value is supplied, the mean will be estimated by (temporarily) setting \eqn{e_n} equal to \eqn{2e_{n-1}}, and a warning message will be generated.
}
  \item{numIterations}{
The number of iterations performed by a binary search that optimizes the CDF to fit the mean.
}
  \item{monoMethod}{
The method for constructing a monotone spline. Must be one of \code{"hyman"} or \code{"monoH.FC"}. The former choice tends to integrate faster and produce smoother density functions. See \code{\link[stats]{splinefun}} for more details.
}
}
\details{
Fits a monotone cubic spline to the points specified by the binned data to produce a smooth cumulative distribution function. The PDF is then obtained by differentiating, so it will be piecewise quadratic and preserve the area of each bin.
}
\value{
Returns a list with the following components.
\item{splinePDF}{A piecewise-quadratic function giving the fitted PDF.}
\item{splineCDF}{A piecewise-cubic function giving the CDF.}
\item{E}{The right-hand endpoint of the support of the PDF.}
\item{shrinkFactor}{If the supplied estimate for the mean is too small to be fitted with our method, the bins edges will be scaled by \code{shrinkFactor}, which will be chosen less than (and close to) 1.}
\item{splineInvCDF}{An approximate inverse of \code{splineCDF}.}
\item{fitWarn}{Flag set to \code{TRUE} if the fitted median falls in the wrong bin.}
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
sb <- stepbins(binedges, bincounts, 76091)
splb <- splinebins(binedges, bincounts, 76091)

plot(splb$splinePDF, 0, 300000, n=500)
plot(sb$stepPDF, do.points=FALSE, col="gray", add=TRUE)
# notice that the curve preserves bin area

library(pracma)
integral(splb$splinePDF, 0, splb$E)
integral(function(x){1-splb$splineCDF(x)}, 0, splb$E) # should be the mean
splb <- splinebins(binedges, bincounts, 76091, numIterations=20)
integral(function(x){1-splb$splineCDF(x)}, 0, splb$E) # closer to given mean
}

