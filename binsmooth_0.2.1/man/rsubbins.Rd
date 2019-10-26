\name{rsubbins}
\alias{rsubbins}
\alias{rsubbinsTail}
\alias{rsubbinsNotail}
\title{
Recursive subdivision PDF and CDF fitted to binned data
}
\description{
Creates a PDF and CDF based on a set of binned data, using recursive subdivision on a step function.
}
\usage{
rsubbins(bEdges, bCounts, m=NULL, eps1 = 0.25, eps2 = 0.75, depth = 3,
        tailShape = c("onebin", "pareto", "exponential"),
        nTail=16, numIterations=20, pIndex=1.160964, tbRatio=0.8)
}
\arguments{
  \item{bEdges}{A vector \eqn{e_1, e_2, \ldots, e_n} giving the right endpoints of each bin. The value in \eqn{e_n} is ignored and assumed to be \code{Inf} or \code{NA}, indicating that the top bin is unbounded. The edges determine \eqn{n} bins on the intervals \eqn{e_{i-1} \le x \le e_i}, where \eqn{e_0} is assumed to be 0.
}
  \item{bCounts}{
A vector \eqn{c_1, c_2, \ldots, c_n} giving the counts for each bin (i.e., the number of data elements in each bin). Assumed to be nonnegative.
}
  \item{m}{
An estimate for the mean of the distribution. If no value is supplied, the mean will be estimated by (temporarily) setting \eqn{e_n} equal to \eqn{2e_{n-1}}, and a warning message will be generated.
}
\item{eps1}{
Parameter controlling how far the edges of the subdivided bins are shifted. Must be between 0 and 0.5.
}
\item{eps2}{
Parameter controlling how wide the middle subdivsion of each bin should be. Must be between 0 and 1.
}
\item{depth}{
Number of times to subdivide the bins.
}
  \item{tailShape}{
Must be one of \code{"onebin"}, \code{"pareto"}, or \code{"exponential"}.
}
  \item{nTail}{
The number of bins to use to form the initial tail, before recursive subdivision.
Ignored if \code{tailShape} equals \code{"onebin"}.
}
  \item{numIterations}{
The number of iterations to optimize the tail to fit the mean. Ignored if
\code{tailShape} equals \code{"onebin"}.
}
  \item{pIndex}{
The Pareto index for the shape of the tail. Defaults to \eqn{\ln(5)/\ln(4)}{ln(5)/ln(4)}.
Ignored unless \code{tailShape} equals \code{"pareto"}.
}
  \item{tbRatio}{
The decay ratio for the tail bins. Ignored unless \code{tailShape} equals \code{"exponential"}.
}
}
\details{
First, a step function PDF is created, as described in \code{\link{stepbins}}. The bins of the resulting PDF are then recursively subdivided and shifted in a manner that preserves the area of the original bins, resulting in a step function with finer bins.

The methods \code{\link{stepbins}} and \code{\link{rsubbins}} are included in this package mainly for the purpose of comparison. For most use cases, \code{\link{splinebins}} will produce more accurate smoothing results.
}
\value{
Returns a list with the following components.
\item{rsubPDF}{A \code{\link[stats]{stepfun}} function giving the fitted PDF.}
\item{rsubCDF}{A piecewise-linear \code{\link[stats]{approxfun}} function giving the CDF.}
\item{E}{The right-hand endpoint of the support of the PDF.}
\item{shrinkFactor}{If the supplied estimate for the mean is too small to be fitted with a step function, the bins edges will be scaled by \code{shrinkFactor}, which will be chosen less than (and close to) 1.}
}
\references{
Paul T. von Hippel, David J. Hunter, McKalie Drown. \emph{Better Estimates from Binned Income Data: Interpolated CDFs and Mean-Matching}, Sociological Science, November 15, 2017. \url{https://www.sociologicalscience.com/articles-v4-26-641/}
}
\author{
David J. Hunter and McKalie Drown
}
\seealso{\code{\link{stepbins}}}
\examples{
# 2005 ACS data from Cook County, Illinois
binedges <- c(10000,15000,20000,25000,30000,35000,40000,45000,
              50000,60000,75000,100000,125000,150000,200000,NA)
bincounts <- c(157532,97369,102673,100888,90835,94191,87688,90481,
               79816,153581,195430,240948,155139,94527,92166,103217)
rsb <- rsubbins(binedges, bincounts, 76091, tailShape="pareto")

plot(rsb$rsubPDF, do.points=FALSE)
plot(rsb$rsubCDF, 0, rsb$E)

library(pracma)
integral(rsb$rsubPDF, 0, rsb$E)
integral(function(x){1-rsb$rsubCDF(x)}, 0, rsb$E) #mean is approximated
}

