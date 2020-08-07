# binsmooth
An R package for smoothing binned data, available on CRAN.

This package provides several methods for generating density functions based on binned data. Methods include step function, recursive
subdivision, and optimized spline. Data are assumed to be nonnegative, 
the top bin is assumed to have no upper bound, but the bin widths need not
be equal. All PDF smoothing methods maintain the areas specified by 
the binned data. (Equivalently, all CDF smoothing methods interpolate 
the points specified by the binned data.) In practice, an estimate for 
the mean of the distribution should be supplied as an optional argument.
Doing so greatly improves the reliability of statistics computed from 
the smoothed density functions. Includes methods for estimating the Gini 
coefficient, the Theil index, percentiles, and random deviates from a 
smoothed distribution. Among the three methods, the optimized spline 
(splinebins) is recommended for most purposes. The percentile and 
random-draw methods should be regarded as experimental, and these methods 
only support splinebins.

These methods are described in the paper 
[Better Estimates from Binned Income Data: Interpolated CDFs and Mean-Matching](https://www.sociologicalscience.com/articles-v4-26-641/),
by Paul T. von Hippel, David J. Hunter, McKalie Drown. *Sociological Science*, November 15, 2017, DOI 10.15195/v4.a26.

You can install the latest CRAN version through R:

```
install.packages("binsmooth")
```
