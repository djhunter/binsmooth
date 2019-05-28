# binsmooth
An R package for smoothing binned data, available on CRAN.

This package provides several methods for generating density functions based on binned data. 
Data are assumed to be nonnegative, but the bin widths need not be uniform, and the top bin may be unbounded. 
All PDF smoothing methods
maintain the areas specified by the binned data. (Equivalently, all CDF
smoothing methods interpolate the points specified by the binned data.) An
estimate for the mean of the distribution may be supplied as an optional
argument, which greatly improves the reliability of statistics computed from
the smoothed density functions. Methods include step function, recursive
subdivision, and optimized spline.

These methods are described in the paper 
[Better Estimates from Binned Income Data: Interpolated CDFs and Mean-Matching](https://www.sociologicalscience.com/articles-v4-26-641/),
by Paul T. von Hippel, David J. Hunter, McKalie Drown. *Sociological Science*, November 15, 2017, DOI 10.15195/v4.a26.

You can install the latest CRAN version through R:

```
install.packages("binsmooth")
```
