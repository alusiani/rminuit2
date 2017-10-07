## rminuit2, numerical multidimensional function minimization with Minuit2

R package with an interface for using the
[Minuit2](https://root.cern.ch/root/html/MATH_MINUIT2_Index.html)
minimization library, which is included in this package.

Minuit2 is a library to perform robust numerical miltidimensional
minimization. It effectively finds the local minimum of
a multi-parameter function, without requiring its gradient.

Minuit2 is a C++ implementation of the Fortran program Minuit, written
by Frank James. It is integrated with the HEP package
[Root](https://root.cern.ch/). Minuit and Minuit2 have been developed
and are mainly used for High Energy Physics (HEP) data analysis.

Assuming that the minimized function is the negative logarithm of a
likelihood (MLL) of measurements of a Physics model, depending on
parameters to be fitted to the measurements, Minuit2 estimates the
optimized parameters' uncertainties and covariance matrix by
numerically computing the Hesse matrix of the second derivatives of
the MLL at the minimum. Minuit2 can also compute the optimized
parameters' asymmetric uncertainties by scanning the function minimum
to find how much a parameter must change to increase the MLL by 1/2.

With respect to other minimization libraries, Minuit2 is very
effective in minimizing functions with many parameters with large
correlations and different scales. Furthermore, it devotes special
care in numerically estimating the parameter's uncertainties.

### Installation and Usage

Download the package tarball and build using R commands, or alternatively instally directly from Github using Hadley Wickham's [devtools](https://github.com/hadley/devtools) package. The R commands are:

```
library(devtools)
install_github("alusiani/rminuit2")
```

The package has been tested so far only on Linux, and is not prepared
for MS Windows. Adaptation to MS Windows is not trivial since Minuit2
is presently compiled with autotool/automake Makefiles.
Contributors are welcome for adaptations for OSs other than Linux.

For usage, please refer to the documentation.

### Dependencies

- [Rcpp](https://github.com/RcppCore/Rcpp) for seamless R and C++ integration

### Authors

Alberto Lusiani

### License

LGPL (>= 2)
