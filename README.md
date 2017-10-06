## rminuit2, numerical multidimensional function minimization with Minuit2

R package with an interface for using the
[Minuit2](https://root.cern.ch/root/html/MATH_MINUIT2_Index.html)
minimization library, which is included in this package.

Minuit2 is a robust numerical algorithm focused on
multidimensional minimization for High Energy Physics (HEP) fitting and maximum
likelihood estimation. It effectively finds the local minimum of
multi-parameter functions, without requiring any derivative of the
function. Under the assumption that the minimized
function is the negative logarithm of a likelihood (MLL) of measurements of
a Physics model, Minuit2 computes estimates of the optimized parameters'
uncertainties and covariance matrix by numerically computing the
Hesse matrix of the second derivatives of the MLL.
Minuit2 can also compute the  optimized parameters' asymmetric
uncertainties by scanning the function minimum to find how much a
parameter must change to increase the MLL by 1/2.

Minuit2 is a C++ implementation of the Fortran program Minuit, written
by Frank James. It is integrated with the HEP package
[Root](https://root.cern.ch/).

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
