## rminuit2, numerical multidimensional function minimization with Minuit2

R package with an interface for using the
[Minuit2](https://root.cern.ch/root/html/MATH_MINUIT2_Index.html)
minimization library, which is included in this package.

Minuit2 is a library to perform robust numerical multidimensional
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

### Non-linear least-square fitting with maximum-likelihood fit

This package also includes the ability to fit the parameters of a
model that aims to describe a set of experimental measurements.  The
model is described with a formula expression similarly to the
[nls](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/nls.html)
function for non-linear least squares fitting. The measurements' error
distribution is modelled as Gaussian and a maximum-likelihood fit is
performed. Errors that depend on the fit parameters can be included
in the model.

### Installation and Usage

Download the package tarball and build using R commands, or alternatively install directly from Github using Hadley Wickham's [devtools](https://github.com/hadley/devtools) package. The R commands are:

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

### Disclaimer

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in
the software.
