## rminuit2, numerical multidimensional function minimization with Minuit2

R package with an interface for using the
[Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html)
minimization library, which is included in this package.

Minuit2 is a library to perform robust numerical multidimensional
minimization. It effectively finds the local minimum of
a multi-parameter function, without requiring its gradient.

Minuit2 is a C++ implementation of the Fortran program Minuit, written
by Frank James. It is integrated with the High Energy Physics (HEP) package
[Root](https://root.cern.ch/) [see Rene Brun and Fons Rademakers,
ROOT - An Object Oriented Data Analysis Framework,
Proceedings AIHENP'96 Workshop, Lausanne, Sep. 1996,
Nucl. Inst. & Meth. in Phys. Res. A 389 (1997) 81-86].
Minuit and Minuit2 have been developed
and are mainly used for HEP data analysis.

In Physics, measurements are supposed to be probabilistically
distributed according to the predictions of a Physics model and the
resolution of the instrumentation used for the
measurements. Typically, the Physics model includes unknown parameters
to be determined via measurements. The likelihood of a set of
measurements is the product of the probability density of each
individual measurement. The optimal Physics model parameters are the
ones that maximize the likelihood or, equivalently, minimize the
negative logarithm of the likelihood (Minus Log-Likelihood, MLL).

When Minuit2 minimizes a MLL function, it determines the optimal
Physics model parameters and it also estimates the
optimized parameters' uncertainties and covariance matrix by
numerically computing the Hesse matrix of the second derivatives of
the MLL at the minimum. Minuit2 can also compute the optimized
parameters' asymmetric uncertainties by scanning the function minimum
to find how much a parameter must change to increase the MLL by 1/2.

With respect to other minimization libraries, Minuit2 is very
effective in minimizing functions with many parameters with large
correlations and different scales. Furthermore, it devotes special
care in numerically estimating the parameter's uncertainties.

### Minimization of a C++ function

Relying on the packages [Rcpp](https://github.com/RcppCore/Rcpp) and [inline](https://CRAN.R-project.org/package=inline) and following the examples n the packages [lbfgs](https://CRAN.R-project.org/package=lbfgs) and [RcppDE](http://cran.r-project.org/web/packages/RcppDE/index.html), it is possible to define the function to be minimized in C++.

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

Download the package tarball and build using R commands, or alternatively install directly from Github using the [remotes](https://CRAN.R-project.org/package=remotes) package. The R commands are:

```
## install.packages("remotes")
remotes::install_github("alusiani/rminuit2")
```

The package has been tested so far only on Linux, and is not prepared
for MS Windows. Adaptation to MS Windows is not trivial since Minuit2
is presently compiled with autotool/automake Makefiles.
Contributors are welcome for adaptations for OSs other than Linux.

For usage, please refer to the documentation.

### Dependencies

- [Rcpp](https://github.com/RcppCore/Rcpp) for seamless R and C++ integration

### Acknowledgements

This package relies on the infrastracture for interfacing C++ offered
by [Rcpp](https://github.com/RcppCore/Rcpp).

To offer the possibility to write in C++ the function to be
minimized, we used code and examples from the packages:

* Antonio Coppola, Brandon Stewart, Naoaki Okazaki, [lbfgs: Efficient L-BFGS and OWL-QN Optimization in R](https://CRAN.R-project.org/package=lbfgs);
* Dirk Eddelbuettel, [RcppDE: Global optimization by Differential Evolution in C++](http://cran.r-project.org/web/packages/RcppDE/index.html).

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
