#' Function Minimization with Minuit2
#'
#' Performs function minimization using
#' \href{https://project-mathlibs.web.cern.ch/project-mathlibs/sw/Minuit2/html/index.html}{Minuit2}
#'
#' @param fn The function to be minimized, with first argument a numeric
#'   vector of the parameters to be optimized. The function should be either an R function
#'   or an external pointer to a C++ function compiled using the
#'   \code{inline} interface. For more info about implementing the objective
#'   function as compiled C++ code, see the lbfgs package vignette
#' @param par numeric vector, initial values of the function parameters
#' @param err numeric vector, expected uncertainty of function parameters
#' @param lower numeric vector, lower bounds for the parameters (default none)
#' @param upper numeric vector, upper bounds for the parameters (default none)
#' @param fix boolean vector, each TRUE element fixes the corresponding parameter
#' @param opt string, pass fit options, default "h" (compute HESSE errors)
#'   \describe{
#'   \item{\code{v}:}{Verbose mode (not yet implemented)}
#'   \item{\code{h}:}{Run Hesse to estimate errors}
#'   \item{\code{m}:}{Get Minos errors}
#'   \item{\code{0}:}{Run Migrad with strategy 0}
#'   \item{\code{1}:}{Run Migrad with strategy 1}
#'   \item{\code{2}:}{Run Migrad with strategy 2}
#'   \item{neither \code{1}, \code{2}, \code{3}:}{Run Migrad with strategy 1, and if it fails run with strategy 2}
#'   }
#' @param envir An R environment containing all extra arguments to be passed
#'   to the objective function, which must be matched exactly. If the objective
#'   function is implemented in C++, extra arguments must be passed using this option,
#'   rather than the \code{...} construct. If the functions are implemented in R, extra
#'   arguments should be passed to them using the \code{...} construct instead.
#' @param ... construct. If the functions are implemented in R, extra arguments
#'   should be passed to them using the \code{...} construct instead.
#' @param maxcalls integer, maximum number of calls, defaults to \code{0} (no limit).
#' @param nsigma numeric, number of standard deviations for Minos errors
#'
#' @return A list with the following components:
#'  \describe{
#'    \item{\code{fval}:}{Value of function at found minimum (1/2 * chi square if the function is the negative log-likelihood of a Gaussian  likelihood.}
#'    \item{\code{Edm}:}{Estimated distance from the value of the function true minimum.}
#'    \item{\code{par}:}{Fitted parameters.}
#'    \item{\code{err}:}{Estimated uncertainties of fitted parameters.}
#'    \item{\code{cov}:}{Covariance matrix of the fitted parameters.}
#'    \item{\code{err_minos_pos}:}{Minos-estimated positive parameters' uncertainties (if Minos errors were requested).}
#'    \item{\code{err_minos_neg}:}{Minos-estimated negative parameters' uncertainties (if Minos errors requested).}
#'    \item{\code{err_minos_pos_valid}:}{boolean vector, TRUE if Minos positive uncertainties are valid (if Minos errors were requested).}
#'    \item{\code{err_minos_neg_valid}:}{boolean vector, TRUE if Minos negative uncertainties are valid (if Minos errors were requested).}
#'    \item{\code{allOK}:}{TRUE if the fit converged and the parameters and their covariance are OK}
#'    \item{\code{MinosErrorsValid}:}{TRUE if the MINOS errors are all valid}
#'    \item{\code{IsValid}:}{TRUE if the fit minimization converged}
#'    \item{\code{IsValidFirstInvocation}:}{TRUE if Minuit strategy 1 succeeded (if it failed Minuit2 strategy 2 is performed).}
#'    \item{\code{IsAboveMaxEdm}:}{TRUE if the estimated distance from the true minimum is above the tolerance.}
#'    \item{\code{HasReachedCallLimit}:}{TRUE if the maximum call limit was exceeded.}
#'    \item{\code{HasValidParameters}:}{TRUE if the fitted parameters are considered valid.}
#'    \item{\code{HasCovariance}:}{TRUE if a covariance matrix is returned.}
#'    \item{\code{HasValidCovariance}:}{TRUE if the estimated covariance matrix is considered valid.}
#'    \item{\code{HasAccurateCovar}:}{TRUE if the accuracy of the estimated covariance matrix is considered valid.}
#'    \item{\code{HasPosDefCovar}:}{TRUE if the numerically computed covariance matrix is positive definite.}
#'    \item{\code{HasMadePosDefCovar}:}{TRUE if the covariance matrix has been adjusted to make it positive definite.}
#'    \item{\code{HesseFailed}:}{TRUE if the numeric computation of the HESSE matrix failed.}
#'  }
#'
#' @author Alberto Lusiani, \email{alusiani@gmail.com}
#'
#' @keywords minimization fitting optimization
#'
#' @examples
#' #
#' # Rosenbrock Banana function
#' #
#' rosenbrock <- function(x) {
#'     x1 <- x[1]
#'     x2 <- x[2]
#'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#' }
#' 
#' # minimize Rosenbrock Banana function fitting its two parameters
#' fit.rc <- rminuit2(rosenbrock, c(a=0.7, b=1.2))
#'
#' # print fitted parameters
#' fit.rc$par
#' 
#' #
#' # simulate model y = a*exp(-x/b)
#' #
#' x = seq(0, 1, length.out=31)
#' y.func = function(x, par) par[1]*exp(-x/par[2])
#' 
#' # simulate data with Gaussian errors for specific model
#' model.par = c(a=2.3, b=0.47)
#' y.err = 0.01
#' y = y.func(x, par=model.par) + rnorm(sd=y.err, n=length(x))
#' 
#' # negative log-likelihood for model
#' halfchisq = function(par, x, y, y.err) {
#'   sum( (y - y.func(x, par))^2 / (2 * y.err^2) )
#' }
#' 
#' # fit model on data, ask to compute Minos errors too
#' fit.rc = rminuit2(halfchisq, c(a=-1, b=10), opt="hm", x=x, y=y, y.err=y.err)
#' 
#' # chi square / number of degrees of freedom
#' cbind(chisq=2*fit.rc$fval, ndof=length(x) - length(model.par))
#' 
#' # fitted parameters and their estimated uncertainties
#' cbind(value=fit.rc$par, error=fit.rc$err, minos_pos=fit.rc$err_minos_pos, minos_neg=fit.rc$err_minos_neg)
#'
#' # parameters' correlation matrix
#' cov2cor(fit.rc$cov)
#'
#' @useDynLib rminuit2
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom methods hasArg
#' @importFrom stats setNames
#'
#' @export
#'
rminuit2 <- function(fn, par, err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
  envir=NULL, ..., maxcalls=0L, nsigma=1) {

  ##--- Initialize environment if NULL
  ## if (!hasArg(envir)) envir <- new.env()
  if (is.null(envir)) envir <- new.env()

  ##--- Set up parameters
  npar <- length(par)
  xtol <- .Machine$double.eps

  if (is.null(err))   err = rep(0.1, npar)
  if (is.null(lower)) lower = rep(-Inf, npar)
  if (is.null(upper)) upper = rep(Inf, npar)
  if (is.null(fix))   fix = rep(0L, npar)
  if (is.null(opt))   opt = ""

  opt = tolower(as.character(opt))
  if (length(opt) < 1) {
    opt = ""
  } else if (length(opt) > 1) {
    opt = opt[1]
  }

  if (length(maxcalls) < 1) {
    maxcalls = 0L
  } else if (length(maxcalls) > 1) {
    maxcalls = maxcalls[1]
  }

  if (length(nsigma) < 1) {
    maxcalls = 1
  } else if (length(maxcalls) > 1) {
    nsigma = nsigma[1]
  }

  if (is.null(names(par))) names(par) = paste0("p", seq(1, npar))
  par = setNames(as.numeric(par), names(par))

  if (length(err) != npar) stop("vector 'err' must have the same length as 'par'")
  if (length(lower) != npar) stop("vector 'lower' must have the same length as 'par'")
  if (length(upper) != npar) stop("vector 'upper' must have the same length as 'par'")
  if (length(fix) != npar) stop("vector 'fix' must have the same length as 'par'")

  ##
  ## prepare args for c++ call preventing invalid input
  ##
  err = as.numeric(err)
  err = ifelse(is.na(err), 0.1, err)
  lower = as.numeric(lower)
  lower = ifelse(is.na(lower), -Inf, lower)
  upper = as.numeric(upper)
  upper = ifelse(is.na(upper), Inf, upper)
  fix = as.integer(fix)
  fix = ifelse(is.na(fix), 0L, fix)
  envir = as.environment(envir)
  maxcalls = as.integer(maxcalls)
  maxcalls = ifelse(is.na(maxcalls), 0L, maxcalls)
  nsigma = as.numeric(nsigma)
  nsigma = ifelse(is.na(nsigma), 1, abs(nsigma))

  ##--- fix errors to zero for fixed parameters
  err[which(fix!=0)] = 0

  ##--- Call main C++ routine
  rc <- .Call('_rminuit2_rminuit2_cpp', PACKAGE = 'rminuit2',
              fn, par, err, lower, upper,
              fix, opt, envir, maxcalls, nsigma)

  ##--- add names for output
  names(rc$par) = names(par)
  names(rc$err) = names(par)

  ##--- add names of Minos errors for output
  if (!is.null(rc$err_minos_pos) && !is.null(rc$err_minos_neg)) {
    names(rc$err_minos_pos) = names(par)
    names(rc$err_minos_neg) = names(par)
    rc$MinosErrorsValid = all(rc$err_minos_pos_valid, rc$err_minos_ned_valid)
    if (!rc$MinosErrorsValid) warning("One or more Minos errors are not valid")
  }

  ##
  ## returned covariance is restricted to just the non-fixed parameters
  ## add zero elements for all fixed parameters
  ##
  npar_fixed = npar - nrow(rc$cov)
  rc$cov = cbind(rc$cov, matrix(0, nrow(rc$cov), npar_fixed))
  rc$cov = rbind(rc$cov, matrix(0, npar_fixed, npar))
  cov.reorder = c( which(fix==0), which(fix!=0) )
  rc$cov[cov.reorder, cov.reorder] = rc$cov

  colnames(rc$cov) = names(par)
  rownames(rc$cov) = names(par)

  ##--- compute overall fit and covariance valid flag
  rc$allOK = all(
    rc$IsValid,
    rc$HasValidParameters,
    rc$HasValidCovariance,
    rc$HasAccurateCovar,
    rc$HasPosDefCovar,
    !rc$HesseFailed,
    rc$HasCovariance,
    !rc$IsAboveMaxEdm,
    !rc$HasReachedCallLimit
  )

  if (!rc$IsValid) warning("migrad failed")
  if (!rc$IsValidFirstInvocation) warning("migrad first invocation failed, use strategy 2")
  if (!rc$HasValidParameters) warning("parameters are not valid")
  if (!rc$HasValidCovariance) warning("covariance matrix is not valid")
  if (!rc$HasAccurateCovar) warning("covariance matrix is not accurate")
  if (!rc$HasPosDefCovar) warning("covariance matrix non positive definite")
  if (rc$HasMadePosDefCovar) warning("covariance matrix was adjusted to be positive definite")
  if (rc$HesseFailed) warning("Hesse procedure to get covariance numerically failed")
  if (!rc$HasCovariance) warning("no covariance matrix computed")
  if (rc$IsAboveMaxEdm) warning("IsAboveMaxEdm")
  if (rc$HasReachedCallLimit) warning("maximum number of calls limit reached")

  return(rc)
}
