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
#'   \itemize{
#'   \item{v}{Verbose mode (not yet implemented)}
#'   \item{h}{Run Hesse to estimate errors}
#'   \item{m}{Get Minos errors}
#'   \item{0}{Run Migrad with strategy 0}
#'   \item{1}{Run Migrad with strategy 1}
#'   \item{2}{Run Migrad with strategy 2}
#'   \item{no 123}{Run Migrad with strategy 1 and if fails 2}
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
#'  \itemize{
#'    \item{value}{The minimized value of the objective function.}
#'    \item{par}{A numerical array with the parameters that minimize the function.}
#'    \item{par_err}{A numerical array with the estimeted parameters' uncertainties.}
#'    \item{par_cov}{A matrix with the estimated covariance matrix of the parameters.}
#'    \item{convergence}{An integer code. Zero indicates that convergence was reached without issues. Negative values indicate errors in the execution of the L-BFGS routine.}
#'    \item{message}{A character object detailing execution errors. This component is only returned if the convergence code is different form zero.}
#'  }
#'
#' @author Alberto Lusiani, \email{alusiani@gmail.com}
#'
#' @keywords minimization fitting optimization
#'
#' @examples
#' # Rosenbrock Banana function
#'
#' rosenbrock <- function(x) {
#'     x1 <- x[1]
#'     x2 <- x[2]
#'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#' }
#'
#' output <- rminuit2(rosenbrock, c(-1.2, 1))
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
  nsigma = ifelse(is.na(nsigma), 0, abs(nsigma))

  ##--- fix errors to zero for fixed parameters
  err[which(fix!=0)] = 0

  ##--- Call main C++ routine
  rc <- .Call('_rminuit2_rminuit2_cpp', PACKAGE = 'rminuit2',
              fn, par, err, lower, upper,
              fix, opt, envir, maxcalls, nsigma)

  names(rc$par) = names(par)
  names(rc$err) = names(par)

  if (!is.null(rc$err_minos_pos)) {
    names(rc$err_minos_pos) = names(par)
  }
  if (!is.null(rc$err_minos_neg)) {
    names(rc$err_minos_neg) = names(par)
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
