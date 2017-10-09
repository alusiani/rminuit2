#' Function Minimization with Minuit2
#'
#' Performs function minimization using
#' \href{https://project-mathlibs.web.cern.ch/project-mathlibs/sw/Minuit2/html/index.html}{Minuit2}
#'
#' @param mll The function to be minimized.
#'   The full potential of this package is attained when the function corresponds to
#'   the negative logarithm of a likelihood or minus-log-likelihood (MLL).
#' @param start numeric vector, initial values of the function parameters
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
#' @param maxcalls integer, maximum number of calls, defaults to \code{0} (no limit).
#' @param nsigma numeric, number of standard deviations for Minos errors
#' @param envir not used. Will possibly be used in the future as the same argument in \code{rminuit2_par}.
#' @param ... extra arguments for the mll function.
#'
#' @return A list with the following components:
#'  \describe{
#'    \item{\code{fval}:}{Value of function at found minimum (1/2 * chi square if the function is the negative log-likelihood of a Gaussian  likelihood.}
#'    \item{\code{Edm}:}{Estimated distance from the value of the function true minimum}
#'    \item{\code{par}:}{Fitted parameters}
#'    \item{\code{err}:}{Estimated uncertainties of fitted parameters}
#'    \item{\code{cov}:}{Covariance matrix of the fitted parameters}
#'    \item{\code{err_minos_pos}:}{Minos-estimated positive parameters' uncertainties (if Minos errors were requested)}
#'    \item{\code{err_minos_neg}:}{Minos-estimated negative parameters' uncertainties (if Minos errors requested)}
#'    \item{\code{err_minos_pos_valid}:}{boolean vector, TRUE if Minos positive uncertainties are valid (if Minos errors were requested)}
#'    \item{\code{err_minos_neg_valid}:}{boolean vector, TRUE if Minos negative uncertainties are valid (if Minos errors were requested)}
#'    \item{\code{allOK}:}{TRUE if the fit converged and the parameters and their covariance are OK}
#'    \item{\code{MinosErrorsValid}:}{TRUE if the MINOS errors are all valid}
#'    \item{\code{IsValid}:}{TRUE if the fit minimization converged}
#'    \item{\code{IsValidFirstInvocation}:}{TRUE if Minuit strategy 1 succeeded (if it failed Minuit2 strategy 2 is performed)}
#'    \item{\code{IsAboveMaxEdm}:}{TRUE if the estimated distance from the true minimum is above the tolerance}
#'    \item{\code{HasReachedCallLimit}:}{TRUE if the maximum call limit was exceeded}
#'    \item{\code{HasValidParameters}:}{TRUE if the fitted parameters are considered valid}
#'    \item{\code{HasCovariance}:}{TRUE if a covariance matrix is returned}
#'    \item{\code{HasValidCovariance}:}{TRUE if the estimated covariance matrix is considered valid}
#'    \item{\code{HasAccurateCovar}:}{TRUE if the accuracy of the estimated covariance matrix is considered valid}
#'    \item{\code{HasPosDefCovar}:}{TRUE if the numerically computed covariance matrix is positive definite}
#'    \item{\code{HasMadePosDefCovar}:}{TRUE if the covariance matrix has been adjusted to make it positive definite}
#'    \item{\code{HesseFailed}:}{TRUE if the numeric computation of the HESSE matrix failed}
#'  }
#'
#' @examples
#'
#' #
#' # Rosenbrock Banana function, to be minimized vs. x, y
#' #
#' rosenbrock <- function(x, y, a, b) {
#'   (a-x)^2 + b*(y-x^2)^2
#' }
#'
#' # minimize Rosenbrock Banana function, also setting parameters a, b
#' fit.rc <- rminuit2(rosenbrock, c(x=0.7, y=1.2), a=1, b=100)
#' fit.rc$par
#'
#' #
#' # simulate model y = a*exp(-x/b)
#' #
#' x = seq(0, 1, length.out=31)
#' y.func = function(x, norm, tau) norm*exp(-x/tau)
#'
#' # simulate data with Gaussian errors for specific model
#' model.par = c(norm=2.3, tau=0.47)
#' y.err = 0.01
#' y = do.call(y.func, c(list(x), model.par)) + rnorm(sd=y.err, n=length(x))
#'
#' # negative log-likelihood for model
#' halfchisq = function(norm, tau, x, y, y.err) {
#'   sum( (y - y.func(x, norm, tau))^2 / (2 * y.err^2) )
#' }
#'
#' # fit model on data, ask to compute Minos errors too
#' fit.rc = rminuit2(halfchisq, c(norm=1, tau=10), opt="hm", x=x, y=y, y.err=y.err)
#'
#' # chi square / number of degrees of freedom
#' cbind(chisq=2*fit.rc$fval, ndof=length(x) - length(model.par))
#'
#' # fitted parameters and their estimated uncertainties
#' cbind(model.value=model.par, value=fit.rc$par, error=fit.rc$err,
#'       minos_pos=fit.rc$err_minos_pos, minos_neg=fit.rc$err_minos_neg)
#'
#' # parameters' correlation matrix
#' cov2cor(fit.rc$cov)
#'
#' @author Alberto Lusiani, \email{alusiani@gmail.com}
#'
#' @keywords minimization fitting optimization
#'
#' @useDynLib rminuit2
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom methods hasArg
#' @importFrom stats setNames
#'
#' @seealso rminuit2_par
#' @export
#'
rminuit2 <- function(mll, start = formals(mll), err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
                     maxcalls=0L, nsigma=1, envir=NULL, ...) {
  call <- match.call()
  mll.args <- formals(mll)

  if (is.numeric(start)) start = as.list(start)
  if (!missing(start) && (!is.list(start) || is.null(names(start))))
    stop("'start' must be a named list or a named numeric vector")

  ##
  ## evaluate expressions in start
  ## catch eval.parent() errors to get meaningful error messages
  ##
  start <- tryCatch(start <- sapply(start, eval.parent),
    error = function(e) {
      for(el in start) {
        rc = tryCatch(eval(el), error=function(e) {NA})
        if (is.na(rc)) {
          stop(paste0("mll function argument '", el, " cannot be evaluated, add it in start or extra args"))
        }
      }
    })

  start.names <- names(start)
  mll.args.match <- match(start.names, names(mll.args))
  not.init = setdiff(names(mll.args), start.names)
  not.init = setdiff(not.init, names(list(...)))
  not.init.flag = sapply(mll.args[not.init], is.symbol)
  if (any(not.init.flag)) {
    stop(paste0("\n  mll function arguments: ",
                paste0("'", not.init[not.init.flag], "'", collapse=", "),
                "\n  not initialized in either the mll function, in start, or in optional parameters"))
  }
  if (anyNA(mll.args.match))
    stop("some named arguments in 'start' are not arguments to the supplied mll function")
  start <- start[order(mll.args.match)]
  start.names <- names(start)
  mll.par <- function(par, ...) {
    par.list <- as.list(par)
    names(par.list) <- start.names
    par.list = c(par.list, list(...))
    do.call("mll", par.list)
  }
  rminuit2_par(mll.par, start=start, err=err, lower=lower, upper=upper, fix=fix, opt=opt,
               maxcalls=maxcalls, nsigma=nsigma, envir=envir, ...)
}

#' Function Minimization with Minuit2
#'
#' Performs function minimization using
#' \href{https://project-mathlibs.web.cern.ch/project-mathlibs/sw/Minuit2/html/index.html}{Minuit2}
#'
#' @inherit rminuit2
#'
#' @param mll The function to be minimized, which must have as first
#'   argument a numeric vector of the parameters to be optimized. Futher
#'   arguments can be specified as optional arguments in \code{rminuit2_par}.
#'
#'   The function can
#'   also be an external pointer to a C++ function compiled using the
#'   \code{inline} interface. For more info about implementing the objective
#'   function as compiled C++ code, see the
#'   \href{https://cran.r-project.org/web/packages/lbfgs/index.html}{lbfgs} package vignette.
#'
#'   The full potential of this package is attained when the function corresponds to
#'   the negative logarithm of a likelihood or minus-log-likelihood (MLL).
#'
#' @param envir An R environment containing all extra arguments to be passed
#'   to the mll function, if the mll function is implemented in C++. Arguments
#'   must be matched exactly. If the mll function is implemented in R
#'   then the extra arguments should be passed to it using the optional arguments
#'   in \code{...} instead.
#'
#' @param ... extra arguments for the mll function, it the mll function is implemented in R.
#'   If the mll function is implemented in C++ then the extra arguments
#'   should be passed to it using the \code{envir} argument instead.
#'
#' @seealso rminuit2
#'
#' @examples
#' #
#' # Rosenbrock Banana function, to be minimized vs. 2 paramaters
#' #
#' rosenbrock <- function(par, a, b) {
#'   x <- par[1]
#'   y <- par[2]
#'   (a-x)^2 + b*(y-x^2)^2
#' }
#'
#' # minimize Rosenbrock Banana function, also setting parameters a, b
#' fit.rc <- rminuit2_par(rosenbrock, c(x=0.7, y=1.2), a=1, b=100)
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
#' fit.rc = rminuit2_par(halfchisq, c(a=1, b=10), opt="hm", x=x, y=y, y.err=y.err)
#'
#' # chi square / number of degrees of freedom
#' cbind(chisq=2*fit.rc$fval, ndof=length(x) - length(model.par))
#'
#' cbind(model.value=model.par, value=fit.rc$par, error=fit.rc$err,
#'       minos_pos=fit.rc$err_minos_pos, minos_neg=fit.rc$err_minos_neg)
#'
#' # parameters' correlation matrix
#' cov2cor(fit.rc$cov)
#'
#' @export
#'
rminuit2_par <- function(mll, start, err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
                         maxcalls=0L, nsigma=1, envir=NULL, ...) {

  ##--- Initialize environment if NULL
  ## if (!hasArg(envir)) envir <- new.env()
  if (is.null(envir)) envir <- new.env()

  ## xtol <- .Machine$double.eps

  npar <- length(start)
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

  if (is.null(names(start))) names(start) = paste0("p", seq(1, npar))
  par.names = names(start)
  start = setNames(as.numeric(start), par.names)

  if (!is.null(names(err))) {
    if (any(! names(err) %in% par.names))
      stop("some named arguments in 'err' are not parameters to be minimized in 'start'")
    err = setNames(ifelse(par.names %in% names(err), err, 0.1), par.names)
  }

  if (!is.null(names(lower))) {
    if (any(! names(lower) %in% par.names))
      stop("some named arguments in 'lower' are not parameters to be minimized in 'start'")
    lower = setNames(ifelse(par.names %in% names(lower), lower, -Inf), par.names)
  }

  if (!is.null(names(upper))) {
    if (any(! names(upper) %in% par.names))
      stop("some named arguments in 'upper' are not parameters to be minimized in 'start'")
    upper = setNames(ifelse(par.names %in% names(upper), upper, Inf), par.names)
  }

  if (!is.null(names(fix))) {
    if (any(! names(fix) %in% par.names))
      stop("some named arguments in 'fix' are not parameters to be minimized in 'start'")
    fix = setNames(ifelse(par.names %in% names(fix), fix, 0), par.names)
  }

  if (length(err) != npar) stop("if unnamed, vector 'err' must have the same length as 'start'")
  if (length(lower) != npar) stop("if unnamed, vector 'lower' must have the same length as 'start'")
  if (length(upper) != npar) stop("if unnamed, vector 'upper' must have the same length as 'start'")
  if (length(fix) != npar) stop("if unnamed, vector 'fix' must have the same length as 'start'")

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
              mll, start, err, lower, upper,
              fix, opt, envir, maxcalls, nsigma)

  ##--- add names for output
  names(rc$par) = par.names
  names(rc$err) = par.names

  ##--- add names of Minos errors for output
  if (!is.null(rc$err_minos_pos) && !is.null(rc$err_minos_neg)) {
    names(rc$err_minos_pos) = par.names
    names(rc$err_minos_neg) = par.names
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

  colnames(rc$cov) = par.names
  rownames(rc$cov) = par.names

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

##
## check if error return code
##
is.error <- function(x) inherits(x, "try-error")

##
## copied from package pryr
## convert argument to environment
##
to_env <- function(x, parent = parent.frame(), quiet = FALSE) {
  if (is.environment(x)) {
    x
  } else if (is.list(x)) {
    list2env(x)
  } else if (is.function(x)) {
    environment(x)
  } else if (length(x) == 1 && is.character(x)) {
    if (!quiet) message("Using environment ", x)
    as.environment(x)
  } else if (length(x) == 1 && is.numeric(x) && x > 0) {
    if (!quiet) message("Using environment ", search()[x])
    as.environment(x)
  } else {
    stop("Input can not be coerced to an environment", call. = FALSE)
  }
}

##
## copied from package pryr
## check if all elements are named with no empty strings
##
all_named <- function(x) {
  if (length(x) == 0) return(TRUE)
  !is.null(names(x)) && all(names(x) != "")
}

##
## copied from package pryr
## make a function
##
make_function <- function (args, body, env = parent.frame())
{
  args <- as.pairlist(args)
  stopifnot(all_named(args), is.language(body))
  env <- to_env(env)
  eval(call("function", args, body), env)
}

##
## given a model formula, assemble minus log-likelihood for gaussian errors
##
rminuit2_make_gaussian_mll <- function(formula, par, data=NULL, weights=NULL, errors=NULL) {
  weights = substitute(weights)
  errors = substitute(errors)

  if (length(formula) == 2) {
    residexpr <- formula[[2]]
  } else if (length(formula) == 3) {
    residexpr <- call("-", formula[[2]], formula[[3]])
  } else stop("Unrecognized formula")

  ##--- name parameters p<n> if unnamed
  if (is.null(names(par)))
    names(par) <- paste0("p", seq_along(par))

  if (is.null(data)) {
    ##--- if no data, get from parent frame
    data <- environment(formula)
  } else {
    data = withCallingHandlers(
      to_env(data, parent = environment(formula)),
      error = function(e) {
        e$message="'data' must be a dataframe, list or environment"
        e$call = sys.call(-2)
        stop(e)
      })
  }

  if (!is.null(errors)) fbody = as.call(c(as.name("/"), residexpr, errors))
  if (!is.null(weights)) fbody = as.call(c(as.name("*"), fbody, weights))
  fbody = as.call(c(as.name("^"), fbody, 2))
  fbody = as.call(c(quote(sum), fbody))
  fbody = as.call(c(as.name("*"), 1/2, fbody))
  fbody = as.call(c(quote(evalq), fbody, quote(localdata)))

  fbody = as.call(c(
    as.name("{"),
    list(quote(if (is.null(names(fpar))) names(fpar) <- names(par)),
         quote(localdata <- list2env(as.list(fpar), parent = data)),
         fbody)))

  make_function(alist(fpar=), fbody)
}

#' Function Minimization with Minuit2
#'
#' Fit data to a model performing minus log-likelihood minimization using
#' \href{https://project-mathlibs.web.cern.ch/project-mathlibs/sw/Minuit2/html/index.html}{Minuit2}
#' and assuming Gaussian uncertainties
#'
#' @inherit rminuit2
#'
#' @param formula expression describing the model to be fitted to data
#'
#' @param start initial values of the model parameters to be fitted
#'
#' @param data list or data.frame containing the data to be fitted to the model
#'
#' @param weights formula corresponding to weights to assign to the data observations
#'
#' @param errors formula corresponding to the uncertaintites of the data observations
#'
#' @param ... extra arguments for the model, weights and error formulas
#'
#' @seealso rminuit2 rminuit2_par
#'
#' @examples
#' #
#' # simulate model y = a*exp(-x/b)
#' #
#' x = seq(0, 1, length.out=31)
#' y.func = function(x, norm, tau) norm*exp(-x/tau)
#'
#' # simulate data with Gaussian errors for specific model
#' model.par = c(norm=2.3, tau=0.47)
#' y.err = 0.01
#' y = do.call(y.func, c(list(x), model.par)) + rnorm(sd=y.err, n=length(x))
#'
#' # fit model on data, ask to compute Minos errors too
#' fit.rc = rminuit2_expr_gaussian(y ~ norm*exp(-x/tau), c(norm=1, tau=10),
#'   data=data.frame(x=x, y=y, y.err=y.err), errors=y.err, opt="hm")
#'
#' # chi square / number of degrees of freedom
#' cbind(chisq=2*fit.rc$fval, ndof=length(x) - length(model.par))
#'
#' # fitted parameters and their estimated uncertainties
#' cbind(model.value=model.par, value=fit.rc$par, error=fit.rc$err,
#'       minos_pos=fit.rc$err_minos_pos, minos_neg=fit.rc$err_minos_neg)
#'
#' # parameters' correlation matrix
#' cov2cor(fit.rc$cov)
#'
#' @export
#'
rminuit2_expr_gaussian = function(formula, start, data=NULL, weights=NULL, errors=NULL,
                                  err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
                                  maxcalls=0L, nsigma=1, envir=NULL, ...) {
  mll = eval(substitute(
    rminuit2_make_gaussian_mll(formula=formula, par=start, data=data, weights=weights, errors=errors)))

  rminuit2_par(mll, start=start, err=err, lower=lower, upper=upper, fix=fix, opt=opt,
               maxcalls=maxcalls, nsigma=nsigma, envir=envir, ...)
}
