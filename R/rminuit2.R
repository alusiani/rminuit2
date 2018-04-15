#' Function Minimization and model fitting with Minuit2
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
#' @param ... extra arguments for the mll function.
#'
#' @return A list with the following components:
#'  \describe{
#'    \item{\code{fval}:}{Value of function at found minimum (1/2 * chi square if the function is the negative log-likelihood of a Gaussian  likelihood.}
#'    \item{\code{Edm}:}{Estimated distance from the value of the function true minimum}
#'    \item{\code{par}:}{Fitted parameters}
#'    \item{\code{err}:}{Estimated uncertainties of fitted parameters}
#'    \item{\code{cov}:}{Covariance matrix of the fitted parameters}
#'    \item{\code{cor}:}{Correlation matrix of the fitted parameters}
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
#' @seealso rminuit2_par
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
#' fit.rc$cor
#'
#' @author Alberto Lusiani, \email{alusiani@gmail.com}
#'
#' @keywords minimization fitting optimization
#'
#' @useDynLib rminuit2
#'
#' @importFrom methods hasArg
#' @importFrom stats setNames cov2cor
#' @importFrom Rcpp sourceCpp evalCpp
#'
#' @export
#'
rminuit2 <- function(mll, start = formals(mll), err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
                     maxcalls=0L, nsigma=1, ...) {
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
  not.init = setdiff(names(mll.args), start.names)
  not.init = setdiff(not.init, names(list(...)))
  not.init.flag = sapply(mll.args[not.init], is.symbol)
  if (any(not.init.flag)) {
    stop(paste0("\n  mll function arguments: ",
                paste0("'", not.init[not.init.flag], "'", collapse=", "),
                "\n  not initialized in either the mll function, in start, or in optional parameters"))
  }
  mll.args.match <- match(start.names, names(mll.args))
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
  rc = rminuit2_par(mll.par, start=start, err=err, lower=lower, upper=upper, fix=fix, opt=opt,
                    maxcalls=maxcalls, nsigma=nsigma, ...)
  invisible(rc)
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
#'   arguments can be specified as optional arguments in \code{rminuit2_par}
#'   or in the environment passed using the \code{envir} argument. 
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
#' @param ... extra arguments for the mll function.
#'   If the mll function is implemented in C++ then the extra arguments
#'   are collected in an environment, which is passed to the function.
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
#' fit.rc$cor
#'
#' @export
#'
rminuit2_par <- function(mll, start, err=NULL, lower=NULL, upper=NULL, fix=NULL, opt="h",
                         maxcalls=0L, nsigma=1, ...) {

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
    fix = setNames(ifelse(par.names %in% names(fix), fix[par.names], 0), par.names)
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
  maxcalls = as.integer(maxcalls)
  maxcalls = ifelse(is.na(maxcalls), 0L, maxcalls)
  nsigma = as.numeric(nsigma)
  nsigma = ifelse(is.na(nsigma), 1, abs(nsigma))

  ##--- fix errors to zero for fixed parameters
  err[which(fix!=0)] = 0

  ##--- when minimizing a C++ function, extra args must be in envir argument
  envir.list = list(...)
  envir = new.env()
  for (n in names(envir.list)) assign(n, envir.list[[n]], envir)

  ##--- Call main C++ routine
  rc <- .Call('_rminuit2_rminuit2_cpp', PACKAGE = 'rminuit2',
              mll, start, err, lower, upper,
              fix, opt, envir, maxcalls, nsigma)

  rm(envir)

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

  ##--- compute correlation matrix
  rc$cor = cov2cor(rc$cov)

  ##--- number of fixed parameter derived from the size of covariance returned by Minuit2
  npar_fixed = npar - nrow(rc$cov)

  ##
  ## returned covariance is restricted to just the non-fixed parameters
  ## add zero elements for all fixed parameters
  ##
  rc$cov = cbind(rc$cov, matrix(0, nrow(rc$cov), npar_fixed))
  rc$cov = rbind(rc$cov, matrix(0, npar_fixed, npar))
  cov.reorder = c( which(fix==0), which(fix!=0) )
  rc$cov[cov.reorder, cov.reorder] = rc$cov

  colnames(rc$cov) = par.names
  rownames(rc$cov) = par.names

  ##
  ## fix correlation matrix for fixed parameters
  ##
  rc$cor = cbind(rc$cor, matrix(0, nrow(rc$cor), npar_fixed))
  rc$cor = rbind(rc$cor, matrix(0, npar_fixed, npar))
  cor.reorder = c( which(fix==0), which(fix!=0) )
  rc$cor[cor.reorder, cor.reorder] = rc$cor
  diag(rc$cor) = 1
 
  colnames(rc$cor) = par.names
  rownames(rc$cor) = par.names

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

  invisible(rc)
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
## make a call
##
make_call <- function (f, ..., .args = list())
{
    if (is.character(f))
        f <- as.name(f)
    as.call(c(f, ..., .args))
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
rminuit2_make_gaussian_mll <- function(formula, par, data=NULL, weights=NULL, errors=NULL, rhs_vars=NULL, ...) {
  weights = substitute(weights)
  errors = substitute(errors)

  if (length(par) == 0) {
    stop("argument 'par' must contain at least one parameter to be fitted")
  }

  if (length(formula) == 2L) {
    residexpr <- formula[[2L]]
    fun_expr = NULL
  } else if (length(formula) == 3L) {
    residexpr <- call("-", formula[[2L]], formula[[3L]])
    fun_expr = formula[[3L]]
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
        e$call = sys.call(-3)
        stop(e)
      })
  }

  ##--- prepare fit parameters as text for args of a function, in a list
  par.txt = mapply(function(n, v) paste(n, "=", v), names(par), par, SIMPLIFY=FALSE)

  extra = list(...)
  if (any(names(extra) %in% names(par))) {
    stop(paste0("Some extra args have the same name as the parameters:\n  ",
               paste("'", names(extra)[names(extra) %in% names(par)], "'", sep=", ")))
  }
  ##--- prepare extra optional args as text for args of a function, in a list
  extra.txt = mapply(function(n, v) paste(n, "=", v), names(extra), extra, SIMPLIFY=FALSE)

  ##--- prepare list with all par and extra, with par grouped in a single par numeric vector
  par.par.extra.txt = c(paste0("par=c(",
                               paste(as.character(par.txt), collapse=", "),
                               ")"),
                        extra.txt)

  ##--- prepare list with all par and extra, with each par in one argument
  par.extra.txt = c(par.txt, extra.txt)
  
  ##
  ## assemble body of mll function using formula, weights and errors expressions
  ##
  fbody = residexpr
  if (!is.null(errors)) fbody = as.call(c(as.name("/"), fbody, errors))
  ##--- function corresponding to pulls
  fbody.pulls = fbody
  fbody = as.call(c(as.name("^"), fbody, 2))
  if (!is.null(weights)) fbody = as.call(c(as.name("*"), fbody, weights))
  fbody = as.call(c(as.name("*"), 1/2, fbody))
  ##--- include sigma term in likelihood (matters only also errors are fitted)
  fbody = as.call(c(as.name("+"), bquote(.(1/2*log(2*pi)) + log(.(errors))), fbody))
  fbody = as.call(c(quote(sum), fbody))

  fbody = as.call(c(
    as.name("{"),
    list(
      bquote(if (is.null(names(par))) names(par) <- .(names(par))),
      quote(mapply(function(name, val) assign(name, val, pos=parent.frame(2)), names(par), par)),
      fbody
    )))

  ##
  ## (memo) mll_fun = evalq(make_function(alist(par=), fbody), envir=data)
  ##
  ## could not find a way with make_function to set arg initialized to named numeric vector
  ## so resort to building text and parsing it
  ##
  data$fbody = fbody
  mll_txt = paste0(
    "make_function(alist(",
    paste(par.par.extra.txt, collapse=", "),
    "), fbody)")
  mll_fun = eval(parse(text=mll_txt), envir=data)
  rm(fbody, envir=data)

  ##
  ## create pulls function
  ##
  fbody = as.call(c(
    as.name("{"),
    list(
      bquote(if (is.null(names(par))) names(par) <- .(names(par))),
      quote(mapply(function(name, val) assign(name, val, pos=parent.frame(2)), names(par), par)),
      fbody.pulls
    )))

  ##
  ## (memo) pulls_fun = evalq(make_function(alist(par=), fbody), envir=data)
  ##
  ## could not find a way with make_function to set arg initialized to named numeric vector
  ## so resort to building text and parsing it
  ##
  data$fbody = fbody
  pulls_txt = paste0(
    "make_function(alist(",
    paste(par.par.extra.txt, collapse=", "),
    "), fbody)")
  pulls_fun = eval(parse(text=pulls_txt), envir=data)
  rm(fbody, envir=data)

  ##
  ## return just mll function when formula is of type "~ <residual formula>"
  ##
  if (is.null(fun_expr))
    return(
      list(
        fun_mll = mll_fun,
        fun_pulls = pulls_fun
      ))

  ##
  ## build function corresponding to formula RHS
  ##
  if (!is.null(rhs_vars)) {
    ##--- take all variables of RHS of formula from the caller
    fun_args = rhs_vars
  } else {
    ##
    ## all.vars(formula) does not return variables that are initialized
    ## brute force procedure to get rid of all (or most) initializations
    ## - one possible failure: variables initialized with expressions containing parenthesis
    ##
    fun_expr_str = gsub("\\s+", " ", paste(deparse(fun_expr), collapse=""))
    fun_expr_str = gsub("\\s*=\\s*[[:alnum:]._+*/^-]+", "", fun_expr_str)
    fun_args = all.vars(parse(text=fun_expr_str))
  }

  fun_args = setdiff(fun_args, c(names(par), names(list(...))))
  fun_args.txt = lapply(fun_args, function(arg) paste0(arg, "="))

  fbody = as.call(c(
    as.name("{"),
    list(
      bquote(if (is.null(names(par))) names(par) <- .(names(par))),
      quote(mapply(function(name, val) assign(name, val, pos=parent.frame(2)), names(par), par)),
      fun_expr)))

  ##
  ## if formula is of type "y ~ f(x, par)" then it is possible to get f(x, par)
  ## therefore return to caller f(x, par) in two formats
  ## - fun_par(x, par, non_fitted_par) where all parameters are passed in a single numeric vector
  ## - fun(x, p1, p2, non_fitted_par) where the parameters are passed one per argument
  ##

  ##
  ## define model function with fit parameters one numeric vector
  ##
  fun_txt = paste0(
    "make_function(alist(",
    paste(c(fun_args.txt, par.par.extra.txt), collapse=", "),
    "), fbody)")
  model_fun_par = eval(parse(text=fun_txt))

  ##
  ## define model function with fit parameters in separate args
  ##
  fun_txt = paste0(
    "make_function(alist(",
    paste(c(fun_args.txt, par.extra.txt), collapse=", "),
    "), fun_expr)")
  model_fun = eval(parse(text=fun_txt))

  ##--- return mll function and model function in two formats
  invisible(list(
    fun_mll = mll_fun,
    fun_pulls = pulls_fun,
    fun = model_fun,
    fun_par = model_fun_par
    ))
}

#' Function Minimization and model fitting with Minuit2
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
#' @param lh desired likelihood for the measurement errors, one of:
#'  \describe{
#'    \item{\code{Gaussian}:}{measurements have Gaussian errors}
#' }
#'
#' @param rhs_vars character vector with the name of all variables in the right-hand
#'   side (RHS) of formula. The code attempts to get the variables automatically from the formula
#'   expression if this argument is not provided, but may fail. For simple formula, most often
#'   there is no need to provide this argument.
#'
#' @param ... extra arguments for the model, weights and error formulas
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
#'    \item{\code{fun}:}{function corresponding to model with parameters initialized to the fitted values}
#'    \item{\code{fun_par}:}{function corresponding to model with parameters initialized to the fitted values, all parameters are passed as a possibly-named numeric vector in the last argument, which is named \code{par}}
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
#' @seealso rminuit2 rminuit2_par
#'
#' @examples
#' #
#' # simulate model y = a*exp(-x/b) + k
#' #
#' x = seq(0, 1, length.out=31)
#' y.func = function(x, norm, tau, const) norm*exp(-x/tau) + const
#'
#' #
#' # simulate data with Gaussian errors for specific model
#' # the model includes norm and tau but not const, which is a fixed parameter
#' #
#' model.par = c(norm=2.3, tau=0.47)
#' model.extra.const = 4.7
#' y.err = 0.01
#' y = do.call(y.func, c(list(x), model.par, const=model.extra.const)) + rnorm(sd=y.err, n=length(x))
#'
#' # fit model on data, ask to compute Minos errors too
#' fit.rc = rminuit2_expr(y ~ norm*exp(-x/tau) + const, c(norm=1, tau=10),
#'   data=data.frame(x=x, y=y, y.err=y.err), errors=y.err,
#'   lh="Gaussian", opt="hm", const=model.extra.const)
#'
#' # chi square / number of degrees of freedom
#' cbind(chisq=fit.rc$chisq, ndof=fit.rc$ndof)
#'
#' # fitted parameters and their estimated uncertainties
#' cbind(model.value=model.par, value=fit.rc$par, error=fit.rc$err,
#'       minos_pos=fit.rc$err_minos_pos, minos_neg=fit.rc$err_minos_neg)
#'
#' # parameters' correlation matrix
#' fit.rc$cor
#'
#' @export
#' 
rminuit2_expr = function(formula, start, data=NULL, weights=NULL, errors=NULL,
                         err=NULL, lower=NULL, upper=NULL, fix=NULL,
                         lh=c("Gaussian"),
                         opt="h", maxcalls=0L, nsigma=1, rhs_vars=NULL, ...) {
  rc = switch(
    tolower(lh[1]),
    
    gaussian = eval(substitute(
      rminuit2_make_gaussian_mll(formula=formula, par=start, data=data, weights=weights, errors=errors, rhs_vars=rhs_vars, ...))),
    
    stop("bad 'lh', must be one of: 'Gaussian' (more will be added)")
  )
  
  rc.fit = rminuit2_par(mll=rc$fun_mll, start=start, err=err, lower=lower, upper=upper, fix=fix, opt=opt,
                        maxcalls=maxcalls, nsigma=nsigma, ...)

  ##--- build function with parameters set to fitted values, all parameters passed in one numeric vector
  fun_args = formals(rc$fun_par)

  formals_txt = paste0(
    "formals(rc$fun_par) = alist(",
    paste0(ifelse(
      names(fun_args) == "par",
      paste0(names(fun_args), "=c(", paste0(names(rc.fit$par), "=", rc.fit$par, collapse=", "), ")"),
      paste0(names(fun_args), "=", fun_args)),
      collapse=", "),
    ")")

  eval(parse(text=formals_txt))

  ##--- build function with parameters set to fitted values, each parameter is passed on one arg
  fun_args = formals(rc$fun)

  formals_txt = paste0(
    "formals(rc$fun) = alist(",
    paste0(ifelse(
      names(fun_args) %in% names(rc.fit$par),
      paste0(names(fun_args), "=", rc.fit$par[names(fun_args)]),
      paste0(names(fun_args), "=", fun_args)),
      collapse=", "),
    ")")

  eval(parse(text=formals_txt))

  ##--- build pulls function with parameters set to fitted values
  fun_args = formals(rc$fun_pulls)

  formals_txt = paste0(
    "formals(rc$fun_pulls) = alist(",
    paste0(ifelse(
      names(fun_args) == "par",
      paste0(names(fun_args), "=c(", paste0(names(rc.fit$par), "=", rc.fit$par, collapse=", "), ")"),
      paste0(names(fun_args), "=", fun_args)),
      collapse=", "),
    ")")

  eval(parse(text=formals_txt))

  ##--- build mll function with parameters set to fitted values
  fun_args = formals(rc$fun_mll)

  formals_txt = paste0(
    "formals(rc$fun_mll) = alist(",
    paste0(ifelse(
      names(fun_args) == "par",
      paste0(names(fun_args), "=c(", paste0(names(rc.fit$par), "=", rc.fit$par, collapse=", "), ")"),
      paste0(names(fun_args), "=", fun_args)),
      collapse=", "),
    ")")

  eval(parse(text=formals_txt))

  ##--- number of observations
  pulls.val = rc$fun_pulls()
  nobs = length(pulls.val)

  ##
  ## chisq, beware it is not 2*fval (mll minimum)
  ## since possibly parameters' dependent errors are included
  ## in the likelihood
  ##
  chisq = sum(pulls.val^2)

  ##--- number of degrees of freedom
  ndof = nobs - length(start) + sum(fix != 0)
  
  invisible(c(
    rc.fit,
    chisq = chisq,
    ndof = ndof,
    nobs = nobs,
    fun_mll = rc$fun_mll,
    fun_pulls = rc$fun_pulls,
    fun = rc$fun,
    fun_par = rc$fun_par
  ))
}
