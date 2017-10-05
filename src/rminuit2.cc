//
// rminuit2
//
// by Alberto Lusiani
// adapted from Minuit2
// Copyright (c) 2017, Alberto Lusiani
// license: GNU LGPL 2.1
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// $Id$

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "rminuit2.h"
#include "evaluate.h"
#include <Rcpp.h>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include "Minuit2/FCNBase.h"

#ifdef  _MSC_VER
  #define inline  __inline
#endif/*_MSC_VER*/

#if     defined(USE_SSE) && defined(__SSE2__) && LBFGS_FLOAT == 64
  /* Use SSE2 optimization for 64bit double precision. */
  #include "arithmetic_sse_double.h"

#elif   defined(USE_SSE) && defined(__SSE__) && LBFGS_FLOAT == 32
  /* Use SSE optimization for 32bit float precision. */
  #include "arithmetic_sse_float.h"

#else
  /* No CPU specific optimization. */
  #include "arithmetic_ansi.h"

#endif

using namespace Rcpp;

using namespace ROOT::Minuit2;

inline bool contained(char cmp, const std::string& str) {
  return str.find(cmp) != std::string::npos;
}

//
// mimimizer function that executes function to minimize
//

class ExtFcn : public FCNBase {

public:
  ExtFcn(Rcpp::EvalBase* const ev_fn) : fExtFcn(ev_fn) {}
  
  virtual ~ExtFcn() {}
  
  virtual double operator()(const std::vector<double>& par) const {
    const Rcpp::NumericVector rc( fExtFcn->eval(Rcpp::wrap(par)) );
#if 0    
    if (rc.size() != 1) {
      //+++ throw exception, function returned vector of length rc.size()
    }
#endif
    return rc.size() == 1 ? rc[0] : NA_REAL;
  }
  
  virtual double Up() const {return 0.5;}
  
private:
  Rcpp::EvalBase* const fExtFcn;
};

//
// rminuit2_cpp
//

//[[Rcpp::export]]
Rcpp::List rminuit2_cpp(
  SEXP fn,
  Rcpp::NumericVector par,
  Rcpp::NumericVector err,
  Rcpp::NumericVector lower,
  Rcpp::NumericVector upper,
  Rcpp::IntegerVector fix,
  Rcpp::StringVector opt,
  SEXP envir,
  Rcpp::IntegerVector maxcalls
  )
{
  //--- Pointer to abstract base classes
  Rcpp::EvalBase *ev_fn(NULL);

  // Assign class based on object type
  if (TYPEOF(fn) == EXTPTRSXP) {
    //
    // Non-standard mode: we are being passed an external pointer
    // So assign pointers using external pointer in calls SEXP
    //
    ev_fn = new Rcpp::EvalCompiled(fn, envir);
  } else {
    //
    // Standard mode: env_ is an env, functions are R objects
    //
    ev_fn = new Rcpp::EvalStandard(fn, envir);    // So assign R functions and environment
  }

  //--- create function to be minimized with Minuit2
  ExtFcn fFCN(ev_fn);

  // int npar(par.size());

  auto dpar( as< std::vector<double> >(par) );
  auto derr( as< std::vector<double> >(err) );
  auto dlower( as< std::vector<double> >(lower) );
  auto dupper( as< std::vector<double> >(upper) );
  auto ifix( as< std::vector<int> >(fix) );
  auto sopt( as< std::string >(opt(0)) );
  int imaxcalls( maxcalls(0) );
  
  // auto par_names( as< std::vector<std::string> >(par.names()) );
  auto par_names = as< Rcpp::StringVector >(par.names());
  
  //--- Minuit2 parameters
  MnUserParameters upar;
  
  //--- init upar names, start values, estimated errors, lower and upper limits
  for(auto i=0; i<par.size(); i++) {
    upar.Add(as< const char * >(par_names(i)), dpar[i], derr[i] );
    if (dlower[i] != R_NegInf) upar.SetLowerLimit(i, dlower[i]);
    if (dupper[i] != R_PosInf) upar.SetUpperLimit(i, dupper[i]);
    // if (dlower[i] != R_NegInf) Rcout << "lower limit par " << i << " " << dlower[i] << std::endl;    
    // if (dupper[i] != R_PosInf) Rcout << "upper limit par " << i << " " << dupper[i] << std::endl;    
    if (ifix[i] != 0) upar.Fix(i);
    // if (ifix[i] != 0) Rcout << "fixed par " << i << std::endl;    
  }

  FunctionMinimum* fminp(0);

  if (contained('0', sopt)) {
    MnMigrad migrad(fFCN, upar, MnStrategy(0));
    fminp = new FunctionMinimum(migrad(imaxcalls));
  }

  if (contained('1', sopt)) {
    MnMigrad migrad(fFCN, fminp ? fminp->UserState() : upar, MnStrategy(1));
    if (fminp) {
      *fminp = migrad(imaxcalls);
    } else {
      fminp = new FunctionMinimum(migrad(imaxcalls));
    }
  }

  if (contained('2', sopt)) {
    MnMigrad migrad(fFCN, fminp ? fminp->UserState() : upar, MnStrategy(2));
    if (fminp) {
      *fminp = migrad(imaxcalls);
    } else {
      fminp = new FunctionMinimum(migrad(imaxcalls));
    }
  }

  auto migrad_first_failed(false);
  if (!contained('0', sopt) && !contained('1', sopt) && !contained('2', sopt)) {
    //--- default strategy
    MnMigrad migrad(fFCN, upar);
    fminp = new FunctionMinimum(migrad(imaxcalls));

    if (!fminp->IsValid()) {
      migrad_first_failed = true;
      //--- try with strategy 2 if failed
      Rcpp::warning("First migrad call failed, try with strategy=2");
      MnMigrad migrad2(fFCN, fminp->UserState(), MnStrategy(2));
      *fminp = migrad2(imaxcalls);
    }

  }
  
  FunctionMinimum& min(*fminp);
  
  if (contained('h', sopt)) {
    //--- run hesse to compute Hessian and errors
    MnHesse hesse;
    hesse(fFCN, min);
  }
  
  // Rcout << "UserParameters " << min.UserParameters() << std::endl;
  // Rcout << "UserCovariance " << min.UserCovariance() << std::endl;
  
  std::vector<double> minos_pos_err;
  std::vector<double> minos_neg_err;
  if (contained('m', sopt)) {
    MnMinos Minos(fFCN, min);
    for(auto i=0; i<par.size(); i++) {
      minos_pos_err.push_back(Minos(i).second);
      minos_neg_err.push_back(Minos(i).first);
    }
  }
  
  // Rcout << "minimum: " << min << std::endl;

  Rcpp::NumericMatrix par_cov(par.size(), par.size());
  for(auto i=0; i<par.size(); i++) {
    for(auto j=0; j<par.size(); j++) {
      par_cov(i, j) = min.UserCovariance()(i, j);
    }
  }
  
  //--- return list
  Rcpp::List rc;

  rc["par"] = min.UserParameters().Params();
  rc["err"] = min.UserParameters().Errors();
  rc["cov"] = par_cov;

  if (contained('m', sopt)) {
    rc["err_minos_pos"] = minos_pos_err;
    rc["err_minos_neg"] = minos_neg_err;
  }

  rc["fval"] = min.Fval();

  rc["IsValid"] = min.IsValid();
  rc["IsValidFirstInvocation"] = !migrad_first_failed;
  rc["HasValidParameters"] = min.HasValidParameters();
  rc["HasValidCovariance"] = min.HasValidCovariance();
  rc["HasAccurateCovar"] = min.HasAccurateCovar();
  rc["HasPosDefCovar"] = min.HasPosDefCovar();
  rc["HasMadePosDefCovar"] = min.HasMadePosDefCovar();
  rc["HesseFailed"] = min.HesseFailed();
  rc["HasCovariance"] = min.HasCovariance();
  rc["IsAboveMaxEdm"] = min.IsAboveMaxEdm();
  rc["HasReachedCallLimit"] = min.HasReachedCallLimit();

  return rc;
}
