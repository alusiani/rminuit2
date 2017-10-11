//
// rminuit2.cc
//
// by Alberto Lusiani
// Copyright (c) 2017, Alberto Lusiani
// license: GNU LGPL 3.0
//
// C++ code for invoking Minuit2 procedures for the R package rminuit2
// - reliaes on the Rcpp package
//

// $Id$

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <memory>

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

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

using namespace Rcpp;

using namespace ROOT::Minuit2;

inline bool contained(char cmp, const std::string& str) {
  return str.find(cmp) != std::string::npos;
}

//
// mimimizer function that executes function to minimize
//

class FcnRcppAdapter : public FCNBase {

public:
  FcnRcppAdapter(Rcpp::EvalBase* const ev_fn, double up = 0.5) :
    fFunc(ev_fn),
    fUp (up)
  {}

  ~FcnRcppAdapter() {}

  double operator()(const std::vector<double>& par) const {
    const Rcpp::NumericVector rc( fFunc->eval(Rcpp::wrap(par)) );
#if 0
    if (rc.size() != 1) {
      //+++ throw exception, function returned vector of length rc.size()
    }
#endif
    return rc.size() == 1 ? rc[0] : NA_REAL;
  }

#if 0
  double operator()(const double* par) const {
    const Rcpp::NumericVector rc( fFunc->eval(Rcpp::wrap(par)) );
#if 0
    if (rc.size() != 1) {
      //+++ throw exception, function returned vector of length rc.size()
    }
#endif
    return rc.size() == 1 ? rc[0] : NA_REAL;
  }
#endif

  double Up() const {return fUp;}

  void SetErrorDef(double up) { fUp = up; }

private:
  Rcpp::EvalBase* const fFunc;
  double fUp;

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
  Rcpp::IntegerVector maxcalls,
  Rcpp::NumericVector nsigma
  )
{
  //--- Pointer to abstract base classes, autodestructs when out of scope
  std::unique_ptr<Rcpp::EvalBase> ev_fn;

  // Assign class based on object type
  if (TYPEOF(fn) == EXTPTRSXP) {
    //
    // Non-standard mode: we are being passed an external pointer
    // So assign pointers using external pointer in calls SEXP
    //
    ev_fn = std::unique_ptr<Rcpp::EvalBase>(new Rcpp::EvalCompiled(fn, envir));
  } else {
    //
    // Standard mode: env_ is an env, functions are R objects
    // So assign R functions and environment
    //
    ev_fn = std::unique_ptr<Rcpp::EvalBase>(new Rcpp::EvalStandard(fn, envir));
  }

  //--- create function to be minimized with Minuit2
  FcnRcppAdapter fFcn(ev_fn.get());

  auto dpar( as< std::vector<double> >(par) );
  auto derr( as< std::vector<double> >(err) );
  auto dlower( as< std::vector<double> >(lower) );
  auto dupper( as< std::vector<double> >(upper) );
  auto ifix( as< std::vector<int> >(fix) );
  auto sopt( as< std::string >(opt(0)) );
  int imaxcalls( maxcalls(0) );
  double dnsigma( nsigma[0] );

  // auto par_names( as< std::vector<std::string> >(par.names()) );
  auto par_names = as< Rcpp::StringVector >(par.names());

  //--- Minuit2 parameters
  MnUserParameters upar;

  //--- init upar names, start values, estimated errors, lower and upper limits
  size_t npar_nonfixed(par.size());

  for(unsigned int i=0; i<par.size(); i++) {
    upar.Add(as< const char * >(par_names(i)), dpar[i], derr[i] );
    if (dlower[i] != R_NegInf) upar.SetLowerLimit(i, dlower[i]);
    if (dupper[i] != R_PosInf) upar.SetUpperLimit(i, dupper[i]);
    // if (dlower[i] != R_NegInf) Rcout << "lower limit par " << i << " " << dlower[i] << std::endl;
    // if (dupper[i] != R_PosInf) Rcout << "upper limit par " << i << " " << dupper[i] << std::endl;
    if (ifix[i] != 0) {
      upar.Fix(i);
      npar_nonfixed--;
    }
    // if (ifix[i] != 0) Rcout << "fixed par " << i << std::endl;
  }

  std::unique_ptr<FunctionMinimum> fminp;

  if (contained('0', sopt)) {
    MnMigrad migrad(fFcn, upar, MnStrategy(0));
    fminp = std::unique_ptr<FunctionMinimum>(new FunctionMinimum(migrad(imaxcalls)));
  }

  if (contained('1', sopt)) {
    MnMigrad migrad(fFcn, fminp ? fminp->UserState() : upar, MnStrategy(1));
    if (fminp != nullptr) {
      *fminp = migrad(imaxcalls);
    } else {
      fminp = std::unique_ptr<FunctionMinimum>(new FunctionMinimum(migrad(imaxcalls)));
    }
  }

  if (contained('2', sopt)) {
    MnMigrad migrad(fFcn, fminp ? fminp->UserState() : upar, MnStrategy(2));
    if (fminp != nullptr) {
      *fminp = migrad(imaxcalls);
    } else {
      fminp = std::unique_ptr<FunctionMinimum>(new FunctionMinimum(migrad(imaxcalls)));
    }
  }

  auto migrad_first_failed(false);
  if (!contained('0', sopt) && !contained('1', sopt) && !contained('2', sopt)) {
    //--- default strategy
    MnMigrad migrad(fFcn, upar);
    fminp = std::unique_ptr<FunctionMinimum>(new FunctionMinimum(migrad(imaxcalls)));

    if (!fminp->IsValid()) {
      migrad_first_failed = true;
      //--- try with strategy 2 if failed
      // Rcpp::warning("First migrad call failed, try with strategy=2");
      MnMigrad migrad2(fFcn, fminp->UserState(), MnStrategy(2));
      *fminp = migrad2(imaxcalls);
    }
  }

  //--- get a reference to the object as syntactic sugar
  FunctionMinimum& min(*fminp);

  if (contained('h', sopt)) {
    //--- run hesse to compute Hessian and errors
    MnHesse hesse;
    hesse(fFcn, min);
  }

  // Rcout << "UserParameters " << min.UserParameters() << std::endl;
  // Rcout << "UserCovariance " << min.UserCovariance() << std::endl;

  std::vector<double> minos_pos_err;
  std::vector<double> minos_neg_err;
  std::vector<bool> minos_pos_err_valid;
  std::vector<bool> minos_neg_err_valid;
  if (contained('m', sopt)) {
    MnMinos Minos(fFcn, min);
    fFcn.SetErrorDef(dnsigma * dnsigma);
    for(unsigned int i=0; i<par.size(); i++) {
      if (ifix[i]) {
        minos_pos_err.push_back(0);
        minos_neg_err.push_back(0);
        minos_pos_err_valid.push_back(true);
        minos_neg_err_valid.push_back(true);
      } else {
        ROOT::Minuit2::MinosError me(Minos.Minos(i));
        minos_pos_err.push_back(me.Upper());
        minos_neg_err.push_back(me.Lower());
        minos_pos_err_valid.push_back(me.UpperValid());
        minos_neg_err_valid.push_back(me.LowerValid());
      }
    }
  }

  //--- get covariance with caution (nrow=0 when invalid)
  Rcpp::NumericMatrix par_cov(npar_nonfixed, npar_nonfixed);
  for(unsigned int i=0; i<npar_nonfixed && i<min.UserCovariance().Nrow(); i++) {
    for(unsigned int j=0; j<npar_nonfixed && j<min.UserCovariance().Nrow(); j++) {
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
    rc["err_minos_pos_valid"] = minos_pos_err_valid;
    rc["err_minos_neg_valid"] = minos_neg_err_valid;
  }

  rc["fval"] = min.Fval();
  rc["Edm"] = min.Edm();

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
