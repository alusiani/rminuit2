//
// rminuit2.cc
//
// Copyright (c) 2017, 2026 Alberto Lusiani
// license: GNU LGPL 3.0
//

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>

#include <Rcpp.h>
#include "./evaluate.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/FCNBase.h"

using namespace Rcpp;
using namespace ROOT::Minuit2;

class FcnRcppAdapter : public ROOT::Minuit2::FCNBase {
public:
  // Store by value: Rcpp::Function is just a wrapper around a SEXP.
  // Copying it by value is cheap and handles internal SEXP protection.
  FcnRcppAdapter(SEXP fn, SEXP env) :
    fEval(fn, env) {}

  // Shallow copy constructor: Let Rcpp handles the 'fn' and 'evaluator' logic.
  FcnRcppAdapter(const FcnRcppAdapter& other) : 
    fEval(other.fEval), fErrorDef(other.fErrorDef) {}

  virtual ~FcnRcppAdapter() {}

  double operator()(const std::vector<double>& par) const {
    // Basic R-check without extra overhead
    if (fEval.getNbEvals() % 500 == 0) R_CheckUserInterrupt();

    // Create the parameter vector for the call
    // Rcpp::NumericVector r_par(par.begin(), par.end());
    Rcpp::NumericVector r_par = Rcpp::wrap(par);

    // Call through the evaluator. 
    return fEval.eval(r_par);
  }

  double Up() const { return fErrorDef; }
  void SetErrorDef(double errorDef) { fErrorDef = errorDef; }

private:
  // fEval is an instance of EvalStandard defined in evaluate.h
  mutable Rcpp::DE::EvalStandard fEval;
  double fErrorDef = 1.0;
};

// [[Rcpp::export]]
List rminuit2_cpp(
  SEXP fn,
  SEXP env,
  NumericVector par,
  NumericVector err,
  NumericVector lower,
  NumericVector upper,
  IntegerVector fix,
  StringVector opt,
  IntegerVector maxcalls,
  NumericVector nsigma
) {
  std::string sopt = as<std::string>(opt[0]);
  int imaxcalls = maxcalls[0];
  double dnsigma = nsigma[0];
  size_t npar = par.size();
  size_t npar_nonfixed = npar;
  StringVector par_names = par.names();

  MnUserParameters upar;
  for(size_t i = 0; i < npar; i++) {
    std::string name = (par_names.size() > i) ? as<std::string>(par_names[i]) : "p" + std::to_string(i);
    upar.Add(name, par[i], err[i]);
    
    if (lower[i] != R_NegInf) upar.SetLowerLimit(i, lower[i]);
    if (upper[i] != R_PosInf) upar.SetUpperLimit(i, upper[i]);
    
    if (fix[i] != 0) {
      upar.Fix(i);
      npar_nonfixed--;
    }
  }

  // Initialize adapter with the provided function and dimension
  FcnRcppAdapter fFcn(fn, env);

  std::unique_ptr<FunctionMinimum> fminp;

  auto run_migrad = [&](const MnUserParameterState& state, int strat = -1) {
    if (strat == -1) {
      // Uses the MnMigrad(FCNBase&, MnUserParameterState&) constructor
      // which defaults to Strategy 1 internally.
      MnMigrad migrad(fFcn, state);
      return migrad(imaxcalls);
    } else {
      // Uses the MnMigrad(FCNBase&, MnUserParameterState&, MnStrategy&) constructor
      MnMigrad migrad(fFcn, state, MnStrategy(strat));
      return migrad(imaxcalls);
    }
  };

  // Execution flow for strategies 0, 1, 2
  if (sopt.find('0') != std::string::npos) {
    fminp = std::make_unique<FunctionMinimum>(run_migrad(upar, 0));
  }
  if (sopt.find('1') != std::string::npos) {
    if (fminp) *fminp = run_migrad(fminp->UserState(), 1);
    else fminp = std::make_unique<FunctionMinimum>(run_migrad(upar, 1));
  }
  if (sopt.find('2') != std::string::npos) {
    if (fminp) *fminp = run_migrad(fminp->UserState(), 2);
    else fminp = std::make_unique<FunctionMinimum>(run_migrad(upar, 2));
  }
  
  bool migrad_first_failed = false;
  if (sopt.find_first_of("012") == std::string::npos) {
    // Now we use the lambda! Passing nothing (or -1) uses the default strategy.
    fminp = std::make_unique<FunctionMinimum>(run_migrad(upar)); 

    if (!fminp->IsValid()) {
      migrad_first_failed = true;
      // Fallback to Strategy 2 if the default failed
      *fminp = run_migrad(fminp->UserState(), 2);
    }
  }

  // Hesse
  if (sopt.find('h') != std::string::npos && fminp) {
    MnHesse hesse;
    hesse(fFcn, *fminp);
  }

  // Minos
  std::vector<double> minos_pos_err, minos_neg_err;
  std::vector<bool> minos_pos_valid, minos_neg_valid;
  if (sopt.find('m') != std::string::npos && fminp) {
    fFcn.SetErrorDef(dnsigma * dnsigma * fFcn.Up());
    MnMinos Minos(fFcn, *fminp);
    
    for(size_t i = 0; i < npar; i++) {
      if (fix[i]) {
        minos_pos_err.push_back(0); minos_neg_err.push_back(0);
        minos_pos_valid.push_back(true); minos_neg_valid.push_back(true);
      } else {
        MinosError me = Minos.Minos(i);
        minos_pos_err.push_back(me.Upper());
        minos_neg_err.push_back(me.Lower());
        minos_pos_valid.push_back(me.UpperValid());
        minos_neg_valid.push_back(me.LowerValid());
      }
    }
  }

  // Covariance Matrix
  NumericMatrix par_cov(npar_nonfixed, npar_nonfixed);
  if (fminp && fminp->HasCovariance()) {
    const MnUserCovariance& cov = fminp->UserCovariance();
    for(size_t i = 0; i < (size_t)npar_nonfixed && i < (size_t)cov.Nrow(); i++) {
      for(size_t j = 0; j < (size_t)npar_nonfixed && j < (size_t)cov.Nrow(); j++) {
        par_cov(i, j) = cov(i, j);
      }
    }
  }

  // Final List Assembly
  List rc = List::create(
    Named("par") = fminp->UserParameters().Params(),
    Named("err") = fminp->UserParameters().Errors(),
    Named("cov") = par_cov,
    Named("fval") = fminp->Fval(),
    Named("Edm") = fminp->Edm(),
    Named("IsValid") = fminp->IsValid(),
    Named("IsValidFirstInvocation") = !migrad_first_failed,
    Named("HasValidParameters") = fminp->HasValidParameters(),
    Named("HasValidCovariance") = fminp->HasValidCovariance(),
    Named("HasAccurateCovar") = fminp->HasAccurateCovar(),
    Named("HasPosDefCovar") = fminp->HasPosDefCovar(),
    Named("HasMadePosDefCovar") = fminp->HasMadePosDefCovar(),
    Named("HesseFailed") = fminp->HesseFailed(),
    Named("HasCovariance") = fminp->HasCovariance(),
    Named("IsAboveMaxEdm") = fminp->IsAboveMaxEdm(),
    Named("HasReachedCallLimit") = fminp->HasReachedCallLimit()
  );

  if (sopt.find('m') != std::string::npos) {
    rc["err_minos_pos"] = wrap(minos_pos_err);
    rc["err_minos_neg"] = wrap(minos_neg_err);
    rc["err_minos_pos_valid"] = wrap(minos_pos_valid);
    rc["err_minos_neg_valid"] = wrap(minos_neg_valid);
  }

  return rc;
}
