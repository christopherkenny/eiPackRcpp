#include <Rcpp.h>
using namespace Rcpp;
using namespace R;

// [[Rcpp::export]]
double beta_ll(double beta, double beta_ref, double alpha, double alpha_ref,
               double tt_ci, double tt_Ci, double beta_X, double beta_ref_X) {
  double ll;

  ll = tt_ci * std::log(beta_X) + tt_Ci * std::log(beta_ref_X) +
    (alpha - 1) * std::log(beta) + (alpha_ref - 1) * std::log(beta_ref);
  return ll;
}

double alpha_ll (double alpha, double sum_alphaminus, double lbmat, double lambda1, double lambda2, int prec){

  double ll;
  ll = (lambda1 - 1) * log(alpha) - lambda2 * alpha +
    prec * (lgammafn((alpha + sum_alphaminus)) - lgammafn(alpha)) + lbmat * alpha;

  return ll;
}

int acc_tog (double prop_ll, double curr_ll){

  int toggle = 0;
  double log_rat;
  log_rat = prop_ll - curr_ll;
  double r = runif(0, 1)[0];
  if (std::log(r) < log_rat) {
    toggle = 1;
  }
  return toggle;
}

NumericMatrix cellcount(NumericVector betas, NumericVector RR, int NG, int NP, int Precincts){

  double tmp = 0.0;

  NumericMatrix ret_val(NG, NP);

  for(int rr = 0; rr < NG; rr++){
    for(int cc = 0; cc < NP; cc++){
      tmp = 0.0;
      for(int ii = 0; ii < Precincts; ii++){
        tmp += betas[rr + NG * cc + NG * NP * ii] * RR[rr + NG*ii];
      }
      ret_val[rr + NG * cc] = tmp;
    }
  }

  return ret_val;
}

NumericMatrix logbetamat(NumericMatrix aa, int NG, int NP, int Precincts){

  double tmp = 0.0;

  NumericMatrix ret_val(NG, NP);

  for(int rr = 0; rr < NG; ++rr){
    for(int cc = 0; cc < NP; ++cc){
      tmp = 0.0;
      for(int ii = 0; ii < Precincts; ++ii){
        tmp += log(aa[rr + NG * cc + NG * NP * ii]);
      }
      ret_val[rr + NG * cc] = tmp;
    }
  }

  return ret_val;
}

//  (NumericVector betaarray, CharacterVector filenames)
// SEXP write_beta (SEXP betaarray, SEXP filenames){
//
//   int ret=0;
//   R_len_t nn, ii;
//   nn = length(filenames);
//   for(ii = 0; ii < nn; ++ii){
//     char tmp[500];
//     sprintf(tmp, "echo \"%.16f\" | gzip >>  %s", REAL(betaarray)[ii], CHAR(STRING_ELT(filenames,ii)));
//     ret=system(tmp);
//   }
//
//
//   R_CheckUserInterrupt();
//
//   return(R_NilValue);
// }
//
//
// SEXP usr_fun_eval(Function usr_fun, SEXP cur_values, SEXP usr_env, int usr_len){
//
//   SEXP R_fcall, usr_fcn_output;
//   if(!isFunction(usr_fun)) error("`usr_fun' must be a function");
//   if(!isEnvironment(usr_env)) error("`usr_env' must be an environment");
//   PROTECT(R_fcall = lang2(usr_fun, R_NilValue));
//   SETCADR(R_fcall, cur_values);
//   PROTECT(usr_fcn_output = allocVector(REALSXP, usr_len));
//   usr_fcn_output = eval(R_fcall, usr_env);
//   UNPROTECT(2);
//   return(usr_fcn_output);
// }
