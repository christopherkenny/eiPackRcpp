#ifndef EIUTIL_H
#define EIUTIL_H

#include <Rcpp.h>
using namespace Rcpp;
using namespace R;

double beta_ll(double beta, double beta_ref, double alpha, double alpha_ref,
               double tt_ci, double tt_Ci, double beta_X, double beta_ref_X);

double alpha_ll(double alpha, double sum_alphaminus, double lbmat, double lambda1,
                double lambda2, int prec);

int acc_tog(double prop_ll, double curr_ll);

NumericMatrix cellcount(NumericVector betas, NumericVector RR, int NG, int NP,
                        int Precincts);

NumericMatrix logbetamat(NumericMatrix aa, int NG, int NP, int Precincts);

#endif
