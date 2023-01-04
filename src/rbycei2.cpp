#include <R.h>
#include "eiutil.h"
#include <Rcpp.h>
using namespace Rcpp;

Rcpp::List rbycei_fcn2(NumericMatrix alphamatrix,
                       NumericMatrix betaarray,
                       NumericVector TT,
                       NumericVector XX,
                       NumericVector tuneA,
                       NumericVector tuneB,
                       int NG,
                       int NP,
                       int Precincts,
                       double Lambda1,
                       double Lambda2,
                       int Sample,
                       int Thin,
                       int Burnin,
                       int Verbose,
                       int Savebeta,
                       int RR,
                       SEXP usr_fcn,
                       SEXP usr_env,
                       int usr_len,
                       CharacterVector betanames
){

  int counter = 0;

  int ii, rr, cc, tt, qq, kk;
  NumericMatrix lbm;
  Rcpp::List output_list;
  SEXP alpha_dim, beta_dim, TT_dim, RR_dim;
  double aprop, acurr, asumm1, lbm_rc, aprop_ll, acurr_ll;
  double bprop, bprop_ref, bcurr, bcurr_ref, tt_ci, tt_Ci, bprop_x, bprop_ref_x, bcurr_x, bcurr_ref_x, ulim, bcurr_ll, bprop_ll;

  int iters = Burnin + Sample * Thin;

  List usr_list = List(4);

  double hldr;
  NumericVector usr(usr_len);
  NumericVector a_acc(NG * NP);
  NumericVector b_acc(NG * (NP - 1) * Precincts);
  NumericMatrix ccount(NG, NP);

  for(qq = 0; qq < NG*NP; ++qq){
    a_acc[qq] = 0;
  }
  for(qq = 0; qq < NG*(NP-1)*Precincts;++qq){
    b_acc[qq] = 0;
  }

  NumericMatrix ccount_draws(Sample , NG*NP);
  NumericMatrix usr_draws(Sample, usr_len);
  NumericMatrix a_draws(Sample, NG*NP);
  NumericMatrix b_draws;


  if(Savebeta == 0){
    b_draws = NumericMatrix(Sample, NG*NP*Precincts);
  }

  GetRNGstate();


  for(kk = 0; kk < iters; ++kk){
    for(ii = 0; ii < Precincts; ++ii){
      for(rr = 0; rr < NG; ++rr){
        for(cc = 0; cc < (NP - 1); ++cc){
          bcurr = betaarray[rr + NG*cc + NG*NP*ii];
          bcurr_ref = betaarray[rr + NG*(NP-1) + NG*NP*ii];
          ulim = bcurr + bcurr_ref;
          bprop = rnorm(bcurr, tuneB[rr + NG*cc + NG*(NP-1)*ii]);
          bprop_ref = ulim - bprop;


          if(bprop > 0 && bprop < ulim){
            tt_ci = TT[cc + NP*ii];
            tt_Ci = TT[(NP - 1) + NP*ii];

            bcurr_x = 0; bcurr_ref_x = 0;
            for(qq = 0; qq < NG; ++qq){
              bcurr_x += betaarray[qq + NG*cc + NG*NP*ii]* XX[qq + NG*ii];
              bcurr_ref_x += betaarray[qq + NG*(NP-1) + NG*NP*ii] * XX[qq + NG*ii];
            }


            bprop_x = bcurr_x - bcurr*XX[rr + NG*ii] + bprop*XX[rr + NG*ii];
            bprop_ref_x =  bcurr_ref_x - bcurr_ref*XX[rr + NG*ii] + bprop_ref*XX[rr + NG*ii];

            bprop_ll = beta_ll(bprop, bprop_ref, alphamatrix[rr + NG*cc], alphamatrix[rr + NG*(NP-1)],tt_ci, tt_Ci, bprop_x, bprop_ref_x);
            bcurr_ll = beta_ll(bcurr, bcurr_ref,  alphamatrix[rr + NG*cc], alphamatrix[rr + NG*(NP-1)], tt_ci, tt_Ci, bcurr_x, bcurr_ref_x);

            /* Rprintf("%f %f %f %f %f\n", bprop, bprop_ref, bprop_x, bprop_ref_x, bprop_ll - bcurr_ll);
             */
            if(acc_tog(bprop_ll, bcurr_ll) == 1){
              betaarray[rr + NG*cc + NG*NP*ii] = bprop;
              betaarray[rr + NG*(NP-1) + NG*NP*ii] = bprop_ref;
              b_acc[rr + NG*cc + NG*(NP - 1)*ii] += 1;

            }
          }
        }
      }
    }


    lbm = logbetamat(betaarray, NG, NP, Precincts);

    for(rr = 0; rr < NG; ++rr){
      for(cc = 0; cc < NP; ++cc){

        acurr = alphamatrix[rr + NG*cc];
        aprop = rnorm(acurr, tuneA[rr + NG*cc]);

        lbm_rc = lbm[rr + NG*cc];
        asumm1 = 0;
        for(tt = 0; tt < NP; ++tt){
          asumm1 += alphamatrix[rr + NG*tt];
        }
        asumm1 = asumm1 - acurr;

        if(aprop > 0){
          aprop_ll = alpha_ll(aprop, asumm1, lbm_rc, Lambda1, Lambda2, Precincts);
          acurr_ll = alpha_ll(acurr, asumm1, lbm_rc, Lambda1, Lambda2, Precincts);
          /*Rprintf("%f %f \n", aprop, aprop_ll - acurr_ll);
           */
          if(acc_tog(aprop_ll, acurr_ll) == 1){
            alphamatrix[rr + NG*cc] = aprop;
            a_acc[rr + NG*cc] += 1;

          }
        }
      }
    }

    if(kk >= Burnin && ((kk % Thin) == 0)){

      usr_list[0] = alphamatrix;
      usr_list[1] = betaarray;
      usr_list[2] = TT;
      usr_list[3] = RR;

      ccount = cellcount(betaarray, RR, NG, NP, Precincts);

      usr = usr_fun_eval(usr_fcn, usr_list, usr_env, usr_len);

      for(qq = 0; qq < usr_len; ++qq){
        usr_draws[counter + qq*Sample] = usr[qq];
      }


      for(qq = 0; qq < NP*NG; ++qq){
        a_draws[counter + qq*Sample] = alphamatrix[qq];
        ccount_draws[counter + qq*Sample] = ccount[qq];
      }

      if(Savebeta == 0){
        for(qq = 0; qq < NP*NG*Precincts; ++qq){
          b_draws[counter + qq*Sample] = betaarray[qq];
        }
      }

      // if(Savebeta == 2){
      //   hldr = write_beta(betaarray, betanames);
      // }

      counter += 1;
    }

    if(Verbose > 0 && kk % Verbose == 0){
      Rprintf("\n MCMC iteration %i of %i \n", kk + 1, iters);
    }
    R_CheckUserInterrupt();
  }


  for(qq = 0; qq < NG*NP; ++qq){
    a_acc[qq] = a_acc[qq]/iters;
  }
  for(qq = 0; qq < NG*(NP-1)*Precincts;++qq){
    b_acc[qq] = b_acc[qq]/iters;
  }


  if(Savebeta==0){
    Rcpp::List output_list = List::create(
      a_draws,
      b_draws,
      a_acc,
      b_acc,
      ccount_draws,
      usr_draws
      );
  }else{
    Rcpp::List output_list = List::create(
      a_draws,
      a_acc,
      b_acc,
      ccount_draws,
      usr_draws
    );
  }

  PutRNGstate();

  return output_list;
}
