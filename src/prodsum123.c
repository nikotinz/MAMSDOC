#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <mvtnormAPI.h>
#include <R_ext/Memory.h>



/* C_MAMS_pnorm is based on the function C_maxabsConditionalPvalue */
/* of the file scr/Distributions.c of the 'party' R package        */
/* of Torsten Hothorn et al. [DLC 202201-007]                      */

double C_MAMS_pnorm(double *lower, double *upper,
     int *infin, double *Sigma, int *Jcut, int *Jfull,
     int *maxpts, double *releps, double *abseps, double *tol) {

    int *nu, *inform, sub, iz, jz, rnd = 0;
    double *delta, *corr, *myerror, *prob, ans;

    /*            */
    /* univariate */
    /*            */

    /*  double pnorm5(double x, double mu, double sigma,  */
    /*                int lower_tail, int log_p)          */
    /*  where lower_tail=1 gives get P[X<x]               */

    if (Jcut[0] == 1){
        if(infin[0] == 0){
            return(pnorm(lower[0], 0.0, 1.0, 1, 0));
        }else{
            return(pnorm(upper[0], 0.0, 1.0, 1, 0) -
                   pnorm(lower[0], 0.0, 1.0, 1, 0));
        }
    }

    /*              */
    /* multivariate */
    /*              */

    nu = R_Calloc(1, int);
    myerror = R_Calloc(1, double);
    prob = R_Calloc(1, double);
    nu[0] = 0;
    inform = R_Calloc(1, int);
    delta = R_Calloc(Jcut[0], double);
    if (Jcut[0] == 2)
         corr = R_Calloc(1, double);
    else
         corr = R_Calloc(Jcut[0] + ((Jcut[0] - 2) * (Jcut[0] - 1))/2, double);


    /* mvtdst assumes the unique elements of the triangular */
    /* covariance matrix to be passes as argument CORREL    */
    for (iz = 0; iz < Jcut[0]; iz++) {
        delta[iz] = 0.0;
        for (jz = 0; jz < iz; jz++) {
            sub = (int) (jz + 1) + (double) ((iz - 1) * iz) / 2 - 1;
            corr[sub] = Sigma[iz*Jfull[0] + jz];
        }
    }

    /* mvtnorm_C_mvtdst (mvtnorm/include/mvtnormAPI.h):                     */
    /*                                                                      */
    /*     N      INTEGER, the number of variables.                         */
    /*     NU     INTEGER, the number of degrees of freedom.                */
    /*            If NU < 1, then an MVN probability is computed.           */
    /*     LOWER  DOUBLE PRECISION, array of lower integration limits.      */
    /*     UPPER  DOUBLE PRECISION, array of upper integration limits.      */
    /*     INFIN  INTEGER, array of integration limits flags:               */
    /*             if INFIN(I) < 0, Ith limits are (-infinity, infinity);   */
    /*             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];   */
    /*             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);    */
    /*             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].    */
    /*     CORREL DOUBLE PRECISION, array of correlation coefficients;      */
    /*            the correlation coefficient in row I column J of the      */
    /*            correlation matrixshould be stored in                     */
    /*               CORREL( J + ((I-2)/ times (I-1))/2 ), for J < I.       */
    /*            The correlation matrix must be positive semi-definite.    */
    /*     DELTA  DOUBLE PRECISION, array of non-centrality parameters.     */
    /*     MAXPTS INTEGER, maximum number of function values allowed. This  */
    /*            parameter can be used to limit the time. A sensible       */
    /*            strategy is to start with MAXPTS = 1000/N, and then       */
    /*            increase MAXPTS if ERROR is too large.                    */
    /*     ABSEPS DOUBLE PRECISION absolute error tolerance.                */
    /*     RELEPS DOUBLE PRECISION relative error tolerance.                */
    /*     ERROR  DOUBLE PRECISION estimated absolute error,                */
    /*            with 99% confidence level.                                */
    /*     VALUE  DOUBLE PRECISION estimated value for the integral         */
    /*     INFORM INTEGER, termination status parameter:                    */
    /*            if INFORM = 0, normal completion with ERROR < EPS;        */
    /*            if INFORM = 1, completion with ERROR > EPS and MAXPTS     */
    /*                           function vaules used; increase MAXPTS to   */
    /*                           decrease ERROR;                            */
    /*            if INFORM = 2, N > 1000 or N < 1.                         */
    /*            if INFORM = 3, correlation matrix not positive            */
    /*                           semi-definite.                             */

    mvtnorm_C_mvtdst(Jcut, nu, lower, upper, infin, corr, delta,
                     maxpts, abseps, releps, myerror, prob, inform, &rnd);

    /* inform == 0 means: everything is OK */
    switch (inform[0]) {
        case 0: break;
        case 1: warning("cmvnorm: completion with ERROR > EPS"); break;
        case 2: warning("cmvnorm: N > 1000 or N < 1");
                prob[0] = 0.0;
                break;
        case 3: warning("cmvnorm: correlation matrix not positive semi-definite");
                prob[0] = 0.0;
                break;
        default: warning("cmvnorm: unknown problem in MVTDST");
                 prob[0] = 0.0;
    }

    /* out */
    ans = prob[0];
    R_Free(corr); R_Free(myerror); R_Free(prob);
    R_Free(nu); R_Free(inform);
    return(ans);
}



/* R_MAMS_pnorm allows direct calls of C_MAMS_pnorm [DLC 202201-007] */

SEXP R_MAMS_pnorm(SEXP lower, SEXP upper, SEXP infin, SEXP Sigma,
                  SEXP Jcut, SEXP Jfull,
                  SEXP maxpts, SEXP releps, SEXP abseps, SEXP tol) {

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 1));

    /* make sure mvtdst has access to RNG */
    GetRNGstate();

    /* argument Jcut plays a key role:                      */
    /* elements of lower, upper, infin and Sigma located    */
    /* in positions >=Jcut are ignored                      */
    REAL(ans)[0] = C_MAMS_pnorm(REAL(lower), REAL(upper), INTEGER(infin),
        REAL(Sigma), INTEGER(Jcut), INTEGER(Jfull),
        INTEGER(maxpts), REAL(releps), REAL(abseps), REAL(tol));

    PutRNGstate();

    UNPROTECT(1);
    return(ans);
}



/* C_prodsum1 based on the function prodsum of  */
/* MAMS::mams version 1.4.2 [DLC 202201-007]    */

SEXP C_prodsum1(SEXP x2, SEXP l2, SEXP u2, SEXP r2, SEXP r02, SEXP r0diff2,
                SEXP Jfull2, SEXP K2, SEXP Sigma2,
                SEXP maxpts2, SEXP releps2, SEXP abseps2, SEXP tol2){

    /* inputs */
    int *Jfull, *K, *maxpts;
    double *x, *l, *u, *r, *r0, *r0diff, *Sigma, *releps, *abseps, *tol;
    Jfull = INTEGER(Jfull2);
    K = INTEGER(K2);
    maxpts = INTEGER(maxpts2);
    x = REAL(x2);
    l = REAL(l2);
    u = REAL(u2);
    r = REAL(r2);
    r0 = REAL(r02);
    r0diff = REAL(r0diff2);
    Sigma = REAL(Sigma2);
    releps = REAL(releps2);
    abseps = REAL(abseps2);
    tol = REAL(tol2);
    /* */
    SEXP ans2, lower2, upper2, L2, U2, Jcut2, infin2;
    int j=0, jj=0, nprotect=0, *Jcut, *infin;
    double *ans, *lower, *upper, *L, *U;
    PROTECT(ans2 = allocVector(REALSXP, 1)); nprotect++;
    ans = REAL(ans2);
    PROTECT(lower2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    lower = REAL(lower2);
    PROTECT(upper2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    upper = REAL(upper2);
    PROTECT(L2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    L = REAL(L2);
    PROTECT(U2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    U = REAL(U2);
    PROTECT(infin2 = allocVector(INTSXP, Jfull[0])); nprotect++;
    infin = INTEGER(infin2);
    PROTECT(Jcut2 = allocVector(INTSXP, 1)); nprotect++;
    Jcut = INTEGER(Jcut2);
    double p0 = 1.0, p1, p2, p3, insum = 0.0;

    /* make sure mvtdst has access to RNG */
    GetRNGstate();

    /* loop */
    p3 = 0.0;
    for(j = 0; j < Jfull[0]; j++){
        p0 *= dnorm(x[j], 0, 1, 0);
        p3 += sqrt(r0diff[j]) * x[j];
        p2  = sqrt(1 + r[j] / r0[j]);
        p1  = sqrt(r[j]) / r0[j] * p3;
        L[j] = p1 + l[j] * p2;
        U[j] = p1 + u[j] * p2;
        for(jj = 0; jj < (j+1); jj++){
            lower[jj] = L[jj];
            if(jj < j){
                upper[jj] = U[jj];
                infin[jj] = 2;
            }else{
                upper[jj] = L[jj];
                infin[jj] = 0;
            }
        }
        Jcut[0] = j+1;
        insum +=  C_MAMS_pnorm(lower, upper, infin, Sigma, Jcut, Jfull,
                               maxpts, releps, abseps, tol);
    }

    /* out */
    ans[0] = p0*pow(insum,K[0]);
    PutRNGstate();

    UNPROTECT(nprotect);
    return(ans2);
}



/* C_prodsum2 based on the function prodsum2 of */
/* MAMS::mams version 2.0.2 [DLC 202201-007]    */

SEXP C_prodsum2(SEXP x2, SEXP r2, SEXP r02, SEXP K2, SEXP u2,
                SEXP delta2, SEXP delta02, SEXP n2, SEXP sig2){

    /* inputs */
    int nprotect=0;
    double *x, *r, *r0, *u, *delta, *delta0, *n, *sig, *K;
    x = REAL(x2);
    r = REAL(r2);
    r0 = REAL(r02);
    u = REAL(u2);
    delta = REAL(delta2);
    delta0 = REAL(delta02);
    n = REAL(n2);
    sig = REAL(sig2);
    K = REAL(K2);
    /* */
    SEXP ans2;
    double *ans, tmp;
    PROTECT(ans2 = allocVector(REALSXP, 1)); nprotect++;
    ans = REAL(ans2);

    /* estimate */
    ans[0]  = dnorm(x[0], 0, 1, 0);
    tmp = x[0] + (delta[0] - delta0[0]) * sqrt(r[0] * n[0]) / sig[0];
    ans[0] *= pow(pnorm(tmp, 0, 1, 1, 0), K[0] - 1.0);
    tmp = sqrt(r0[0] / r[0]) *
          (x[0] + delta[0] * sqrt(r[0] * n[0]) / sig[0] - u[0] *
           sqrt(1 + r[0] / r0[0]));
    ans[0] *= pnorm(tmp, 0, 1, 1, 0);

    /* out */

    UNPROTECT(nprotect);
    return(ans2);
}



/* C_prodsum3 based on the function prodsum of   */
/* MAMS::mams version 2.0.2 [DLC/NAM 202403-010] */

SEXP C_prodsum3(SEXP x2, SEXP l2, SEXP u2, SEXP r2, SEXP r02, SEXP r0diff2,
                SEXP Jfull2, SEXP K2, SEXP delta2, SEXP delta02,
                SEXP n2, SEXP sig2, SEXP Sigma2, SEXP SigmaJ2,
                SEXP maxpts2, SEXP releps2, SEXP abseps2, SEXP tol2){

    /* inputs */
    int *Jfull, *K, *maxpts;
    double *x, *l, *u, *r, *r0, *r0diff, *Sigma, *releps, *abseps, *tol;
    double *delta, *delta0, *n, *sig, *SigmaJ;
    Jfull = INTEGER(Jfull2);
    K = INTEGER(K2);
    maxpts = INTEGER(maxpts2);
    x = REAL(x2);
    l = REAL(l2);
    u = REAL(u2);
    r = REAL(r2);
    r0 = REAL(r02);
    r0diff = REAL(r0diff2);
    Sigma = REAL(Sigma2);
    delta = REAL(delta2);
    delta0 = REAL(delta02);
    n = REAL(n2);
    sig = REAL(sig2);
    SigmaJ = REAL(SigmaJ2);
    releps = REAL(releps2);
    abseps = REAL(abseps2);
    tol = REAL(tol2);
    /* */
    SEXP ans2, lower2, upper2, L2, U2, Jcut2, infin2;
    int j=0, jj=0, nprotect=0, *Jcut, *infin;
    double *ans, *lower, *upper, *L, *U;
    PROTECT(ans2 = allocVector(REALSXP, 1)); nprotect++;
    ans = REAL(ans2);
    PROTECT(lower2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    lower = REAL(lower2);
    PROTECT(upper2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    upper = REAL(upper2);
    PROTECT(L2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    L = REAL(L2);
    PROTECT(U2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    U = REAL(U2);
    PROTECT(infin2 = allocVector(INTSXP, Jfull[0])); nprotect++;
    infin = INTEGER(infin2);
    PROTECT(Jcut2 = allocVector(INTSXP, 1)); nprotect++;
    Jcut = INTEGER(Jcut2);
    double p0 = 1.0, p1, p2, p3, p4, insum = 0.0;

    /* make sure mvtdst has access to RNG */
    GetRNGstate();

    /* part I */

    p3 = 0.0;
    for(j = 0; j < Jfull[0]; j++){
        p0 *= dnorm(x[j], 0, 1, 0);
        p3 += sqrt(r0diff[j]) * x[j];
        p4  = delta0[0] * sqrt(r[j] * n[0]) / sig[0];
        p2  = sqrt(1 + r[j] / r0[j]);
        p1  = sqrt(r[j]) / r0[j] * p3;
        L[j] = p1 + l[j] * p2 - p4;
        U[j] = p1 + u[j] * p2 - p4;
        for(jj = 0; jj < (j+1); jj++){
            lower[jj] = L[jj];
            if(jj < j){
                upper[jj] = U[jj];
                infin[jj] = 2;
            }else{
                if(jj == (Jfull[0]-1)){
                    upper[jj] = x[jj] + (delta[0] - delta0[0]) *
                                sqrt(r[jj] * n[0]) / sig[0];
                }else{
                    upper[jj] = L[jj];
                }
                infin[jj] = 0;
            }
        }
        Jcut[0] = j+1;
        insum +=  C_MAMS_pnorm(lower, upper, infin, Sigma, Jcut, Jfull,
                               maxpts, releps, abseps, tol);
    }

    ans[0] = p0*pow(insum,K[0]-1);

    /* part II */

    p3 = 0.0;
    Jcut[0] = Jfull[0]-1;
    for(j = 0; j < Jcut[0]; j++){
        p3 += sqrt(r0diff[j]) * x[j];
        p4  = delta[0] * sqrt(r[j] * n[0]) / sig[0] +
              sqrt(r[j] / r[Jcut[0]]) * x[Jcut[0]];
        p2  = sqrt(1 + r[j] / r0[j]);
        p1  = sqrt(r[Jcut[0]] / (r[Jcut[0]] - r[j]));
        L[j] = p1 * (sqrt(r[j]) / r0[j] * p3 + l[j] * p2 - p4);
        U[j] = p1 * (sqrt(r[j]) / r0[j] * p3 + u[j] * p2 - p4);
        infin[j] = 2.0;
    }
    ans[0] *= C_MAMS_pnorm(L, U, infin, SigmaJ, Jcut, Jcut,
                           maxpts, releps, abseps, tol);
    p4 = (r0[Jcut[0]] / sqrt(r[Jcut[0]]) *
           (x[Jcut[0]] + delta[0] * sqrt(r[Jcut[0]] * n[0]) / sig[0] -
            u[Jcut[0]] * sqrt(1 + r[Jcut[0]] / r0[Jcut[0]])) - p3) /
         sqrt(r0diff[Jcut[0]]);
    ans[0] *= pnorm(p4, 0, 1, 1, 0);

    if (isnan(ans[0]) || isinf(ans[0])) {
        ans[0] = 0.0;
    }
    /* out */
    PutRNGstate();

    UNPROTECT(nprotect);
    return(ans2);
}



/* C_prodsum1_nb based on the function prodsum of  */
/* MAMS::new.bounds version 1.4.2 [DLC 202201-007] */

SEXP C_prodsum1_nb(SEXP x2, SEXP l2, SEXP u2, SEXP R2, SEXP r02, SEXP r0diff2,
                   SEXP Jfull2, SEXP K2, SEXP Sigma2,
                   SEXP maxpts2, SEXP releps2, SEXP abseps2, SEXP tol2){

    /* inputs */
    int *Jfull, *K, *maxpts;
    double *x, *l, *u, *R, *r0, *r0diff, *Sigma, *releps, *abseps, *tol;
    Jfull = INTEGER(Jfull2);
    K = INTEGER(K2);
    maxpts = INTEGER(maxpts2);
    x = REAL(x2);
    l = REAL(l2);
    u = REAL(u2);
    R = REAL(R2);
    r0 = REAL(r02);
    r0diff = REAL(r0diff2);
    Sigma = REAL(Sigma2);
    releps = REAL(releps2);
    abseps = REAL(abseps2);
    tol = REAL(tol2);
    /* */
    SEXP ans2, lower2, upper2, L2, U2, Jcut2, r2, infin2, Sigmak2;
    int k=0, j=0, jj=0, pos, nprotect=0, *Jcut, *infin;
    double *ans, *lower, *upper, *L, *U, *r, *Sigmak;
    PROTECT(ans2 = allocVector(REALSXP, 1)); nprotect++;
    ans = REAL(ans2);
    PROTECT(lower2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    lower = REAL(lower2);
    PROTECT(upper2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    upper = REAL(upper2);
    PROTECT(L2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    L = REAL(L2);
    PROTECT(U2 = allocVector(REALSXP, Jfull[0])); nprotect++;
    U = REAL(U2);
    PROTECT(r2 = allocVector(REALSXP, 1)); nprotect++;
    r = REAL(r2);
    PROTECT(Sigmak2 = allocVector(REALSXP, Jfull[0]*Jfull[0])); nprotect++;
    Sigmak = REAL(Sigmak2);
    PROTECT(infin2 = allocVector(INTSXP, Jfull[0])); nprotect++;
    infin = INTEGER(infin2);
    PROTECT(Jcut2 = allocVector(INTSXP, 1)); nprotect++;
    Jcut = INTEGER(Jcut2);
    double p1, p2, p3, insum;
    ans[0] = 1.0;

    /* make sure mvtdst has access to RNG */
    GetRNGstate();
    for(j = 0; j < Jfull[0]; j++){
        ans[0] *= dnorm(x[j], 0, 1, 0);
    }

    /* outer loop */

    for(k = 0; k < K[0]; k++){

        /* Sigma of arm of interest*/
        for(j = 0; j < Jfull[0]; j++){
            for(jj = 0; jj < (j+1); jj++){
                pos = jj * Jfull[0] + j;
                Sigmak[pos] = Sigma[pos + k * Jfull[0] * Jfull[0]];
                Sigmak[j * Jfull[0] + jj] = Sigmak[pos];
            }
        }

        /* inner loop */

        p3 = 0.0;
        insum = 0.0;
        for(j = 0; j < Jfull[0]; j++){
            r[0] = R[Jfull[0]*k + j];
            p3  += sqrt(r0diff[j]) * x[j];
            p2   = sqrt(1 + r[0] / r0[j]);
            p1   = sqrt(r[0]) / r0[j] * p3;
            L[j] = p1 + l[j] * p2;
            U[j] = p1 + u[j] * p2;
            for(jj = 0; jj < (j+1); jj++){
                lower[jj] = L[jj];
                if(jj < j){
                    upper[jj] = U[jj];
                    infin[jj] = 2;
                }else{
                    upper[jj] = L[jj];
                    infin[jj] = 0;
                }
            }
            Jcut[0] = j+1;
            insum +=  C_MAMS_pnorm(lower, upper, infin, Sigmak, Jcut, Jfull,
                                   maxpts, releps, abseps, tol);
        }
        ans[0] *= insum;

    }

    /* out */
    PutRNGstate();

    UNPROTECT(nprotect);
    return(ans2);
}


