
/**
    Some convenience functions
    *\file Convenience.c
    *\author $Author: hothorn $
    *\date $Date: 2008-10-15 11:04:00 +0200 (Wed, 15 Oct 2008) $
*/
                
#include "party.h"


/**
    Linear statistic of x, y, and weights 
    and its conditional expectation and covariance \n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param cexpcovinf logical: recompute exp and cov of the influence fct
    *\param expcovinf an object of class `ExpectCovarInfluence'
    *\param ans return value; an object of class `LinStatExpectCovar'
*/
                                                
void C_LinStatExpCov(const double *x, const int p,
                     const double *y, const int q,
                     const double *weights, const int n,
                     const int cexpcovinf, SEXP expcovinf, SEXP ans) {

    C_LinearStatistic(x, p, y, q, weights, n, 
                      REAL(GET_SLOT(ans, PL2_linearstatisticSym)));
    if (cexpcovinf)
        C_ExpectCovarInfluence(y, q, weights, n, expcovinf);
    C_ExpectCovarLinearStatistic(x, p, y, q, weights, n, 
                                 expcovinf, ans);
}


/**
    Moore-Penrose inverse of the covariance matrix \n
    *\param linexpcov an object of class `LinStatExpectCovarMPinv'
    *\param tol tolerance
*/

void C_LinStatExpCovMPinv(SEXP linexpcov, double tol) {
    C_MPinv(GET_SLOT(linexpcov, PL2_covarianceSym), tol, 
            GET_SLOT(linexpcov, PL2_svdmemSym), linexpcov);
}


/**
    Compute test statistic
    *\param linexpcov an object of class `LinStatExpectCovar'
    *\param type integer, 1 (maxabs) or 2 (quadform)
    *\param tol tolerance
*/

double C_TestStatistic(const SEXP linexpcov, const int type, const double tol) {

    int pq;
    double ans = 0.0;
    
    pq = get_dimension(linexpcov);

    switch(type) {
        /* maxabs-type test statistic */
        case 1:
            ans = C_maxabsTestStatistic(
                REAL(GET_SLOT(linexpcov, PL2_linearstatisticSym)),
                REAL(GET_SLOT(linexpcov, PL2_expectationSym)),
                REAL(GET_SLOT(linexpcov, PL2_covarianceSym)),
                pq, tol);
            break;
        /* quadform-type test statistic */
        case 2:
            ans = C_quadformTestStatistic(
                REAL(GET_SLOT(linexpcov, PL2_linearstatisticSym)), 
                REAL(GET_SLOT(linexpcov, PL2_expectationSym)),
                REAL(GET_SLOT(linexpcov, PL2_MPinvSym)), pq);
            break;
        default: error("C_TestStatistic: undefined value for type argument");
    }
    return(ans);
}


/**
    Compute asymptotic conditional P-value
    *\param tstat test statistic
    *\param linexpcov an object of class `LinStatExpectCovar'
    *\param type integer, 1 (maxabs) or 2 (quadform)
    *\param tol tolerance
    *\param maxpts argument to C_maxabsConditionalPvalue
    *\param releps argument to C_maxabsConditionalPvalue
    *\param abseps argument to C_maxabsConditionalPvalue
*/

double C_ConditionalPvalue(const double tstat, SEXP linexpcov,
                           const int type, double tol,
                           int *maxpts, double *releps, double *abseps) {
                           
    int pq;
    double ans = 1.0;
    
    pq = get_dimension(linexpcov);

    switch(type) {
        /* maxabs-type test statistic */
        case MAXABS:
            ans = C_maxabsConditionalPvalue(tstat,
                REAL(GET_SLOT(linexpcov, PL2_covarianceSym)),
                pq, maxpts, releps, abseps, &tol);
            break;
        /* quadform-type test statistic */
        case QUADFORM:
            /* var = 0 => rank = 0 */
            if (REAL(GET_SLOT(linexpcov, PL2_rankSym))[0] > 0.5)
                ans = C_quadformConditionalPvalue(tstat, 
                    REAL(GET_SLOT(linexpcov, PL2_rankSym))[0]);
            break;
        default: error("C_ConditionalPvalue: undefined value for type argument");
    }
    return(ans);
}


/**
    extract the (first) response variable from a learning sample
    *\param learnsample an object of class `LearningSample'
*/

SEXP R_get_response(SEXP learnsample) {
    return(VECTOR_ELT(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                               PL2_variablesSym), 0));
}


/**
    change the values of the response variable in a learning sample
    *\param learnsample an object of class `LearningSample'
    *\param y a REAL with new values
*/

void R_set_response(SEXP learnsample, SEXP y) {

    double *v, *t, *j, *dy, *p;
    int i, n;
    
    n = LENGTH(y);
    dy = REAL(y);
    
    if (LENGTH(R_get_response(learnsample)) != n)
        error("lengths of arguments don't match");
    
    v = REAL(VECTOR_ELT(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                                 PL2_variablesSym), 0));
    t = REAL(VECTOR_ELT(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                                 PL2_transformationsSym), 0));
    j = REAL(get_test_trafo(GET_SLOT(learnsample, PL2_responsesSym)));
    p = REAL(get_predict_trafo(GET_SLOT(learnsample, PL2_responsesSym)));
    
    for (i = 0; i < n; i++) {
        v[i] = dy[i];
        t[i] = dy[i];
        j[i] = dy[i];
        p[i] = dy[i];
    }
}
