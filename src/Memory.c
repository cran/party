
#include "party.h"

SEXP new_ExpectCovarInfluence(int q) {

    SEXP ans, expect, covar, sw;

    PROTECT(ans = party_NEW_OBJECT("ExpectCovarInfluence"));
    SET_SLOT(ans, PL2_expectationSym,
             expect = PROTECT(allocVector(REALSXP, q)));
    for (int i = 0; i < q; i++)
        REAL(expect)[i] = 0.0;
    SET_SLOT(ans, PL2_covarianceSym,
             covar = PROTECT(allocMatrix(REALSXP, q, q)));
    for (int i = 0; i < q * q; i++)
        REAL(covar)[i] = 0.0;
    SET_SLOT(ans, PL2_sumweightsSym, 
             sw = PROTECT(allocVector(REALSXP, 1)));
    REAL(sw)[0] = log(1);
    SET_SLOT(ans, PL2_dimensionSym, PROTECT(ScalarInteger(q)));
    UNPROTECT(5);
    return(ans);
}

SEXP new_LinStatExpectCovar(int p, int q) {

    SEXP ans, expect, covar, linearstatistic;

    PROTECT(ans = party_NEW_OBJECT("LinStatExpectCovar"));

    SET_SLOT(ans, PL2_expectationSym, expect = PROTECT(allocVector(REALSXP, p * q)));
    for (int i = 0; i < p * q; i++)
        REAL(expect)[i] = 0.0;

    SET_SLOT(ans, PL2_covarianceSym, covar = PROTECT(allocMatrix(REALSXP, p * q, p * q)));
    for (int i = 0; i < p * q * p * q; i++)
        REAL(covar)[i] = 0.0;

    SET_SLOT(ans, PL2_dimensionSym, PROTECT(ScalarInteger(p * q)));

    SET_SLOT(ans, PL2_linearstatisticSym, linearstatistic = PROTECT(allocVector(REALSXP, p * q)));
    for (int i = 0; i < p * q; i++)
        REAL(linearstatistic)[i] = 0.0;

    SET_SLOT(ans, PL2_expcovinfSym, PROTECT(new_ExpectCovarInfluence(q)));  

    UNPROTECT(6);
    return(ans);
}


SEXP new_svd_mem(int p) {

    SEXP ans, u, v, s;

    PROTECT(ans = party_NEW_OBJECT("svd_mem"));

    SET_SLOT(ans, PL2_pSym, PROTECT(ScalarInteger(p)));
    SET_SLOT(ans, PL2_methodSym, PROTECT(mkString("dgesdd")));
    SET_SLOT(ans, PL2_jobuSym, PROTECT(mkString("S")));
    SET_SLOT(ans, PL2_jobvSym, PROTECT(mkString("")));

    SET_SLOT(ans, PL2_uSym, u = PROTECT(allocMatrix(REALSXP, p, p)));
    for (int i = 0; i < p * p; i++)
        REAL(u)[i] = 0.0;
    SET_SLOT(ans, PL2_vSym, v = PROTECT(allocMatrix(REALSXP, p, p)));
    for (int i = 0; i < p * p; i++)
        REAL(v)[i] = 0.0;
    SET_SLOT(ans, PL2_sSym, s = PROTECT(allocVector(REALSXP, p)));
    for (int i = 0; i < p; i++)
        REAL(s)[i] = 0.0;

    UNPROTECT(8);
    return(ans);
}


SEXP new_LinStatExpectCovarMPinv(int p, int q) {

    SEXP ans, expect, covar, linearstatistic, MPinv;

    PROTECT(ans = party_NEW_OBJECT("LinStatExpectCovarMPinv"));

    SET_SLOT(ans, PL2_expectationSym, expect = PROTECT(allocVector(REALSXP, p * q)));
    for (int i = 0; i < p * q; i++)
        REAL(expect)[i] = 0.0;

    SET_SLOT(ans, PL2_covarianceSym, covar = PROTECT(allocMatrix(REALSXP, p * q, p * q)));
    for (int i = 0; i < p * q * p * q; i++)
        REAL(covar)[i] = 0.0;

    SET_SLOT(ans, PL2_dimensionSym, PROTECT(ScalarInteger(p * q)));

    SET_SLOT(ans, PL2_linearstatisticSym, linearstatistic = PROTECT(allocVector(REALSXP, p * q)));
    for (int i = 0; i < p * q; i++)
        REAL(linearstatistic)[i] = 0.0;

    SET_SLOT(ans, PL2_MPinvSym, MPinv = PROTECT(allocMatrix(REALSXP, p * q, p * q)));
    for (int i = 0; i < p * q * p * q; i++)
        REAL(MPinv)[i] = 0.0;

    SET_SLOT(ans, PL2_rankSym, PROTECT(ScalarReal(0.0)));

    SET_SLOT(ans, PL2_svdmemSym, PROTECT(new_svd_mem(p * q)));  

    SET_SLOT(ans, PL2_expcovinfSym, PROTECT(new_ExpectCovarInfluence(q)));  

    UNPROTECT(9);
    return(ans);
}


SEXP ctree_memory (SEXP object, SEXP MP_INV) {

    SEXP ans, weights, splitstatistics, dontuse, dontusetmp, varmemory;
    int q, p, nobs, ninputs;

    q = ncol(get_test_trafo(GET_SLOT(object, PL2_responsesSym)));
    
    ninputs = get_ninputs(object);
    nobs = get_nobs(object);

    ans = PROTECT(party_NEW_OBJECT("TreeFitMemory"));
    SET_SLOT(ans, PL2_expcovinfSym, PROTECT(new_ExpectCovarInfluence(q)));
    SET_SLOT(ans, PL2_expcovinfssSym, PROTECT(new_ExpectCovarInfluence(1)));
    SET_SLOT(ans, PL2_linexpcov2sampleSym, PROTECT(new_LinStatExpectCovar(1, q)));

    SET_SLOT(ans, PL2_weightsSym, weights = PROTECT(allocVector(REALSXP, nobs)));
    for (int i = 0; i < nobs; i++)
        REAL(weights)[i] = 0.0;
    SET_SLOT(ans, PL2_splitstatisticsSym, splitstatistics = PROTECT(allocVector(REALSXP, nobs)));
    for (int i = 0; i < nobs; i++)
        REAL(splitstatistics)[i] = 0.0;
    SET_SLOT(ans, PL2_dontuseSym, dontuse = PROTECT(allocVector(LGLSXP, ninputs)));
    for (int i = 0; i < ninputs; i++)
        LOGICAL(dontuse)[i] = 0.0;
    SET_SLOT(ans, PL2_dontusetmpSym, dontusetmp = PROTECT(allocVector(LGLSXP, ninputs)));
    for (int i = 0; i < ninputs; i++)
        LOGICAL(dontusetmp)[i] = 0.0;


    varmemory = PROTECT(allocVector(VECSXP, ninputs));
    
    for (int i = 0; i < ninputs; i++) {
    
        p = ncol(get_transformation(GET_SLOT(object, PL2_inputsSym), i + 1));

        if (LOGICAL(MP_INV)[0]) {        
            SET_VECTOR_ELT(varmemory, i, new_LinStatExpectCovarMPinv(p, q));
        } else {
           SET_VECTOR_ELT(varmemory, i, new_LinStatExpectCovar(p, q));
        }
    }

    SET_SLOT(ans, PL2_varmemorySym, varmemory);

    UNPROTECT(9);
    return(ans);
}

