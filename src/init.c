
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP party_init();
extern SEXP R_Ensemble(SEXP, SEXP, SEXP);
extern SEXP R_Ensemble_weights(SEXP, SEXP, SEXP);
extern SEXP R_ExpectCovarInfluence(SEXP, SEXP);
extern SEXP R_ExpectCovarLinearStatistic(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_get_nodebynum(SEXP, SEXP);
extern SEXP R_get_nodeID(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_getpredictions(SEXP, SEXP);
extern SEXP R_LinearStatistic(SEXP, SEXP, SEXP);
extern SEXP R_maxabsTestStatistic(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_MPinv(SEXP, SEXP, SEXP);
extern SEXP R_PermutedLinearStatistic(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_predict(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_predictRF_weights(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_proximity(SEXP);
extern SEXP R_quadformTestStatistic(SEXP, SEXP, SEXP);
extern SEXP R_split(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_splitcategorical(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_svd(SEXP, SEXP);
extern SEXP R_TreeGrow(SEXP, SEXP, SEXP);
extern SEXP R_quadformConditionalPvalue(SEXP, SEXP);
extern SEXP R_maxabsConditionalPvalue(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_IndependenceTest(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_Node(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_kronecker(SEXP, SEXP);
extern SEXP R_max(SEXP);
extern SEXP R_abs(SEXP);
extern SEXP R_matprod(SEXP, SEXP);
extern SEXP R_matprodT(SEXP, SEXP);
extern SEXP R_permute(SEXP);
extern SEXP R_rsubset(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"party_init",                   (DL_FUNC) &party_init,                   0},
    {"R_Ensemble",                   (DL_FUNC) &R_Ensemble,                   3},
    {"R_Ensemble_weights",           (DL_FUNC) &R_Ensemble_weights,           3},
    {"R_ExpectCovarInfluence",       (DL_FUNC) &R_ExpectCovarInfluence,       2},
    {"R_ExpectCovarLinearStatistic", (DL_FUNC) &R_ExpectCovarLinearStatistic, 4},
    {"R_get_nodebynum",              (DL_FUNC) &R_get_nodebynum,              2},
    {"R_get_nodeID",                 (DL_FUNC) &R_get_nodeID,                 4},
    {"R_getpredictions",             (DL_FUNC) &R_getpredictions,             2},
    {"R_LinearStatistic",            (DL_FUNC) &R_LinearStatistic,            3},
    {"R_maxabsTestStatistic",        (DL_FUNC) &R_maxabsTestStatistic,        4},
    {"R_MPinv",                      (DL_FUNC) &R_MPinv,                      3},
    {"R_PermutedLinearStatistic",    (DL_FUNC) &R_PermutedLinearStatistic,    4},
    {"R_predict",                    (DL_FUNC) &R_predict,                    4},
    {"R_predictRF_weights",          (DL_FUNC) &R_predictRF_weights,          6},
    {"R_proximity",                  (DL_FUNC) &R_proximity,                  1},
    {"R_quadformTestStatistic",      (DL_FUNC) &R_quadformTestStatistic,      3},
    {"R_split",                      (DL_FUNC) &R_split,                      7},
    {"R_splitcategorical",           (DL_FUNC) &R_splitcategorical,           8},
    {"R_svd",                        (DL_FUNC) &R_svd,                        2},
    {"R_TreeGrow",                   (DL_FUNC) &R_TreeGrow,                   3},
    {"R_quadformConditionalPvalue",  (DL_FUNC) &R_quadformConditionalPvalue,  2},
    {"R_maxabsConditionalPvalue",    (DL_FUNC) &R_maxabsConditionalPvalue,    6},
    {"R_IndependenceTest",           (DL_FUNC) &R_IndependenceTest,           5},
    {"R_Node",                       (DL_FUNC) &R_Node,                       4},
    {"R_kronecker",                  (DL_FUNC) &R_kronecker,                  2},
    {"R_abs",                        (DL_FUNC) &R_abs,                        1},
    {"R_max",                        (DL_FUNC) &R_max,                        1},
    {"R_matprod",                    (DL_FUNC) &R_matprod,                    2},
    {"R_matprodT",                   (DL_FUNC) &R_matprodT,                   2},
    {"R_permute",                    (DL_FUNC) &R_permute,                    1},
    {"R_rsubset",                    (DL_FUNC) &R_rsubset,                    2},
    {NULL, NULL, 0}
};

void R_init_party(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
