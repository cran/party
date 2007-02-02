
/**
    Some commonly needed utility functions.
    *\file Utils.c
    *\author $Author: hothorn $
    *\date $Date: 2007-02-02 11:22:45 +0100 (Fri, 02 Feb 2007) $
*/
                
#include "party.h"
                
                
/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans) {

    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
                }
            }
        }
    }
}  


/**
    R-interface to C_kronecker
    *\param A matrix
    *\param B matrix
*/
                
SEXP R_kronecker (SEXP A, SEXP B) {

    /*  The Kronecker product, a real (mr x ns) matrix */
    SEXP ans; 
    int *adim, *bdim;

    if (!isReal(A) || !isReal(B)) 
        error("R_kronecker: A and B are not of type REALSXP");

    if (isMatrix(A)) {
        adim = INTEGER(getAttrib(A, R_DimSymbol));
    } else {
        /* assume row vectors */
        adim = Calloc(2, int);
        adim[0] = 1;
        adim[1] = LENGTH(A);
    }
    
    if (isMatrix(B)) {
        bdim = INTEGER(getAttrib(B, R_DimSymbol));
    } else {
        /* assume row vectors */
        bdim = Calloc(2, int);
        bdim[0] = 1;
        bdim[1] = LENGTH(B);
    }

    PROTECT(ans = allocMatrix(REALSXP, 
                              adim[0] * bdim[0], 
                              adim[1] * bdim[1]));
    C_kronecker(REAL(A), adim[0], adim[1], 
                REAL(B), bdim[0], bdim[1], REAL(ans));
    if (!isMatrix(A)) Free(adim); 
    if (!isMatrix(B)) Free(bdim);
    UNPROTECT(1);
    return(ans);
}


/**
    C- and R-interface to La_svd (R/src/main/lapack.c)
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

SEXP CR_svd (SEXP x, SEXP svdmem) {

    int p, i;
    double *du, *dv;

    if (!isMatrix(x) || !isReal(x))
        error("x is not a real matrix");

    du = REAL(GET_SLOT(svdmem, PL2_uSym));
    dv = REAL(GET_SLOT(svdmem, PL2_vSym));
    p = INTEGER(GET_SLOT(svdmem, PL2_pSym))[0];
    for (i = 0; i < p*p; i++) {
        du[i] = 0.0;
        dv[i] = 0.0;
    }
    SET_SLOT(svdmem, PL2_svdSym, La_svd(GET_SLOT(svdmem, PL2_jobuSym), 
        GET_SLOT(svdmem, PL2_jobvSym), x, GET_SLOT(svdmem, PL2_sSym), 
        GET_SLOT(svdmem, PL2_uSym), GET_SLOT(svdmem, PL2_vSym), 
        GET_SLOT(svdmem, PL2_methodSym)));
    return(R_NilValue);
}


/**
    Moore-Penrose inverse of a matrix
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
    *\param ans return value; an object of class `ExpectCovarMPinv'
*/

void C_MPinv (SEXP x, double tol, SEXP svdmem, SEXP ans) {

    SEXP svdx, d, u, vt, dummy;
    int i, j, p, k, *positive;
    double *dd, *du, *dvt, *dMPinv;
    double *drank;
    
    drank = REAL(GET_SLOT(ans, PL2_rankSym));
    dMPinv = REAL(GET_SLOT(ans, PL2_MPinvSym));

    dummy = CR_svd(x, svdmem);
    svdx = GET_SLOT(svdmem, PL2_svdSym);
    d = VECTOR_ELT(svdx, 0);
    dd = REAL(d);
    u = VECTOR_ELT(svdx, 1);
    du = REAL(u);
    vt = VECTOR_ELT(svdx, 2);
    dvt = REAL(vt);
    p = LENGTH(d);

    if (tol * dd[0] > tol) tol = tol * dd[0];

    positive = Calloc(p, int); 
    
    drank[0] = 0.0;
    for (i = 0; i < p; i++) {
        if (dd[i] > tol) {
            positive[i] = 1;
            drank[0] += 1.0;
        } 
    }
    
    for (j = 0; j < p; j++) {
        if (positive[j]) {
            for (i = 0; i < p; i++)
                du[j * p + i] *= (1 / dd[j]);
        }
    }
    
    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++) {
            dMPinv[j * p + i] = 0.0;
            for (k = 0; k < p; k++) {
                if (positive[k])
                    dMPinv[j * p + i] += dvt[i * p + k] * du[p * k + j]; 
            }
        }
    }

    Free(positive);
}

/**
    R-interface to C_MPinv 
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_MPinv (SEXP x, SEXP tol, SEXP svdmem) {

    SEXP ans;
    int p;

    if (!isMatrix(x) || !isReal(x))
        error("R_MPinv: x is not a real matrix");

    if (nrow(x) != ncol(x)) 
        error("R_MPinv: x is not a square matrix");

    if (!isReal(tol) || LENGTH(tol) != 1)
        error("R_MPinv: tol is not a scalar real");
    
    p = nrow(x);
    if (p != INTEGER(GET_SLOT(svdmem, PL2_pSym))[0])
        error("R_MPinv: dimensions don't match");

    PROTECT(ans = NEW_OBJECT(MAKE_CLASS("LinStatExpectCovarMPinv")));
    SET_SLOT(ans, PL2_MPinvSym, PROTECT(allocMatrix(REALSXP, p, p)));
    SET_SLOT(ans, PL2_rankSym, PROTECT(allocVector(REALSXP, 1)));
    
    C_MPinv(x, REAL(tol)[0], svdmem, ans);
    
    UNPROTECT(3);
    return(ans);
}

/**
    the maximum of a double vector
    *\param x vector
    *\param n its length
*/


double C_max(const double *x, const int n) {
   double tmp = 0.0;
   int i;
   
   for (i = 0; i < n; i++) {
       if (x[i] > tmp) tmp = x[i];
   }
   return(tmp);
}


/**
    R-interface to C_max
    *\param x numeric vector
*/

SEXP R_max(SEXP x) {

    SEXP ans;
    int n;
    
    if (!isReal(x)) 
        error("R_max: x is not of type REALSXP");
    n = LENGTH(x);
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = C_max(REAL(x), n);
    UNPROTECT(1);
    return(ans);
}


/**
    absolute value 
    *\param x numeric vector
    *\param n length(x)
*/

void C_abs(double *x, int n) {

    int i;
    for (i = 0; i < n; i++) x[i] = fabs(x[i]);
}


/**
    R-interface to C_abs
    *\param x numeric vector
*/

SEXP R_abs(SEXP x) {

    SEXP ans;
    int n;
    
    if (!isReal(x)) 
        error("R_max: x is not of type REALSXP");
    n = LENGTH(x);
    PROTECT(ans = duplicate(x));
    C_abs(REAL(ans), n);
    UNPROTECT(1);
    return(ans);
}


/**
    matrix product x %*% y
    *\param x a matrix
    *\param nrx number of rows of x
    *\param ncx number of cols of x
    *\param y a matrix
    *\param nry number of rows of y
    *\param ncy number of cols of y
    *\param z a matrix of dimension nrx x ncy
*/

void C_matprod(double *x, int nrx, int ncx,
               double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    int i;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
	                x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}


/**
    R-interface to C_matprod
    *\param x a matrix
    *\param y a matrix
*/

SEXP R_matprod(SEXP x, SEXP y) {

    SEXP ans;
    
    int nrx, ncx, nry, ncy;
    
    nrx = nrow(x);
    ncx = ncol(x);
    nry = nrow(y);
    ncy = ncol(y);

    if (ncx != nry)
        error("R_matprod: dimensions don't match");
    PROTECT(ans = allocMatrix(REALSXP, nrx, ncy));
    C_matprod(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/**
    matrix product x %*% t(y)
    *\param x a matrix
    *\param nrx number of rows of x
    *\param ncx number of cols of x
    *\param y a matrix
    *\param nry number of rows of y
    *\param ncy number of cols of y
    *\param z a matrix of dimension nrx x ncy
*/

void C_matprodT(double *x, int nrx, int ncx,
                double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    int i;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncy, &one,
	                x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
}


/**
    R-interface to C_matprodT
    *\param x a matrix
    *\param y a matrix
*/

SEXP R_matprodT(SEXP x, SEXP y) {

    SEXP ans;
    int nrx, ncx, nry, ncy;
    
    nrx = nrow(x);
    ncx = ncol(x);
    nry = nrow(y);
    ncy = ncol(y);

    if (ncx != ncy)
        error("R_matprod: dimensions don't match");
    PROTECT(ans = allocMatrix(REALSXP, nrx, nry));
    C_matprodT(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans));
    UNPROTECT(1);
    return(ans);
}


/**
    compute a permutation of a (random subset of) 0:(m-1)
    *\param x an integer vector of length m
    *\param m integer
    *\param k integer
    *\param ans an integer vector of length k
*/    

void C_SampleNoReplace(int *x, int m, int k, int *ans) {
     
    int i, j, n = m;

    for (i = 0; i < m; i++)
        x[i] = i;
    for (i = 0; i < k; i++) {
        j = n * unif_rand(); 
        ans[i] = x[j];
        x[j] = x[--n];  
    }
}


/**
    R-interface to C_SampleNoReplace: the permutation case
    *\param m integer
*/    

SEXP R_permute(SEXP m) {
    
    SEXP x, ans;
    int n;
    
    n = INTEGER(m)[0];
    PROTECT(x = allocVector(INTSXP, n));
    PROTECT(ans = allocVector(INTSXP, n));
    C_SampleNoReplace(INTEGER(x), n, n, INTEGER(ans));
    UNPROTECT(2);
    return(ans);
}


/**
    R-interface to C_SampleNoReplace: the subset case
    *\param m integer
    *\param k integer
*/    

SEXP R_rsubset(SEXP m, SEXP k) {
    
    SEXP x, ans;
    int n, j;
    
    n = INTEGER(m)[0];
    j = INTEGER(k)[0];
    PROTECT(x = allocVector(INTSXP, n));
    PROTECT(ans = allocVector(INTSXP, j));
    C_SampleNoReplace(INTEGER(x), n, j, INTEGER(ans));
    UNPROTECT(2);
    return(ans);
}


/**
    determine if i is element of the integer vector set
    *\param i an integer
    *\param iset a pointer to an integer vector
    *\param p length(iset)
*/

int i_in_set(int i, int *iset, int p) {

    int j, is = 0;
        
    if (p == 0) return(0);
                    
    for (j = 0; j < p; j++) {
        if (iset[j] == i) {  
            is = 1;
            break; 
        }
    }
    return(is);
}

int C_i_in_set(int i, SEXP set) {
    if (LENGTH(set) > 0)
        return(i_in_set(i, INTEGER(set), LENGTH(set)));
    else 
        return(0);
}
    
int nrow(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
}

int ncol(SEXP x) {
    return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
}

/* compute index of variable with smallest p-value 
   (and largest test statistic in case two or more p-values coincide -- 
    should not happen anymore since we use 1 - (1 - p)^k for Bonferroni adjustment)
*/
int C_whichmax(double *pvalue, double *teststat, int ninputs) {

    int ans = -1, j;
    double tmppval = 0.0, tmptstat = 0.0;
       
    /* <FIXME> can we switch to the log scale here? </FIXME> */

    tmppval = 0.0;
    tmptstat = 0.0;
    for (j = 0; j < ninputs; j++) {
        if (pvalue[j] > tmppval) {
            ans = j;
            tmppval = pvalue[j];
            tmptstat = teststat[j];
        } else {
            if (pvalue[j] == tmppval && teststat[j] > tmptstat) {  
                ans = j;
                tmppval = pvalue[j];
                tmptstat = teststat[j];
            }
        }
    }
    return(ans);
}

SEXP R_whichmax(SEXP x, SEXP y) {
    SEXP ans;
    
    if (LENGTH(x) != LENGTH(y)) error("different length");
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = C_whichmax(REAL(x), REAL(y), LENGTH(x));
    UNPROTECT(1);
    return(ans);
}

SEXP R_listplus(SEXP a, SEXP b, SEXP which) {

    int na, nb, i, j, *iwhich;
    double *dae, *dbe;
    SEXP ae, be;

    na = LENGTH(a);
    nb = LENGTH(b);
    if (na != nb) error("a and b are of different length");
    
    iwhich = LOGICAL(which);
    
    for (i = 0; i < na; i++) {
        if (iwhich[i]) continue;
        
        ae = VECTOR_ELT(a, i);
        be = VECTOR_ELT(b, i);

        if (LENGTH(ae) != LENGTH(be)) 
            error("elements %d are of different length", i);
            
        if (!isReal(ae) || !isReal(be))
            error("elements %d are not of type double", i);
            
        dae = REAL(ae);
        dbe = REAL(be);
        for (j = 0; j < LENGTH(ae); j++) 
            dae[j] += dbe[j];
    }
    return(a);
}

SEXP R_modify_response(SEXP x, SEXP vf) {

    double *src, *tar;
    int i, n;
    
    src = REAL(x);
    n = LENGTH(x);

    tar = REAL(get_transformation(vf, 1));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_test_trafo(vf));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_predict_trafo(vf));
    for (i = 0; i < n; i++)
        tar[i] = src[i];

    tar = REAL(get_variable(vf, 1));
    for (i = 0; i < n; i++)
        tar[i] = src[i];
                                          
    return(R_NilValue);
}

double F77_SUB(unifrnd)(void) { return unif_rand(); }
