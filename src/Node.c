
/**
    Node computations
    *\file Node.c
    *\author $Author: hothorn $
    *\date $Date: 2007-01-15 11:24:41 +0100 (Mon, 15 Jan 2007) $
*/
                
#include "party.h"


/**
    Compute prediction of a node
    *\param y the response variable (raw numeric values or dummy encoded factor)
    *\param n number of observations
    *\param q number of columns of y
    *\param weights case weights
    *\param sweights sum of case weights
    *\param ans return value; the q-dimensional predictions
*/
        
void C_prediction(const double *y, int n, int q, const double *weights, 
                  const double sweights, double *ans) {

    int i, j, jn;
    
    for (j = 0; j < q; j++) {
        ans[j] = 0.0;
        jn = j * n;
        for (i = 0; i < n; i++) 
            ans[j] += weights[i] * y[jn + i];
        ans[j] = ans[j] / sweights;
    }
}


/**
    The main function for all node computations
    *\param node an initialized node (an S3 object!)
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param TERMINAL logical indicating if this node will
                     be a terminal node
*/

void C_Node(SEXP node, SEXP learnsample, SEXP weights, 
            SEXP fitmem, SEXP controls, int TERMINAL) {
    
    int nobs, ninputs, jselect, q, j, k, i;
    double mincriterion, sweights, *dprediction;
    double *teststat, *pvalue, smax, cutpoint = 0.0, maxstat = 0.0;
    double *standstat, *splitstat;
    SEXP responses, inputs, x, expcovinf, thisweights, linexpcov;
    SEXP varctrl, splitctrl, gtctrl, tgctrl, split, jointy;
    double *dxtransf, *dweights;
    int *itable;
    
    nobs = get_nobs(learnsample);
    ninputs = get_ninputs(learnsample);
    varctrl = get_varctrl(controls);
    splitctrl = get_splitctrl(controls);
    gtctrl = get_gtctrl(controls);
    tgctrl = get_tgctrl(controls);
    mincriterion = get_mincriterion(gtctrl);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    jointy = get_jointtransf(responses);
    q = ncol(jointy);

    /* <FIXME> we compute C_GlobalTest even for TERMINAL nodes! </FIXME> */

    /* compute the test statistics and the node criteria for each input */        
    C_GlobalTest(learnsample, weights, fitmem, varctrl,
                 gtctrl, get_minsplit(splitctrl), 
                 REAL(S3get_teststat(node)), REAL(S3get_criterion(node)));
    
    /* sum of weights: C_GlobalTest did nothing if sweights < mincriterion */
    sweights = REAL(GET_SLOT(GET_SLOT(fitmem, PL2_expcovinfSym), 
                             PL2_sumweightsSym))[0];

    /* compute the prediction of this node */
    dprediction = REAL(S3get_prediction(node));

    /* <FIXME> feed raw numeric values OR dummy encoded factors as y 
       Problem: what happens for survival times ? */
    C_prediction(REAL(jointy), nobs, q, REAL(weights), 
                     sweights, dprediction);
    /* </FIXME> */

    teststat = REAL(S3get_teststat(node));
    pvalue = REAL(S3get_criterion(node));

    /* try the two out of ninputs best inputs variables */
    /* <FIXME> be more flexible and add a parameter controlling
               the number of inputs tried </FIXME> */
    for (j = 0; j < 2; j++) {

        smax = C_max(pvalue, ninputs);
        REAL(S3get_maxcriterion(node))[0] = smax;
    
        /* if the global null hypothesis was rejected */
        if (smax > mincriterion && !TERMINAL) {

            /* the input variable with largest association to the response */
            jselect = C_whichmax(pvalue, teststat, ninputs) + 1;

            /* get the raw numeric values or the codings of a factor */
            x = get_variable(inputs, jselect);
            if (has_missings(inputs, jselect)) {
                expcovinf = GET_SLOT(get_varmemory(fitmem, jselect), 
                                    PL2_expcovinfSym);
                thisweights = get_weights(fitmem, jselect);
            } else {
                expcovinf = GET_SLOT(fitmem, PL2_expcovinfSym);
                thisweights = weights;
            }

            /* <FIXME> handle ordered factors separatly??? </FIXME> */
            if (!is_nominal(inputs, jselect)) {
            
                /* search for a split in a ordered variable x */
                split = S3get_primarysplit(node);
                
                /* check if the n-vector of splitstatistics 
                   should be returned for each primary split */
                if (get_savesplitstats(tgctrl)) {
                    C_init_orderedsplit(split, nobs);
                    splitstat = REAL(S3get_splitstatistics(split));
                } else {
                    C_init_orderedsplit(split, 0);
                    splitstat = REAL(get_splitstatistics(fitmem));
                }

                C_split(REAL(x), 1, REAL(jointy), q, REAL(weights), nobs,
                        INTEGER(get_ordering(inputs, jselect)), splitctrl, 
                        GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                        expcovinf, REAL(S3get_splitpoint(split)), &maxstat,
                        splitstat);
                S3set_variableID(split, jselect);
             } else {
           
                 /* search of a set of levels (split) in a numeric variable x */
                 split = S3get_primarysplit(node);
                 
                /* check if the n-vector of splitstatistics 
                   should be returned for each primary split */
                if (get_savesplitstats(tgctrl)) {
                    C_init_nominalsplit(split, 
                        LENGTH(get_levels(inputs, jselect)), 
                        nobs);
                    splitstat = REAL(S3get_splitstatistics(split));
                } else {
                    C_init_nominalsplit(split, 
                        LENGTH(get_levels(inputs, jselect)), 
                        0);
                    splitstat = REAL(get_splitstatistics(fitmem));
                }
          
                 linexpcov = get_varmemory(fitmem, jselect);
                 standstat = Calloc(get_dimension(linexpcov), double);
                 C_standardize(REAL(GET_SLOT(linexpcov, 
                                             PL2_linearstatisticSym)),
                               REAL(GET_SLOT(linexpcov, PL2_expectationSym)),
                               REAL(GET_SLOT(linexpcov, PL2_covarianceSym)),
                               get_dimension(linexpcov), get_tol(splitctrl), 
                               standstat);
 
                 C_splitcategorical(INTEGER(x), 
                                    LENGTH(get_levels(inputs, jselect)), 
                                    REAL(jointy), q, REAL(weights), 
                                    nobs, standstat, splitctrl, 
                                    GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                                    expcovinf, &cutpoint, 
                                    INTEGER(S3get_splitpoint(split)),
                                    &maxstat, splitstat);

                 /* compute which levels of a factor are available in this node 
                    (for printing) later on. A real `table' for this node would
                    induce too much overhead here. Maybe later. */
                    
                 itable = INTEGER(S3get_table(split));
                 dxtransf = REAL(get_transformation(inputs, jselect));
                 dweights = REAL(thisweights);
                 for (k = 0; k < LENGTH(get_levels(inputs, jselect)); k++) {
                     itable[k] = 0;
                     for (i = 0; i < nobs; i++) {
                         if (dxtransf[k * nobs + i] * dweights[i] > 0) {
                             itable[k] = 1;
                             continue;
                         }
                     }
                 }

                 Free(standstat);
            }
            if (maxstat == 0) {
                warning("no admissible split found\n");
            
                if (j == 1) {          
                    S3set_nodeterminal(node);
                } else {
                    /* <FIXME> why? </FIXME> */
                    pvalue[jselect - 1] = 0.0;
                }
            } else {
                S3set_variableID(split, jselect);
                break;
            }
        } else {
            S3set_nodeterminal(node);
            break;
        }
    }
}       


/**
    R-interface to C_Node
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
*/

SEXP R_Node(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls) {
            
     SEXP ans;
     
     PROTECT(ans = allocVector(VECSXP, NODE_LENGTH));
     C_init_node(ans, get_nobs(learnsample), get_ninputs(learnsample), 
                 get_maxsurrogate(get_splitctrl(controls)),
                 ncol(get_jointtransf(GET_SLOT(learnsample, PL2_responsesSym))));

     C_Node(ans, learnsample, weights, fitmem, controls, 0);
     UNPROTECT(1);
     return(ans);
}
