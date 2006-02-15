
/**
    Random forest with conditional inference trees
    *\file RandomForest.c
    *\author $Author: hothorn $
    *\date $Date: 2006-02-14 11:35:22 +0100 (Di, 14 Feb 2006) $
*/

#include "party.h"


/**
    An experimental implementation of random forest like algorithms \n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param ntrees number of trees grown
*/


SEXP R_Ensemble(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls, SEXP ntrees) {
            
     SEXP nweights, tree, where, ans;
     double *dnweights, *dweights, sw = 0.0, *prob;
     int nobs, i, b, B , nodenum = 1, *iweights, *iwhere;
     
     B = INTEGER(ntrees)[0];
     nobs = get_nobs(learnsample);
     PROTECT(ans = allocVector(VECSXP, B));

     iweights = Calloc(nobs, int);
     prob = Calloc(nobs, double);
     dweights = REAL(weights);

     for (i = 0; i < nobs; i++)
         sw += dweights[i];
     for (i = 0; i < nobs; i++)
         prob[i] = dweights[i]/sw;

     for (b  = 0; b < B; b++) {
         SET_VECTOR_ELT(ans, b, tree = allocVector(VECSXP, NODE_LENGTH + 1));
         SET_VECTOR_ELT(tree, NODE_LENGTH, where = allocVector(INTSXP, nobs));
         iwhere = INTEGER(where);
         for (i = 0; i < nobs; i++) iwhere[i] = 0;
     
         C_init_node(tree, nobs, get_ninputs(learnsample), 
                     get_maxsurrogate(get_splitctrl(controls)),
                     ncol(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                          PL2_jointtransfSym)));

         /* weights for a bootstrap sample */
         GetRNGstate();
         rmultinom((int) sw, prob, nobs, iweights);
         PutRNGstate();

         nweights = S3get_nodeweights(tree);
         dnweights = REAL(nweights);
         for (i = 0; i < nobs; i++) dnweights[i] = (double) iweights[i];
     
         C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum, 1);
         nodenum = 1;
     }
     Free(prob); Free(iweights);
     UNPROTECT(1);
     return(ans);
}
