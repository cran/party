
/**
    Random forest with conditional inference trees
    *\file RandomForest.c
    *\author $Author: thothorn $
    *\date $Date: 2017-02-06 17:53:27 +0100 (Mon, 06 Feb 2017) $
*/

#include "party.h"

/**
    An experimental implementation of random forest like algorithms \n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param bwhere integer matrix (n x ntree) for terminal node numbers
    *\param bweights double matrix (n x ntree) for bootstrap case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
*/


SEXP R_Ensemble(SEXP learnsample, SEXP weights, SEXP controls) {
            
     SEXP ans, nweights, tree, where, ensemble, bw, fitmem, bwhere, bweights;
     double *dnweights, *dweights, sw = 0.0, *prob, tmp;
     int nobs, i, b, B , nodenum = 1, *iweights, *iweightstmp, 
         *iwhere, replace, fraction, wgrzero = 0, realweights = 0;
     int j, k, l, swi = 0;
     
     B = get_ntree(controls);
     nobs = get_nobs(learnsample);
     
     PROTECT(ans = NEW_OBJECT(MAKE_CLASS("RandomForest")));
     PROTECT(bwhere = allocVector(VECSXP, B));
     PROTECT(bweights = allocVector(VECSXP, B));
     PROTECT(ensemble = allocVector(VECSXP, B));
     PROTECT(fitmem = ctree_memory(learnsample, PROTECT(ScalarLogical(1))));

     SET_SLOT(ans, PL2_ensembleSym, ensemble);
     SET_SLOT(ans, PL2_whereSym, bwhere);
     SET_SLOT(ans, PL2_weightsSym, bweights);

     iweights = Calloc(nobs, int);
     iweightstmp = Calloc(nobs, int);
     prob = Calloc(nobs, double);
     dweights = REAL(weights);

     for (i = 0; i < nobs; i++) {
         /* sum of weights */
         sw += dweights[i];
         /* number of weights > 0 */
         if (dweights[i] > 0) wgrzero++;
         /* case weights or real weights? */
         if (dweights[i] - ftrunc(dweights[i]) > 0) 
             realweights = 1;
     }
     for (i = 0; i < nobs; i++)
         prob[i] = dweights[i]/sw;
     swi = (int) ftrunc(sw);

     replace = get_replace(controls);
     /* fraction of number of obs with weight > 0 */
     if (realweights) {
         /* fraction of number of obs with weight > 0 for real weights*/
         tmp = (get_fraction(controls) * wgrzero);
     } else {
         /* fraction of sum of weights for case weights */
         tmp = (get_fraction(controls) * sw);
     }
     fraction = (int) ftrunc(tmp);
     if (ftrunc(tmp) < tmp) fraction++;

     if (!replace) {
         if (fraction < 10)
             error("fraction of %f is too small", fraction);
     }

     /* <FIXME> can we call those guys ONCE? what about the deeper
         calls??? </FIXME> */
     GetRNGstate();
  
     if (get_trace(controls))
         Rprintf("\n");
     for (b  = 0; b < B; b++) {
         SET_VECTOR_ELT(ensemble, b, tree = allocVector(VECSXP, NODE_LENGTH + 1));
         SET_VECTOR_ELT(bwhere, b, where = allocVector(INTSXP, nobs));
         SET_VECTOR_ELT(bweights, b, bw = allocVector(REALSXP, nobs));
         
         iwhere = INTEGER(where);
         for (i = 0; i < nobs; i++) iwhere[i] = 0;
     
         C_init_node(tree, nobs, get_ninputs(learnsample), 
                     get_maxsurrogate(get_splitctrl(controls)),
                     ncol(get_predict_trafo(GET_SLOT(learnsample, 
                                                   PL2_responsesSym))));

         /* generate altered weights for perturbation */
         if (replace) {
             /* weights for a bootstrap sample */
             rmultinom(swi, prob, nobs, iweights);
         } else {
             /* weights for sample splitting */
             C_SampleSplitting(nobs, prob, iweights, fraction);
         }

         nweights = S3get_nodeweights(tree);
         dnweights = REAL(nweights);
         for (i = 0; i < nobs; i++) {
             REAL(bw)[i] = (double) iweights[i];
             dnweights[i] = REAL(bw)[i];
         }
     
         C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum, 1);
         nodenum = 1;
         C_remove_weights(tree, 0);

         /* compute terminal node ids for all observations (not just those with weight > 0) */
         for (i = 0; i < nobs; i++)
            iwhere[i] = C_get_nodeID(tree, GET_SLOT(learnsample, PL2_inputsSym), 0.0, i, -1);
         
         if (get_trace(controls)) {
             /* progress bar; inspired by 
             http://avinashjoshi.co.in/2009/10/13/creating-a-progress-bar-in-c/ */
             Rprintf("[");
             /* Print the = until the current percentage */
             l = (int) ceil( ((double) b * 50.0) / B);
             for (j = 0; j < l; j++)
                 Rprintf("=");
             Rprintf(">");
             for (k = j; k < 50; k++)
                 Rprintf(" ");
             Rprintf("]");
             /* % completed */
                 Rprintf(" %3d%% completed", j * 2);
             /* To delete the previous line */
             Rprintf("\r");
             /* Flush all char in buffer */
             /* fflush(stdout); */
         }
     }
     if (get_trace(controls))
         Rprintf("\n");

     PutRNGstate();

     Free(prob); Free(iweights); Free(iweightstmp);
     UNPROTECT(6);
     return(ans);
}

/**
    An experimental implementation of random forest like algorithms \n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param bwhere integer matrix (n x ntree) for terminal node numbers
    *\param bweights double matrix (n x ntree) for bootstrap case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
*/


SEXP R_Ensemble_weights(SEXP learnsample, SEXP bweights, 
                        SEXP controls) {
            
     SEXP nweights, tree, where, ensemble, fitmem, ans, bwhere;
     double *dnweights, *dweights;
     int nobs, i, b, B , nodenum = 1, *iwhere;
     int j, k, l;
     
     B = get_ntree(controls);
     nobs = get_nobs(learnsample);
     
     PROTECT(ans = NEW_OBJECT(MAKE_CLASS("RandomForest")));
     PROTECT(bwhere = allocVector(VECSXP, B));
     PROTECT(ensemble = allocVector(VECSXP, B));
     PROTECT(fitmem = ctree_memory(learnsample, PROTECT(ScalarLogical(1))));

     SET_SLOT(ans, PL2_ensembleSym, ensemble);
     SET_SLOT(ans, PL2_whereSym, bwhere);
     SET_SLOT(ans, PL2_weightsSym, bweights);

     /* <FIXME> can we call those guys ONCE? what about the deeper
         calls??? </FIXME> */
     GetRNGstate();
  
     if (get_trace(controls))
         Rprintf("\n");
     for (b  = 0; b < B; b++) {
         SET_VECTOR_ELT(ensemble, b, tree = allocVector(VECSXP, NODE_LENGTH + 1));
         SET_VECTOR_ELT(bwhere, b, where = allocVector(INTSXP, nobs));
         
         iwhere = INTEGER(where);
         for (i = 0; i < nobs; i++) iwhere[i] = 0;
     
         C_init_node(tree, nobs, get_ninputs(learnsample), 
                     get_maxsurrogate(get_splitctrl(controls)),
                     ncol(get_predict_trafo(GET_SLOT(learnsample, 
                                                   PL2_responsesSym))));

         nweights = S3get_nodeweights(tree);
         dnweights = REAL(nweights);
         dweights = REAL(VECTOR_ELT(bweights, b));
         for (i = 0; i < nobs; i++) {
             dnweights[i] = dweights[i];
         }
     
         C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum, 1);
         nodenum = 1;
         C_remove_weights(tree, 0);

         /* compute terminal node ids for all observations (not just those with weight > 0) */
         for (i = 0; i < nobs; i++)
            iwhere[i] = C_get_nodeID(tree, GET_SLOT(learnsample, PL2_inputsSym), 0.0, i, -1);
         
         if (get_trace(controls)) {
             /* progress bar; inspired by 
             http://avinashjoshi.co.in/2009/10/13/creating-a-progress-bar-in-c/ */
             Rprintf("[");
             /* Print the = until the current percentage */
             l = (int) ceil( ((double) b * 50.0) / B);
             for (j = 0; j < l; j++)
                 Rprintf("=");
             Rprintf(">");
             for (k = j; k < 50; k++)
                 Rprintf(" ");
             Rprintf("]");
             /* % completed */
                 Rprintf(" %3d%% completed", j * 2);
             /* To delete the previous line */
             Rprintf("\r");
             /* Flush all char in buffer */
             /* fflush(stdout); */
         }
     }
     if (get_trace(controls))
         Rprintf("\n");

     PutRNGstate();

     UNPROTECT(5);
     return(ans);
}
