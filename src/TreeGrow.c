
/**
    The tree growing recursion
    *\file TreeGrow.c
    *\author $Author: thothorn $
    *\date $Date: 2026-03-25 13:58:15 +0100 (Wed, 25 Mar 2026) $
*/

#include "party.h"


/**
    The main tree growing function, handles the recursion. \n
    *\param node  a list representing the current node
    *\param learnsample an object of class `LearningSample'
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param where a pointer to an integer vector of n-elements
    *\param nodenum a pointer to a integer vector of length 1
    *\param depth an integer giving the depth of the current node
*/

void C_TreeGrow(SEXP node, SEXP learnsample, SEXP fitmem, 
                SEXP controls, int *where, int *nodenum, int depth) {

    SEXP weights, tc, sc;
    int nobs, i, stop;
    double *dweights;
    
    weights = S3get_nodeweights(node);
    
    /* stop if either stumps have been requested or 
       the maximum depth is exceeded */
    PROTECT(tc = get_tgctrl(controls));
    PROTECT(sc = get_splitctrl(controls));
    stop = (nodenum[0] == 2 || nodenum[0] == 3) && 
           get_stump(tc);
    stop = stop || !check_depth(tc, depth);
    
    if (stop)
        C_Node(node, learnsample, weights, fitmem, controls, 1, depth);
    else
        C_Node(node, learnsample, weights, fitmem, controls, 0, depth);
    
    S3set_nodeID(node, nodenum[0]);    
    
    if (!S3get_nodeterminal(node)) {

        C_splitnode(node, learnsample, controls);

        /* determine surrogate splits and split missing values */
        if (get_maxsurrogate(sc) > 0) {
            C_surrogates(node, learnsample, weights, controls, fitmem);
            C_splitsurrogate(node, learnsample);
        }
            
        nodenum[0] += 1;
        C_TreeGrow(S3get_leftnode(node), learnsample, fitmem, 
                   controls, where, nodenum, depth + 1);

        nodenum[0] += 1;                                      
        C_TreeGrow(S3get_rightnode(node), learnsample, fitmem, 
                   controls, where, nodenum, depth + 1);
                   
    } else {
        dweights = REAL(weights);
        nobs = get_nobs(learnsample);
        for (i = 0; i < nobs; i++)
            if (dweights[i] > 0) where[i] = nodenum[0];
    } 
    UNPROTECT(2);
}


/**
    R-interface to C_TreeGrow\n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
    *\param where a vector of node indices for each observation
*/

SEXP R_TreeGrow(SEXP learnsample, SEXP weights, SEXP controls) {
            
     SEXP ans, tree, where, nweights, fitmem, sc, responses;
     double *dnweights, *dweights;
     int nobs, i, nodenum = 1, *iwhere;


     GetRNGstate();
     
     PROTECT(fitmem = ctree_memory(learnsample, PROTECT(ScalarLogical(1))));
     PROTECT(sc = get_splitctrl(controls));
     PROTECT(responses = GET_SLOT(learnsample, PL2_responsesSym));
     nobs = get_nobs(learnsample);
     PROTECT(ans = allocVector(VECSXP, 2));
     SET_VECTOR_ELT(ans, 0, where = allocVector(INTSXP, nobs));
     iwhere = INTEGER(where);
     for (int i = 0; i < nobs; i++) iwhere[i] = 0;
     SET_VECTOR_ELT(ans, 1, tree = allocVector(VECSXP, NODE_LENGTH));
     C_init_node(tree, nobs, get_ninputs(learnsample), get_maxsurrogate(sc),
                 ncol(PROTECT(get_predict_trafo(responses))));

     nweights = S3get_nodeweights(tree);
     dnweights = REAL(nweights);
     dweights = REAL(weights);
     for (i = 0; i < nobs; i++) dnweights[i] = dweights[i];
     
     C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum, 1);

     if (get_remove_weights(controls))
         C_remove_weights(tree, 0);
     
     PutRNGstate();
     
     UNPROTECT(6);
     return(ans);
}
