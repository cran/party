
/**
    Node splitting and prediction
    *\file Predict.c
    *\author $Author: hothorn $
    *\date $Date: 2007-01-15 11:24:41 +0100 (Mon, 15 Jan 2007) $
*/
                
#include "party.h"


/**
    Split a node according to a splitting rule \n
    *\param node the current node with primary split specified
    *\param learnsample learning sample
    *\param control an object of class `TreeControl'
    *\todo outplace the splitting since there are at least 3 functions
           with nearly identical code
*/
                
void C_splitnode(SEXP node, SEXP learnsample, SEXP control) {

    SEXP weights, leftnode, rightnode, split;
    SEXP responses, inputs, whichNA;
    double cutpoint, *dx, *dweights, *leftweights, *rightweights;
    double sleft = 0.0, sright = 0.0;
    int *ix, *levelset, *iwhichNA;
    int nobs, i, nna;
                    
    weights = S3get_nodeweights(node);
    dweights = REAL(weights);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    nobs = get_nobs(learnsample);
            
    /* set up memory for the left daughter */
    SET_VECTOR_ELT(node, S3_LEFT, leftnode = allocVector(VECSXP, NODE_LENGTH));
    C_init_node(leftnode, nobs, 
        get_ninputs(learnsample), get_maxsurrogate(get_splitctrl(control)),
        ncol(get_jointtransf(GET_SLOT(learnsample, PL2_responsesSym))));
    leftweights = REAL(S3get_nodeweights(leftnode));

    /* set up memory for the right daughter */
    SET_VECTOR_ELT(node, S3_RIGHT, 
                   rightnode = allocVector(VECSXP, NODE_LENGTH));
    C_init_node(rightnode, nobs, 
        get_ninputs(learnsample), get_maxsurrogate(get_splitctrl(control)),
        ncol(get_jointtransf(GET_SLOT(learnsample, PL2_responsesSym))));
    rightweights = REAL(S3get_nodeweights(rightnode));

    /* split according to the primary split */
    split = S3get_primarysplit(node);
    if (has_missings(inputs, S3get_variableID(split))) {
        whichNA = get_missings(inputs, S3get_variableID(split));
        iwhichNA = INTEGER(whichNA);
        nna = LENGTH(whichNA);
    } else {
        nna = 0;
        whichNA = R_NilValue;
        iwhichNA = NULL;
    }
    
    if (S3is_ordered(split)) {
        cutpoint = REAL(S3get_splitpoint(split))[0];
        dx = REAL(get_variable(inputs, S3get_variableID(split)));
        for (i = 0; i < nobs; i++) {
            if (nna > 0) {
                if (i_in_set(i + 1, iwhichNA, nna)) continue;
            }
            if (dx[i] <= cutpoint) 
                leftweights[i] = dweights[i]; 
            else 
                leftweights[i] = 0.0;
            rightweights[i] = dweights[i] - leftweights[i];
            sleft += leftweights[i];
            sright += rightweights[i];
        }
    } else {
        levelset = INTEGER(S3get_splitpoint(split));
        ix = INTEGER(get_variable(inputs, S3get_variableID(split)));

        for (i = 0; i < nobs; i++) {
            if (nna > 0) {
                if (i_in_set(i + 1, iwhichNA, nna)) continue;
            }
            if (levelset[ix[i] - 1])
                leftweights[i] = dweights[i];
            else 
                leftweights[i] = 0.0;
            rightweights[i] = dweights[i] - leftweights[i];
            sleft += leftweights[i];
            sright += rightweights[i];
        }
    }
    
    /* for the moment: NA's go with majority */
    if (nna > 0) {
        for (i = 0; i < nna; i++) {
            if (sleft > sright) {
                leftweights[iwhichNA[i] - 1] = dweights[iwhichNA[i] - 1];
                rightweights[iwhichNA[i] - 1] = 0.0;
            } else {
                rightweights[iwhichNA[i] - 1] = dweights[iwhichNA[i] - 1];
                leftweights[iwhichNA[i] - 1] = 0.0;
            }
        }
    }
}


/**
    Get the terminal node for obs. number `numobs' of `newinputs' \n
    *\param subtree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param numobs observation number
    *\todo handle surrogate splits
*/

SEXP C_get_node(SEXP subtree, SEXP newinputs, 
                double mincriterion, int numobs) {

    SEXP split, whichNA, weights, ssplit, surrsplit;
    double cutpoint, x, *dweights, swleft, swright;
    int level, *levelset, i, ns;

    if (S3get_nodeterminal(subtree) || 
        REAL(S3get_maxcriterion(subtree))[0] < mincriterion) 
        return(subtree);
    
    split = S3get_primarysplit(subtree);

    /* missing values. Maybe store the proportions left / 
       right in each node? */
    if (has_missings(newinputs, S3get_variableID(split))) {
        whichNA = get_missings(newinputs, S3get_variableID(split));
    
        if (C_i_in_set(numobs, whichNA)) {
        
            surrsplit = S3get_surrogatesplits(subtree);
            ns = 0;
            i = numobs;      

            /* try to find a surrogate split */
            while(TRUE) {
    
                if (ns >= LENGTH(surrsplit)) break;
            
                ssplit = VECTOR_ELT(surrsplit, ns);
                if (has_missings(newinputs, S3get_variableID(ssplit))) {
                    if (INTEGER(get_missings(newinputs, 
                                             S3get_variableID(ssplit)))[i]) {
                        ns++;
                        continue;
                    }
                }

                cutpoint = REAL(S3get_splitpoint(ssplit))[0];
                x = REAL(get_variable(newinputs, S3get_variableID(ssplit)))[i];
                     
                if (S3get_toleft(ssplit)) {
                    if (x <= cutpoint) {
                        return(C_get_node(S3get_leftnode(subtree),
                                          newinputs, mincriterion, numobs));
                    } else {
                        return(C_get_node(S3get_rightnode(subtree),
                               newinputs, mincriterion, numobs));
                    }
                } else {
                    if (x <= cutpoint) {
                        return(C_get_node(S3get_rightnode(subtree),
                                          newinputs, mincriterion, numobs));
                    } else {
                        return(C_get_node(S3get_leftnode(subtree),
                               newinputs, mincriterion, numobs));
                    }
                }
                break;
            }

            /* if this was not successful, we go with the majority */
            weights = S3get_nodeweights(S3get_leftnode(subtree));
            dweights = REAL(weights);
            swleft = 0.0;
            for (i = 0; i < LENGTH(weights); i++)
                swleft += dweights[i];
            weights = S3get_nodeweights(S3get_rightnode(subtree));
            dweights = REAL(weights);
            swright = 0.0;
            for (i = 0; i < LENGTH(weights); i++)
                swright += dweights[i];
            if (swleft > swright) {
                return(C_get_node(S3get_leftnode(subtree), 
                                  newinputs, mincriterion, numobs));
            } else {
                return(C_get_node(S3get_rightnode(subtree), 
                                  newinputs, mincriterion, numobs));
            }
        }
    }
    
    if (S3is_ordered(split)) {
        cutpoint = REAL(S3get_splitpoint(split))[0];
        x = REAL(get_variable(newinputs, 
                     S3get_variableID(split)))[numobs];
        if (x <= cutpoint) {
            return(C_get_node(S3get_leftnode(subtree), 
                              newinputs, mincriterion, numobs));
        } else {
            return(C_get_node(S3get_rightnode(subtree), 
                              newinputs, mincriterion, numobs));
        }
    } else {
        levelset = INTEGER(S3get_splitpoint(split));
        level = INTEGER(get_variable(newinputs, 
                            S3get_variableID(split)))[numobs];
        /* level is in 1, ..., K */
        if (levelset[level - 1]) {
            return(C_get_node(S3get_leftnode(subtree), newinputs, 
                              mincriterion, numobs));
        } else {
            return(C_get_node(S3get_rightnode(subtree), newinputs, 
                              mincriterion, numobs));
        }
    }
}


/**
    R-Interface to C_get_node \n
    *\param subtree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param numobs observation number
*/

SEXP R_get_node(SEXP subtree, SEXP newinputs, SEXP mincriterion, 
                SEXP numobs) {
    return(C_get_node(subtree, newinputs, REAL(mincriterion)[0],
                      INTEGER(numobs)[0] - 1));
}


/**
    Get the node with nodeID `nodenum' \n
    *\param subtree a tree
    *\param nodenum a nodeID
*/

SEXP C_get_nodebynum(SEXP subtree, int nodenum) {
    
    if (nodenum == S3get_nodeID(subtree)) return(subtree);

    if (S3get_nodeterminal(subtree)) 
        error("no node with number %d\n", nodenum);

    if (nodenum < S3get_nodeID(S3get_rightnode(subtree))) {
        return(C_get_nodebynum(S3get_leftnode(subtree), nodenum));
    } else {
        return(C_get_nodebynum(S3get_rightnode(subtree), nodenum));
    }
}


/**
    R-Interface to C_get_nodenum \n
    *\param subtree a tree
    *\param nodenum a nodeID
*/

SEXP R_get_nodebynum(SEXP subtree, SEXP nodenum) {
    return(C_get_nodebynum(subtree, INTEGER(nodenum)[0]));
}


/**
    Get the prediction of a new observation\n
    *\param subtree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param numobs observation number
*/

SEXP C_get_prediction(SEXP subtree, SEXP newinputs, 
                      double mincriterion, int numobs) {
    return(S3get_prediction(C_get_node(subtree, newinputs, 
                            mincriterion, numobs)));
}


/**
    Get the weights for a new observation \n
    *\param subtree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param numobs observation number
*/

SEXP C_get_nodeweights(SEXP subtree, SEXP newinputs, 
                       double mincriterion, int numobs) {
    return(S3get_nodeweights(C_get_node(subtree, newinputs, 
                             mincriterion, numobs)));
}


/**
    Get the nodeID for a new observation \n
    *\param subtree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param numobs observation number
*/

int C_get_nodeID(SEXP subtree, SEXP newinputs,
                  double mincriterion, int numobs) {
     return(S3get_nodeID(C_get_node(subtree, newinputs, 
            mincriterion, numobs)));
}


/**
    R-Interface to C_get_nodeID \n
    *\param tree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
*/

SEXP R_get_nodeID(SEXP tree, SEXP newinputs, SEXP mincriterion) {

    SEXP ans;
    int nobs, i, *dans;
            
    nobs = get_nobs(newinputs);
    PROTECT(ans = allocVector(INTSXP, nobs));
    dans = INTEGER(ans);
    for (i = 0; i < nobs; i++)
         dans[i] = C_get_nodeID(tree, newinputs, REAL(mincriterion)[0], i);
    UNPROTECT(1);
    return(ans);
}


/**
    Get all predictions for `newinputs' \n
    *\param tree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param ans return value
*/

void C_predict(SEXP tree, SEXP newinputs, double mincriterion, SEXP ans) {
    
    int nobs, i;
    
    nobs = get_nobs(newinputs);    
    if (LENGTH(ans) != nobs) 
        error("ans is not of length %d\n", nobs);
        
    for (i = 0; i < nobs; i++)
        SET_VECTOR_ELT(ans, i, C_get_prediction(tree, newinputs, 
                       mincriterion, i));
}


/**
    R-Interface to C_predict \n
    *\param tree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
*/

SEXP R_predict(SEXP tree, SEXP newinputs, SEXP mincriterion) {

    SEXP ans;
    int nobs;
    
    nobs = get_nobs(newinputs);
    PROTECT(ans = allocVector(VECSXP, nobs));
    C_predict(tree, newinputs, REAL(mincriterion)[0], ans);
    UNPROTECT(1);
    return(ans);
}


/**
    Get the predictions from `where' nodes\n
    *\param tree a tree
    *\param where vector of nodeID's
    *\param ans return value
*/

void C_getpredictions(SEXP tree, SEXP where, SEXP ans) {

    int nobs, i, *iwhere;
    
    nobs = LENGTH(where);
    iwhere = INTEGER(where);
    if (LENGTH(ans) != nobs)
        error("ans is not of length %d\n", nobs);
        
    for (i = 0; i < nobs; i++)
        SET_VECTOR_ELT(ans, i, S3get_prediction(
            C_get_nodebynum(tree, iwhere[i])));
}


/**
    R-Interface to C_getpredictions\n
    *\param tree a tree
    *\param where vector of nodeID's
*/
            
SEXP R_getpredictions(SEXP tree, SEXP where) {

    SEXP ans;
    int nobs;
            
    nobs = LENGTH(where);
    PROTECT(ans = allocVector(VECSXP, nobs));
    C_getpredictions(tree, where, ans);
    UNPROTECT(1);
    return(ans);
}                        


/**
    Get the weights from `where' nodes\n
    *\param tree a tree
    *\param where vector of nodeID's
    *\param ans return value
*/

void C_getweights(SEXP tree, SEXP where, SEXP ans) {

    int nobs, i, *iwhere;
    
    nobs = LENGTH(where);
    iwhere = INTEGER(where);
    if (LENGTH(ans) != nobs)
        error("ans is not of length %d\n", nobs);
        
    for (i = 0; i < nobs; i++)
        SET_VECTOR_ELT(ans, i, S3get_nodeweights(
            C_get_nodebynum(tree, iwhere[i])));
}


/**
    R-Interface to C_getweigts \n
    *\param tree a tree
    *\param where vector of nodeID's
*/

SEXP R_getweights(SEXP tree, SEXP where) {

    SEXP ans;
    int nobs;
            
    nobs = LENGTH(where);
    PROTECT(ans = allocVector(VECSXP, nobs));
    C_getweights(tree, where, ans);
    UNPROTECT(1);
    return(ans);
}                        


/**
    Get the weights for all observations in  `newinputs'
    *\param tree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param ans return value
*/

void C_weights(SEXP tree, SEXP newinputs, 
               double mincriterion, SEXP ans) {
    
    int nobs, i;
    
    nobs = get_nobs(newinputs);    
    if (LENGTH(ans) != nobs) 
        error("ans is not of length %d\n", nobs);
        
    for (i = 0; i < nobs; i++)
        SET_VECTOR_ELT(ans, i, C_get_nodeweights(tree, newinputs, 
                       mincriterion, i));
}


/**
    R-Interface to C_weights \n
    *\param tree a tree
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
*/

SEXP R_weights(SEXP tree, SEXP newinputs, SEXP mincriterion) {

    SEXP ans;
    int nobs;
    
    nobs = get_nobs(newinputs);
    PROTECT(ans = allocVector(VECSXP, nobs));
    C_weights(tree, newinputs, REAL(mincriterion)[0], ans);
    UNPROTECT(1);
    return(ans);
}


/**
    Predictions from RandomForest objects
    *\param forest a list of trees
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param oobpred a logical indicating out-of-bag predictions
*/

SEXP R_predictRF(SEXP forest, SEXP newinputs, SEXP mincriterion, SEXP oobpred) {

    SEXP ans, tmp, tree;
    int ntrees, nobs, i, b, j, q, iwhere, oob = 0, count = 0;
    
    if (LOGICAL(oobpred)[0]) oob = 1;
    
    nobs = get_nobs(newinputs);
    ntrees = LENGTH(forest);
    q = LENGTH(S3get_prediction(
                   C_get_nodebynum(VECTOR_ELT(forest, 0), 1)));

    if (oob) {
        if (LENGTH(S3get_nodeweights(
                       C_get_nodebynum(VECTOR_ELT(forest, 0), 1))) != nobs)
            error("number of observations don't match");
    }    
    
    PROTECT(ans = allocVector(VECSXP, nobs));
    
    for (i = 0; i < nobs; i++) {
        count = 0;
        SET_VECTOR_ELT(ans, i, allocVector(REALSXP, q));
        for (j = 0; j < q; j++)
                    REAL(VECTOR_ELT(ans, i))[j] = 0.0;
        for (b = 0; b < ntrees; b++) {
            tree = VECTOR_ELT(forest, b);

            if (oob && 
                REAL(S3get_nodeweights(C_get_nodebynum(tree, 1)))[i] > 0.0) 
                continue;

            iwhere = C_get_nodeID(tree, newinputs, REAL(mincriterion)[0], i);
            tmp = S3get_prediction(C_get_nodebynum(tree, iwhere));
            for (j = 0; j < q; j++)
                REAL(VECTOR_ELT(ans, i))[j] += REAL(tmp)[j];
            count++;
        }
        if (count == 0) 
            error("cannot compute out-of-bag predictions for obs ", i + 1);
        for (j = 0; j < q; j++)
            REAL(VECTOR_ELT(ans, i))[j] = REAL(VECTOR_ELT(ans, i))[j] / count;
    }
    UNPROTECT(1);
    return(ans);
}

/**
    Predictions from RandomForest objects, based in total weights
    *\param forest a list of trees
    *\param response a matrix of (transformed) response values
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param oobpred a logical indicating out-of-bag predictions
*/

SEXP R_predictRF2(SEXP forest, SEXP response, SEXP newinputs, 
                  SEXP mincriterion, SEXP oobpred) {

    SEXP ans, tmp, tree, w;
    int ntrees, nobs, i, b, j, q, n, iwhere, oob = 0;
    double *dtmp, *dw, sumw = 0.0;

    if (LOGICAL(oobpred)[0]) oob = 1;
    
    nobs = get_nobs(newinputs);
    ntrees = LENGTH(forest);
    n = nrow(response);
    q = ncol(response);

    if (oob) {
        if (n != nobs)
            error("number of observations don't match");
    }    
    
    PROTECT(ans = allocVector(VECSXP, nobs));
    PROTECT(w = allocMatrix(REALSXP, 1, n));
    dw = REAL(w);
    
    for (i = 0; i < nobs; i++) {

        SET_VECTOR_ELT(ans, i, allocVector(REALSXP, q));
        for (j = 0; j < n; j++)
            dw[j] = 0.0;

        for (b = 0; b < ntrees; b++) {
            tree = VECTOR_ELT(forest, b);

            if (oob && 
                REAL(S3get_nodeweights(C_get_nodebynum(tree, 1)))[i] > 0.0) 
                continue;

            iwhere = C_get_nodeID(tree, newinputs, REAL(mincriterion)[0], i);
            tmp = S3get_nodeweights(C_get_nodebynum(tree, iwhere));
            dtmp = REAL(tmp);
            
            for (j = 0; j < n; j++)
                dw[j] += dtmp[j];
        }
        
        C_matprod(dw, 1, n, REAL(response), n, q, REAL(VECTOR_ELT(ans, i)));

        sumw = 0.0;
        for (j = 0; j < n; j++)
            sumw += dw[j];

        for (j = 0; j < q; j++)
            REAL(VECTOR_ELT(ans, i))[j] = REAL(VECTOR_ELT(ans, i))[j] / sumw;
    }
    UNPROTECT(2);
    return(ans);
}

/**
    Predictions weights from RandomForest objects
    *\param forest a list of trees
    *\param newinputs an object of class `VariableFrame'
    *\param mincriterion overwrites mincriterion used for tree growing
    *\param oobpred a logical indicating out-of-bag predictions
*/

SEXP R_predictRF_weights(SEXP forest, SEXP newinputs, SEXP mincriterion, SEXP oobpred) {

    SEXP ans, tree, bw;
    int ntrees, nobs, i, b, j, q, iwhere, oob = 0, count = 0, ntrain;
    double *dtmp;
    
    if (LOGICAL(oobpred)[0]) oob = 1;
    
    nobs = get_nobs(newinputs);
    ntrees = LENGTH(forest);
    q = LENGTH(S3get_prediction(
                   C_get_nodebynum(VECTOR_ELT(forest, 0), 1)));

    if (oob) {
        if (LENGTH(S3get_nodeweights(
                       C_get_nodebynum(VECTOR_ELT(forest, 0), 1))) != nobs)
            error("number of observations don't match");
    }    
    
    tree = VECTOR_ELT(forest, 0);
    ntrain = LENGTH(S3get_nodeweights(C_get_nodebynum(tree, 1)));
    
    PROTECT(ans = allocVector(VECSXP, nobs));
    
    for (i = 0; i < nobs; i++) {
        count = 0;
        SET_VECTOR_ELT(ans, i, bw = allocVector(REALSXP, ntrain));
        for (j = 0; j < ntrain; j++)
            REAL(bw)[j] = 0.0;
        for (b = 0; b < ntrees; b++) {
            tree = VECTOR_ELT(forest, b);

            if (oob && 
                REAL(S3get_nodeweights(C_get_nodebynum(tree, 1)))[i] > 0.0) 
                continue;

            iwhere = C_get_nodeID(tree, newinputs, REAL(mincriterion)[0], i);
            dtmp = REAL(S3get_nodeweights(C_get_nodebynum(tree, iwhere)));
            for (j = 0; j < ntrain; j++)
                REAL(bw)[j] += dtmp[j];
            count++;
        }
        if (count == 0) 
            error("cannot compute out-of-bag predictions for obs ", i + 1);
    }
    UNPROTECT(1);
    return(ans);
}
