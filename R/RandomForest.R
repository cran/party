
# $Id: RandomForest.R 3661 2007-07-23 09:44:30Z hothorn $

### the fitting procedure
cforestfit <- function(object, controls, weights = NULL, fitmem = NULL, ...) {

    if (!extends(class(object), "LearningSample"))
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "ForestControl"))
        stop(sQuote("controls"), " is not of class ", sQuote("ForestControl"))

    if (is.null(fitmem)) 
        fitmem <- ctree_memory(object, TRUE)
    if (!extends(class(fitmem), "TreeFitMemory"))
        stop(sQuote("fitmem"), " is not of class ", sQuote("TreeFitMemory"))

    if (is.null(weights))
        weights <- object@weights
    storage.mode(weights) <- "double"
    if (length(weights) != object@nobs || storage.mode(weights) != "double")
        stop(sQuote("weights"), " are not a double vector of ", 
             object@nobs, " elements")

    ### grow the tree
    bweights <- vector(mode = "list", length = controls@ntree)
    bwhere <- vector(mode = "list", length = controls@ntree)
    ensemble <- .Call("R_Ensemble", object, weights, bwhere, bweights, fitmem, controls,
                      PACKAGE = "party")

    ### prepare the returned object
    RET <- new("RandomForest")
    RET@ensemble <- ensemble
    RET@where <- bwhere
    RET@weights <- bweights
    RET@responses <- object@responses
    if (inherits(object, "LearningSampleFormula"))
        RET@data <- object@menv

    ### (estimated) conditional distribution of the response given the
    ### covariates
    RET@cond_distr_response <- function(newdata = NULL, mincriterion = 0, ...) { 
        
        pw <- RET@prediction_weights(newdata = newdata, mincriterion =
                                     mincriterion, ...)

        response <- object@responses

        ### survival: estimated Kaplan-Meier
        if (any(response@is_censored)) {
            resp <- response@variables[[1]]
            RET <- lapply(pw, function(w) 
                survival:::survfit(resp, weights = w, subset = w > 0))
            return(RET)
        }

        ### classification: estimated class probabilities
        ### regression: the means, not really a distribution
        RET <- lapply(pw, function(w) w %*% response@predict_trafo / sum(w))
        return(RET)
    }

    ### predict in the response space, always!
    RET@predict_response <- function(newdata = NULL, mincriterion = 0, ...) { 

        cdresp <- RET@cond_distr_response(newdata = newdata, 
                                          mincriterion = mincriterion, ...)

        response <- object@responses
        ### classification: classes
        if (all(response@is_nominal || response@is_ordinal)) {
            lev <- levels(response@variables[[1]])
            RET <- factor(lev[unlist(lapply(cdresp, which.max))],
                          levels = levels(response@variables[[1]]))
            return(RET)
        }

        ### survival: median survival time
        if (any(response@is_censored)) {
            RET <- sapply(cdresp, mst)
            return(RET)
        }

        ### regression: mean (median would be possible)
        RET <- matrix(unlist(cdresp),
                      nrow = length(cdresp), byrow = TRUE)
        ### <FIXME> what about multivariate responses?
        colnames(RET) <- names(response@variables)
        ### </FIXME>
        return(RET)
    }

    RET@prediction_weights <- function(newdata = NULL, 
                                       mincriterion = 0, OOB = FALSE) {

        newinp <- newinputs(object, newdata)

        return(.Call("R_predictRF_weights", ensemble, bwhere, bweights, 
                     newinp, mincriterion, OOB && is.null(newdata), 
                     PACKAGE = "party"))
    }
    return(RET)
}

### the unfitted forest, an object of class `StatModel'
### see package `modeltools'
RandomForest <- new("StatModel",
                    capabilities = new("StatModelCapabilities"),
                    name = "random forest",
                    dpp = ctreedpp,
                    fit = cforestfit,
                    predict = function(object, ...) 
                        object@predict_response(...))

cforest_control <- function(teststat = "max", 
                            testtype = "Teststatistic",
                            mincriterion = qnorm(0.9),
                            savesplitstats = FALSE,
                            ntree = 500, mtry = 5, replace = TRUE, 
                            fraction = 0.632, ...) {
    RET <- ctree_control(teststat = teststat, testtype = testtype,
                         mincriterion = mincriterion, 
                         savesplitstats = savesplitstats, 
                         mtry = mtry, ...)
    class(RET) <- "ForestControl"
    RET@ntree <- as.integer(ntree)
    RET@replace <- replace
    RET@fraction <- as.double(fraction)
    if (!validObject(RET))
        stop("RET is not a valid object of class", class(RET))
    RET
}
    
### the top-level convenience function
cforest <- function(formula, data = list(), subset = NULL, weights = NULL, 
                    controls = cforest_control(),
                    xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL) {

    ### setup learning sample
    ls <- dpp(RandomForest, formula, data, subset, xtrafo = xtrafo, 
              ytrafo = ytrafo, scores = scores)

    ### setup memory
    fitmem <- ctree_memory(ls, TRUE)

    ### fit and return a conditional tree
    fit(RandomForest, ls, controls = controls, weights = weights, 
        fitmem = fitmem)
}

###
### variable importance for `cforest'
###
### see ?importance (in `randomForest'), too
###
###

### extract ID of _all_ variables the tree uses for splitting
varIDs <- function(node) {

    v <- c()
    foo <- function(node) {
        if (node[[4]]) return(NULL)
        v <<- c(v, node[[5]][[1]])
        foo(node[[8]])
        foo(node[[9]])
    }
    foo(node)
    return(v)
}

varimp <- function(x, mincriterion = 0.0) {

    inputs <- x@data@get("input")
    response <- x@responses
    if (length(response@variables) != 1)
        stop("cannot compute variable importance measure for multivariate response")
    y <- x@responses@variables[[1]]
    inp <- initVariableFrame(inputs, trafo = NULL) 
    if (!all(complete.cases(inp@variables)))
        stop("cannot compute variable importance measure with missing values")
    tmp <- inp
    ### jt <- response@jointtransf

    CLASS <- all(response@is_nominal || response@is_ordinal)
    if (CLASS) {
        error <- function(x, oob) 
            mean((levels(y)[sapply(x, which.max)] != y)[oob])
    } else {
        error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
    }

    perror <- matrix(0, nrow = length(x@ensemble), ncol = ncol(inputs))
    colnames(perror) <- colnames(inputs)

    for (b in 1:length(x@ensemble)) {

        tree <- x@ensemble[[b]]
        oob <- x@weights[[b]] == 0

        p <- .Call("R_predict", tree, inp, mincriterion,
                   PACKAGE = "party")
        eoob <- error(p, oob)

        for (j in unique(varIDs(tree))) {
            perm <- sample(which(oob))
            tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
            ### <FIXME> check for NA's
            ### tmp@whichNA[[j]] <- which(is.na(tmp@variables[[j]]))
            ### </FIXME>

            p <- .Call("R_predict", tree, tmp, mincriterion,
                       PACKAGE = "party")

            perror[b, j] <- (error(p, oob) - eoob)
            tmp <- inp
        }
    }
    data.frame("MeanDecreaseAccuracy" = colMeans(perror), 
               "Standard Deviation" = apply(perror, 2, sd))
}
