
# $Id: RandomForest.R 2538 2006-04-04 15:47:15Z hothorn $

### the fitting procedure
cforestfit <- function(object, controls, weights = NULL, fitmem = NULL, 
                       ntree = 500, ...) {

    if (!extends(class(object), "LearningSample"))
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "TreeControl"))
        stop(sQuote("controls"), " is not of class ", sQuote("TreeControl"))

    if (is.null(fitmem)) 
        fitmem <- ctree_memory(object, TRUE)
    if (!extends(class(fitmem), "TreeFitMemory"))
        stop(sQuote("fitmem"), " is not of class ", sQuote("TreeFitMemory"))

    if (is.null(weights))
        weights <- object@weights
    if (length(weights) != object@nobs || storage.mode(weights) != "double")
        stop(sQuote("weights"), " are not a double vector of ", 
             object@nobs, " elements")

    where <- rep(0, object@nobs)
    storage.mode(where) <- "integer"

    ### grow the tree
    ensemble <- .Call("R_Ensemble", object, weights, fitmem, controls,
                      as.integer(ntree), PACKAGE = "party")

    ### prepare the returned object
    RET <- new("RandomForest")
    RET@ensemble <- ensemble
    RET@responses <- object@responses
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
        RET <- lapply(pw, function(w) w %*% response@jointtransf / sum(w))
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
        RET <- unlist(cdresp)
        RET <- matrix(unlist(RET),
                      nrow = length(RET), byrow = TRUE)
        ### <FIXME> what about multivariate responses?
        if (response@ninputs == 1)
            colnames(RET) <- names(response@variables)
        ### </FIXME>
        return(RET)
    }

    RET@prediction_weights <- function(newdata = NULL, 
                                       mincriterion = 0, OOB = FALSE) {

        if (is.null(newdata)) {
            newinp <- object@inputs
        } else {
            newinp <- object@menv@get("input", data = newdata)
            newinp <- initVariableFrame(newinp, trafo = NULL)
        }

        return(.Call("R_predictRF_weights", ensemble, newinp, mincriterion,
                     OOB && is.null(newdata), PACKAGE = "party"))
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

### the top-level convenience function
cforest <- function(formula, data = list(), subset = NULL, weights = NULL, 
                    controls = ctree_control(teststat = "max", 
                                             testtype = "Teststatistic", 
                                             mincriterion = qnorm(0.9), mtry = 5, 
                                             savesplitstats = FALSE),
                    xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL, ntree = 500) {

    ### setup learning sample
    ls <- dpp(RandomForest, formula, data, subset, xtrafo = xtrafo, 
              ytrafo = ytrafo, scores = scores)

    ### setup memory
    fitmem <- ctree_memory(ls, TRUE)

    ### fit and return a conditional tree
    fit(RandomForest, ls, controls = controls, weights = weights, 
        fitmem = fitmem, ntree = ntree)
}
