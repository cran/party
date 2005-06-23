
# $Id: Variables.R,v 1.5 2005/06/22 09:32:30 hothorn Exp $

initVariableFrame.df <- function(obj, trafo = NULL, scores = NULL, ...) {

    RET <- new("VariableFrame", nrow(obj), ncol(obj))
    
    is_ordinal <- sapply(obj, is.ordered)
    is_nominal <- sapply(obj, is.factor) & !is_ordinal

    jt <- c()

    ### assign user-specified scores to variables in `obj'
    if (!is.null(scores)) {
        if (!is.list(scores) || is.null(names(scores)))
            stop(sQuote("scores"), " is not a named list")
        scores <- scores[names(scores) %in% colnames(obj)]
    }
    if (!is.null(scores)) {
        tmp <- sapply(names(scores), function(n) {
            if (!(is.factor(obj[[n]]) && is.ordered(obj[[n]])) || 
                nlevels(obj[[n]]) != length(scores[[n]]))
                stop("cannot assign scores to variable ", sQuote(n))
            attr(obj[[n]], "scores") <- scores[[n]]
        })
    }

    ### speedup: extract the slots and re-assign them afterwards
    scores <- RET@scores
    levels <- RET@levels
    variables <- obj
    transformations <- RET@transformations
    ordering <- RET@ordering
    has_missings <- RET@has_missings
    whichNA <- RET@whichNA
 
    for (j in 1:ncol(obj)) {
        x <- variables[[j]]
        if (is_ordinal[j]) {
            sc <- attr(x, "scores")
            if (is.null(sc))
                sc <- 1:nlevels(x)
            storage.mode(sc) <- "double"
            scores[[j]] <- sc
        }

        if (is.factor(x)) {
            xt <- sapply(levels(x), function(l) as.numeric(x == l))
            storage.mode(xt) <- "double"
            levels[[j]] <- attr(x, "levels")

            ### <FIXME> storage mode of nominal and ordered 
            ### factors are different!!!       
            if (!is_ordinal[j]) {
                storage.mode(RET@variables[[j]]) <- "integer"
                tordering <- NULL
            } else {
                storage.mode(variables[[j]]) <- "double"
                tordering <- order(as.numeric(x))
            }
            ### </FIXME>
        } else {
            if (!is.null(trafo)) 
                xt <- matrix(trafo(x), ncol = 1)
            else 
                xt <- matrix(x, ncol = 1)
            storage.mode(xt) <- "double"
            tordering <- order(x)
            storage.mode(tordering) <- "integer"
            storage.mode(variables[[j]]) <- "double"
        }
        nas <- is.na(x)
        xt[nas,] <- 0
        transformations[[j]] <- xt
        ordering[[j]] <- tordering
        has_missings[j] <- any(nas)
        whichNA[[j]] <- which(nas)
        ### this is suboptimal
        jt <- cbind(jt, xt)
    }            
    RET@jointtransf <- jt
    RET@is_nominal <- is_nominal
    RET@is_ordinal <- is_ordinal
    RET@is_censored <- rep(FALSE, ncol(obj))
    RET@variables <- variables
    RET@scores <- scores
    RET@levels <- levels
    RET@transformations <- transformations
    RET@ordering <- ordering
    RET@has_missings <- has_missings
    RET@whichNA <- whichNA
    RET
}

setGeneric(name = "initVariableFrame",
           def = function(obj, trafo = NULL, ...)
               standardGeneric("initVariableFrame")
)

setClassUnion("function_OR_NULL", c("function", "NULL"))

setMethod("initVariableFrame", 
    signature = c("data.frame", "function_OR_NULL"), 
    definition = initVariableFrame.df
)

initVariableFrame.Surv <- function(obj, trafo = logrank_trafo, ...) {

    if (is.null(trafo)) trafo <- logrank_trafo
    RET <- new("VariableFrame", nrow(obj), as.integer(1))
    RET@variables <- list(obj)

    RET@is_nominal <- FALSE
    RET@is_ordinal <- FALSE
    
    RET@transformations <- list(matrix(trafo(obj), ncol = 1))

    RET@ordering <- list(order(RET@transformations[[1]]))
    RET@has_missings <- any(is.na(obj))
    RET@whichNA <- list(which(is.na(obj)))
    RET@jointtransf <- RET@transformations[[1]]
    RET@is_censored <- TRUE
    RET
}

setMethod("initVariableFrame", signature = c("Surv", "function_OR_NULL"),
    definition = initVariableFrame.Surv
)

initVariableFrame.matrix <- function(obj, trafo = NULL, ...) {

    RET <- new("VariableFrame", nrow(obj), ncol(obj))
    objDF <- as.data.frame(obj)
    class(objDF) <- "list"
    RET@variables <- objDF

    RET@is_nominal <- rep(FALSE, ncol(obj))
    RET@is_ordinal <- rep(FALSE, ncol(obj))
   
    if (!is.null(trafo)) 
        RET@transformations <- lapply(objDF, function(x) 
                                      matrix(trafo(x), ncol = 1))
    else
        RET@transformations <- lapply(objDF, function(x) matrix(x, ncol = 1))

    RET@ordering <- lapply(RET@transformations, order)
    RET@has_missings <- unlist(lapply(RET@transformations, function(x) 
                               any(is.na(x))))
    if (any(is.na(obj)))
        RET@whichNA <- lapply(objDF, function(x) which(is.na(x)))
    RET@jointtransf <- matrix(unlist(RET@transformations), 
                              nrow = nrow(obj), byrow = TRUE)
    RET@is_censored <- rep(FALSE, ncol(obj))
    RET
}

setMethod("initVariableFrame", signature = c("matrix", "function_OR_NULL"),
    definition = initVariableFrame.matrix
)


setGeneric(name = "response",
           def = function(object, ...)
               standardGeneric("response")
)

setMethod("response",
    signature = "BinaryTree",
    definition = function(object) object@responses@variables
)
