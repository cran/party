
# $Id: Variables.R,v 1.10 2006/01/17 15:28:01 hothorn Exp $

ptrafo <- function(data, numeric_trafo = id_trafo, 
   factor_trafo = function(x) model.matrix(~ x - 1), 
   surv_trafo = logrank_trafo, var_trafo = NULL)

    trafo(data = data, numeric_trafo = numeric_trafo, factor_trafo =
          factor_trafo, surv_trafo = surv_trafo, var_trafo = var_trafo)


initVariableFrame.df <- function(obj, trafo = ptrafo, scores = NULL, ...) {

    RET <- new("VariableFrame", nrow(obj), ncol(obj))
    
    is_ordinal <- sapply(obj, is.ordered)
    is_nominal <- sapply(obj, is.factor) & !is_ordinal

    ### assign user-specified scores to variables in `obj'
    if (!is.null(scores)) {
        if (!is.list(scores) || is.null(names(scores)))
            stop(sQuote("scores"), " is not a named list")
        scores <- scores[names(scores) %in% colnames(obj)]
    }
    if (!is.null(scores)) {
        for (n in names(scores)) {
            if (!(is.factor(obj[[n]]) && is.ordered(obj[[n]])) || 
                nlevels(obj[[n]]) != length(scores[[n]]))
                stop("cannot assign scores to variable ", sQuote(n))
            attr(obj[[n]], "scores") <- scores[[n]]
        }
    }

    RET@scores <- lapply(obj, function(x) {
        sc <- NULL
        if (is.ordered(x)) {
            sc <- attr(x, "scores")
            if (is.null(sc)) sc <- 1:nlevels(x)
            storage.mode(sc) <- "double"
        }
        sc
    })

    ### transformations
    jt <- trafo(obj)

    ### for each variable
    xt <- vector(mode = "list", length = ncol(obj))
    for (i in 1:ncol(obj))
        xt[[i]] <- jt[,attr(jt, "assign") == i, drop = FALSE]

    ### ordering
    ordering <- lapply(obj, function(x) {
        if (is.factor(x) && !is.ordered(x)) return(NULL)
        if (is.ordered(x)) return(as.integer(order(as.numeric(x))))
        as.integer(order(x))
    })

    ### div.
    levels <- lapply(obj, function(x) if(is.factor(x)) levels(x))
    whichNA <- lapply(obj, function(x) which(is.na(x)))
    has_missings <- sapply(obj, function(x) any(is.na(x)))
    censored <- sapply(obj, function(x) inherits(x, "Surv"))

    ### some "handwork" 
    for (j in 1:ncol(obj)) {
        x <- obj[[j]]

        if (censored[j]) 
            ordering[[j]] <- as.integer(order(xt[[j]]))

        if (is.factor(x)) {

            ### <FIXME> storage mode of nominal and ordered 
            ### factors are different!!!       
            if (!is_ordinal[j]) {
                storage.mode(obj[[j]]) <- "integer"
            } else {
                storage.mode(obj[[j]]) <- "double"
            }
            ### </FIXME>
        } else {
            storage.mode(obj[[j]]) <- "double"
        }
        nas <- is.na(x)
        xt[[j]][nas] <- 0
        jt[nas, attr(jt, "assign") == j] <- 0
    }            

    RET@jointtransf <- jt
    RET@transformations <- xt
    RET@is_nominal <- is_nominal
    RET@is_ordinal <- is_ordinal
    RET@is_censored <- censored
    RET@variables <- obj
    RET@levels <- levels
    RET@ordering <- ordering
    RET@has_missings <- has_missings
    RET@whichNA <- whichNA
    RET
}

setGeneric(name = "initVariableFrame",
           def = function(obj, ...)
               standardGeneric("initVariableFrame")
)

setMethod("initVariableFrame", 
    signature = "data.frame",
    definition = initVariableFrame.df
)

setGeneric(name = "response",
           def = function(object, ...)
               standardGeneric("response")
)

setMethod("response",
    signature = "BinaryTree",
    definition = function(object) object@responses@variables
)

get_variables <- function(x)
    x@variables

    