
# $Id: Variables.R 3259 2007-02-02 10:22:45Z hothorn $

### factor handling
ff_trafo <- function(x) {
    ### temporarily define `na.pass' as na.action
    opt <- options()
    on.exit(options(opt))
    options(na.action = na.pass)
    ### construct design matrix _without_ intercept
    mm <- model.matrix(~ x - 1)
    colnames(mm) <- levels(x)  
    ### remove unused levels   
    mm <- mm[,colSums(mm, na.rm = TRUE) > 0,drop = FALSE]
    return(mm)
}

ptrafo <- function(data, numeric_trafo = id_trafo, 
   factor_trafo = ff_trafo, surv_trafo = logrank_trafo, var_trafo = NULL)

    trafo(data = data, numeric_trafo = numeric_trafo, factor_trafo =
          factor_trafo, surv_trafo = surv_trafo, var_trafo = var_trafo)


initVariableFrame.df <- function(obj, trafo = ptrafo, scores = NULL, response = FALSE, ...) {

    if (response) {
        RET <- new("ResponseFrame", nrow(obj), ncol(obj))
        tmp <- lapply(obj, function(x) {
            if (is.factor(x)) return(ff_trafo(x))
            ### FIXME
            if (inherits(x, "Surv")) return(logrank_trafo(x))
            return(x)
        })
        RET@predict_trafo <- as.matrix(as.data.frame(tmp))
        storage.mode(RET@predict_trafo) <- "double"
    } else {
        RET <- new("VariableFrame", nrow(obj), ncol(obj))
    }
    
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
            if (any(order(scores[[n]]) != 1:length(scores[[n]])))
                stop("scores are not increasingly ordered")
            attr(obj[[n]], "scores") <- scores[[n]]
        }
    }

    RET@scores <- lapply(obj, function(x) {
        sc <- 0
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
    rm(jt)

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
            ### replace ordinal factors by their numeric scores
            if (is_ordinal[j]) {
                xt[[j]] <- matrix(RET@scores[[j]][obj[[j]]], ncol = 1)
                storage.mode(xt[[j]]) <- "double"
                ### R 2.5.0 does not allow to change the storage mode of factors
                class(obj[[j]]) <- "was_ordered"
                storage.mode(obj[[j]]) <- "double"
            } else {
                storage.mode(obj[[j]]) <- "integer"
            }
        } else {
            storage.mode(obj[[j]]) <- "double"
        }
        nas <- is.na(x)
        xt[[j]][nas] <- 0
    }            

    RET@transformations <- xt
    RET@is_nominal <- is_nominal
    RET@is_ordinal <- is_ordinal
    RET@is_censored <- censored
    RET@variables <- obj
    RET@levels <- levels
    RET@ordering <- ordering
    RET@has_missings <- has_missings
    RET@whichNA <- whichNA

    if (response) {
        RET@test_trafo <- as.matrix(as.data.frame(xt))
        storage.mode(RET@test_trafo) <- "double"
    }
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

    