
# $Id: zInitMethods.R 82 2005-06-09 13:40:51Z hothorn $

setMethod(f = "initialize", signature = "ExpectCovar",
    definition = function(.Object, pq = 1) {
        pq <- as.integer(pq)
        .Object@expectation <- rep(0, pq)
        .Object@covariance <- matrix(0, nrow = pq, ncol = pq)
        .Object@dimension  <- pq
        .Object
    }
)


setMethod(f = "initialize", signature = "ExpectCovarInfluence",
    definition = function(.Object, q) {
        .Object@expectation <- rep(0, q)
        .Object@covariance <- matrix(0, nrow = q, ncol = q)
        .Object@dimension  <- q
        .Object@sumweights <- 0.0
        .Object
    }
)


setMethod(f = "initialize", signature = "LinStatExpectCovar",
    definition = function(.Object, p, q) {
        .Object@expectation <- rep(0, p*q)
        .Object@covariance <- matrix(0, nrow = p*q, ncol = p*q)
        .Object@dimension  <- p*q
        .Object@linearstatistic <- rep(0, p*q)
        .Object@expcovinf <- new("ExpectCovarInfluence", q)
        .Object
    }
)

setMethod(f = "initialize", signature = "LinStatExpectCovarMPinv",
    definition = function(.Object, p, q) {
        .Object@expectation <- rep(0, p*q)
        .Object@covariance <- matrix(0, nrow = p*q, ncol = p*q)
        .Object@MPinv <- matrix(0, nrow = p*q, ncol = p*q)
        .Object@rank <- as.double(0.0)
        .Object@svdmem <- new("svd_mem", p*q)
        .Object@dimension  <- p*q
        .Object@linearstatistic <- rep(0, p*q)
        .Object@expcovinf <- new("ExpectCovarInfluence", q)
        .Object
    }
)



setMethod(f = "initialize", signature = "svd_mem",
    definition = function(.Object, p) {
        if (p <= 0) stop(sQuote("p"), " is not a positive integer")
        .Object@p <- p
        .Object@method <- "dgesdd"
        .Object@jobu <- "S"
        .Object@jobv <- ""
        .Object@u <- matrix(0, nrow = p, ncol = p)
        .Object@v <- matrix(0, nrow = p, ncol = p)
        .Object@s <- rep(0, p)
        .Object
    }
)

setMethod(f = "initialize", signature = "VariableFrame",
    definition = function(.Object, nobs, ninputs) {
        if (nobs <= 0 || ninputs <= 0) stop(sQuote("nobs"), " or ", 
            sQuote("ninputs"), " is not a positive integer")
        .Object@nobs <- nobs
        .Object@ninputs <- ninputs
        .Object@variables <- vector(mode = "list", length = ninputs)
        .Object@transformations <- vector(mode = "list", length = ninputs)
        .Object@is_nominal <- vector(mode = "logical", length = ninputs)
        .Object@is_ordinal <- vector(mode = "logical", length = ninputs)
        .Object@ordering <- vector(mode = "list", length = ninputs)
        .Object@levels <- vector(mode = "list", length = ninputs)
        .Object@scores <- vector(mode = "list", length = ninputs)
        .Object@has_missings <- vector(mode = "logical", length = ninputs)
        .Object@whichNA <- vector(mode = "list", length = ninputs)
        # .Object@jointtransf <- matrix()
        .Object
    }
)
