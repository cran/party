
# $Id: Memory.R,v 1.5 2005/06/23 08:28:42 hothorn Exp $

ctree_memory <- function(object, MPinv = FALSE) {

    inputs <- object@inputs
    responses <- object@responses
    q <- ncol(responses@jointtransf)
    nobs <- inputs@nobs
    ninputs <- inputs@ninputs
    RET <- new("TreeFitMemory")
    expcovinf <- new("ExpectCovarInfluence", q)
    expcovinfss <- new("ExpectCovarInfluence", as.integer(1))
    linexpcov2sample  <- new("LinStatExpectCovar", as.integer(1), q) 
    weights <- vector(mode = "list", length = ninputs)
    splitstatistics <- as.double(rep(0, nobs))
    varmemory <- vector(mode = "list", length = ninputs)
    varMmemory <- vector(mode = "list", length = ninputs)
    Mscorematrices <- vector(mode = "list", length = ninputs)
    dontuse <- rep(FALSE, ninputs)
    dontusetmp <- rep(FALSE, ninputs)

    yORDINAL <- (responses@is_ordinal[1] && responses@ninputs == 1)

    for (j in 1:inputs@ninputs) {

        xORDINAL <- inputs@is_ordinal[j]
        p <- ncol(inputs@transformations[[j]])

        if (inputs@has_missings[j])
            weights[[j]] <- rep(0, nobs)

        if (xORDINAL && !yORDINAL)
            Mscorematrices[[j]] <- .Call("R_scmatleft",
                inputs@scores[[j]], p * q, PACKAGE = "party")
        if (xORDINAL && yORDINAL)
            ### grrr: class(kronecker(1:4, 1:3)) == "array"
            Mscorematrices[[j]] <- matrix(kronecker(responses@scores[[1]], 
                inputs@scores[[j]]), nrow = 1)

        if (!xORDINAL && yORDINAL)
            Mscorematrices[[j]] <- .Call("R_scmatright", 
                responses@scores[[1]], p * q, PACKAGE = "party")
        
        if (!xORDINAL && !yORDINAL) {
            if (MPinv) {
                varmemory[[j]] <- new("LinStatExpectCovarMPinv", p, q)
            } else {
                varmemory[[j]] <- new("LinStatExpectCovar", p, q)
            }
        } else {
            varmemory[[j]] <- new("LinStatExpectCovar", p, q)
            Mp <- nrow(Mscorematrices[[j]])
            if (MPinv) {
                varMmemory[[j]] <- new("LinStatExpectCovarMPinv", Mp,
                                           as.integer(1))
            } else {
                varMmemory[[j]] <- new("LinStatExpectCovar", Mp,
                                          as.integer(1))
            }
        }
    }

    RET@expcovinf <- expcovinf
    RET@expcovinfss <- expcovinfss
    RET@linexpcov2sample  <- linexpcov2sample
    RET@weights <- weights
    RET@varmemory <- varmemory
    RET@varMmemory <- varMmemory
    RET@Mscorematrices <- Mscorematrices
    RET@dontuse <- dontuse
    RET@dontusetmp <- dontusetmp
    RET@splitstatistics <- splitstatistics
    RET
}
