
# $Id: Predict.R,v 1.7 2005/06/22 09:32:30 hothorn Exp $

predict.BinaryTree <- function(object, ...) {
    conditionalTree@predict(object, ...)
}

setGeneric("treeresponse", function(object, ...) 
           standardGeneric("treeresponse"))

setMethod("treeresponse", signature = signature(object = "BinaryTree"),
    definition = function(object, newdata = NULL, ...)   
        object@cond_distr_response(newdata = newdata, ...)
)


setGeneric("weights", function(object, ...) standardGeneric("weights"))

setMethod("weights", signature = signature(object = "BinaryTree"),
    definition = function(object, newdata = NULL, ...)
        object@prediction_weights(newdata = newdata, ...)
)


setGeneric("where", function(object, ...) standardGeneric("where"))

setMethod("where", signature = signature(object = "BinaryTree"),
    definition = function(object, newdata = NULL, ...)
        object@get_where(newdata = newdata, ...)
)


setGeneric("nodes", function(object, where, ...) standardGeneric("nodes"))

setMethod("nodes", signature = signature(object = "BinaryTree", 
                                         where = "integer"),
    definition = function(object, where, ...)
        lapply(where, function(i) .Call("R_get_nodebynum", object@tree, i))
)

setMethod("nodes", signature = signature(object = "BinaryTree", 
                                         where = "numeric"),
    definition = function(object, where, ...)
        nodes(object, as.integer(where))
)
