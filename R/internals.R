
party_intern <- function(..., fun = c("R_TreeGrow", "R_get_nodeID",
    "R_getpredictions", "initVariableFrame", "ctreedpp", "newinputs")) {

    fun <- match.arg(fun)
    do.call(fun, list(...))
}
