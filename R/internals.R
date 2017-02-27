
party_intern <- function(..., fun = c("R_TreeGrow", "R_get_nodeID",
    "R_getpredictions", "initVariableFrame", "ctreedpp", "newinputs")) {

    fun <- match.arg(fun)
    fun <- gsub("^R", "\\.R", fun) ### needed to add . to R_xyz functions
    do.call(fun, list(...))
}
