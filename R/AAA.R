
# $Id: AAA.R,v 1.7 2005/06/16 11:09:59 hothorn Exp $

.onLoad <- function(lib, pkg) {
    if (!require("methods")) stop("cannot load methods")
    if (!require("grid")) stop("cannot load grid")
    if (!require("survival")) stop("cannot load survival")
    if (!require("modeltools")) stop("cannot load modeltools")
    # if (!require("coin")) stop("cannot load coin")
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
