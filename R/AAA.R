
# $Id: AAA.R 2466 2006-02-14 18:55:44Z hothorn $

.onLoad <- function(lib, pkg) {
    if (!require("methods")) stop("cannot load methods")
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
