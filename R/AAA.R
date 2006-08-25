
# $Id: AAA.R 2466 2006-02-14 18:55:44Z hothorn $

.onLoad <- function(lib, pkg) {
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
