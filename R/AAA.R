
# $Id: AAA.R 630 2017-02-27 14:58:59Z thothorn $

.onLoad <- function(lib, pkg) {
    GCtorture <<- FALSE
    .Call(party_init)
    return(TRUE)
}
