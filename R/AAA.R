
# $Id: AAA.R 681 2023-11-27 11:55:24Z thothorn $

.onLoad <- function(lib, pkg) {
    ### GCtorture <<- FALSE
    .Call(party_init)
    return(TRUE)
}
