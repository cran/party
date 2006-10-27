
library("party")

### get rid of the NAMESPACE
nsparty <- attach(NULL, name="ns-party")
.Internal(lib.fixup(asNamespace("party"), nsparty))

splitctrl <- new("SplitControl")
splitctrl@minprob <- 0.1
splitctrl@minsplit <- as.integer(5)
xf <- gl(5, 100/5)
weights <- rep(1, length(xf))
n <- length(xf)

ret <- rep(TRUE, 100)
for (i in 1:length(ret)) {
    yf <- gl(4, n/4)[sample(1:length(xf))]
    # yf <- rnorm(n)
    split <- Split(xf, yf, weights, splitctrl)
    lev1 <- levels(xf)[!as.logical(split[[4]])]
    lev2 <- show(maxstat_test(yf ~ xf))$estimate[[1]]

    ret[i] <- (all(lev1 %in% lev2) && all(lev2 %in% lev1)) || 
               all(c(lev1, lev2) %in% levels(xf))
}

all(ret)

