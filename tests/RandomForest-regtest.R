
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(ipred))
    stop("cannot load package ipred")
if (!require(coin))
    stop("cannot load package coin")
gctorture(on = GCtorture)

data("GlaucomaM", package = "ipred")
rf <- cforest(Class ~ ., data = GlaucomaM, control = cforest_control(ntree = 100))
stopifnot(mean(GlaucomaM$Class != predict(rf)) < 
          mean(GlaucomaM$Class != predict(rf, OOB = TRUE)))

data("GBSG2", package = "ipred")
rfS <- cforest(Surv(time, cens) ~ ., data = GBSG2, control = cforest_control(ntree = 100))
treeresponse(rfS, newdata = GBSG2[1:2,])

### give it a try, at least
varimp(rf)