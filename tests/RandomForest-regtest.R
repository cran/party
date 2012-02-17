
set.seed(290875)
library("party")
if (!require("ipred"))
    stop("cannot load package ipred")
if (!require("coin"))
    stop("cannot load package coin")

data("GlaucomaM", package = "ipred")
rf <- cforest(Class ~ ., data = GlaucomaM, control = cforest_unbiased(ntree = 100))
stopifnot(mean(GlaucomaM$Class != predict(rf)) < 
          mean(GlaucomaM$Class != predict(rf, OOB = TRUE)))

data("GBSG2", package = "ipred")
rfS <- cforest(Surv(time, cens) ~ ., data = GBSG2, control = cforest_unbiased(ntree = 100))
treeresponse(rfS, newdata = GBSG2[1:2,])

### give it a try, at least
varimp(rf, pre1.0_0 = TRUE)

P <- proximity(rf)
stopifnot(max(abs(P - t(P))) == 0)

P[1:10,1:10]

### variable importances
a <- cforest(Species ~ ., data = iris,
             control = cforest_unbiased(mtry = 2, ntree = 10))
varimp(a, pre1.0_0 = TRUE)
varimp(a, conditional = TRUE)

airq <- subset(airquality, complete.cases(airquality))
a <- cforest(Ozone ~ ., data = airq,
             control = cforest_unbiased(mtry = 2, ntree = 10))
varimp(a, pre1.0_0 = TRUE)   
varimp(a, conditional = TRUE)

a <- cforest(ME ~ ., data = mammoexp, control = cforest_classical(ntree = 10))
varimp(a, pre1.0_0 = TRUE)   
varimp(a, conditional = TRUE)

stopifnot(all.equal(unique(sapply(a@weights, sum)), nrow(mammoexp)))

