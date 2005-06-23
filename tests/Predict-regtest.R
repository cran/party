
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(ipred))
    stop("cannot load package ipred")
if (!require(coin))
    stop("cannot load package coin")
gctorture(on = GCtorture)

### get rid of the NAMESPACE
load(file.path(.find.package("party"), "R", "all.rda"))

data(treepipit, package = "coin")
ct <- ctree(counts ~ ., data = treepipit)
stopifnot(isequal(predict(ct), predict(ct, newdata = treepipit)))


data(GlaucomaM, package = "ipred")
ct <- ctree(Class ~ ., data = GlaucomaM)
stopifnot(isequal(predict(ct), predict(ct, newdata = GlaucomaM)))
stopifnot(isequal(predict(ct, type = "prob"), predict(ct, type = "prob", 
                  newdata = GlaucomaM)))

data(GBSG2, package = "ipred")  

GBSG2tree <- ctree(Surv(time, cens) ~ ., data = GBSG2)
stopifnot(isequal(GBSG2tree@predict_response(), 
          GBSG2tree@predict_response(newdata = GBSG2)))
stopifnot(isequal(GBSG2tree@cond_distr_response(), 
          GBSG2tree@cond_distr_response(newdata = GBSG2)))

data(mammoexp)
attr(mammoexp$ME, "scores") <- 1:3   
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
names(mammoexp)[names(mammoexp) == "SYMPT"] <- "symptoms"
names(mammoexp)[names(mammoexp) == "PB"] <- "benefit"

names(mammoexp)
mtree <- ctree(ME ~ ., data = mammoexp)
stopifnot(isequal(predict(mtree), predict(mtree, newdata = mammoexp)))
stopifnot(isequal(predict(mtree), predict(mtree, newdata = mammoexp)))
