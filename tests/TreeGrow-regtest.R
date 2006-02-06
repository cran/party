
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(ipred))
    stop("cannot load package ipred")
if (!require(coin))
    stop("cannot load package coin")
gctorture(on = GCtorture)

### get rid of the NAMESPACE
nsparty <- attach(NULL, name="ns-party")
.Internal(lib.fixup(asNamespace("party"), nsparty))

gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)

data(GlaucomaM, package = "ipred")
gtree <- ctree(Class ~ ., data = GlaucomaM)
tree <- gtree@tree
stopifnot(isequal(tree[[5]][[3]], 0.059))
predict(gtree)

# print(tree)

stump <- ctree(Class ~ ., data = GlaucomaM, 
               control = ctree_control(stump = TRUE))
print(stump)

data(treepipit, package = "coin")

tr <- ctree(counts ~ ., data = treepipit)
tr
plot(tr)


data(GlaucomaM, package = "ipred")

tr <- ctree(Class ~ ., data = GlaucomaM)
tr
plot(tr)

data(GBSG2, package = "ipred")  

GBSG2tree <- ctree(Surv(time, cens) ~ ., data = GBSG2)
GBSG2tree
plot(GBSG2tree)
plot(GBSG2tree, terminal_panel = node_surv(GBSG2tree))
survfit(Surv(time, cens) ~ as.factor(GBSG2tree@where), data = GBSG2)
names(GBSG2)

tr <- ctree(Surv(time, cens) ~ ., data = GBSG2, 
            control = ctree_control(teststattype = "maxabs", 
                                    testtype = "Raw"))
tr
plot(tr)

data(mammoexp)
attr(mammoexp$ME, "scores") <- 1:3   
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
names(mammoexp)[names(mammoexp) == "SYMPT"] <- "symptoms"
names(mammoexp)[names(mammoexp) == "PB"] <- "benefit"

names(mammoexp)
tr <- ctree(ME ~ ., data = mammoexp)
tr
plot(tr)

