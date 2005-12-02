
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(ipred))
    stop("cannot load package ipred")
gctorture(on = GCtorture)

### get rid of the NAMESPACE
load(file.path(.find.package("party"), "R", "all.rda"))

gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)

data(GlaucomaM, package = "ipred")
inp <- initVariableFrame(GlaucomaM[,-63,drop = FALSE], trafo = NULL) #, fun = rank)
resp <- initVariableFrame(GlaucomaM[,"Class",drop = FALSE], trafo = NULL)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- ctree_memory(ls, TRUE)
ctrl <- ctree_control()
node <- .Call("R_Node", ls, ls@weights, tm, ctrl)
stopifnot(isequal(node[[5]][[3]], 0.059))

### and now with ranked inputs -> Wilcoxon-Mann-Whitney tests
inp <- initVariableFrame(GlaucomaM[,-63,drop = FALSE], trafo = function(data)
ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(GlaucomaM[,"Class",drop = FALSE], trafo = NULL)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- ctree_memory(ls, TRUE)
ctrl <- ctree_control()
node <- .Call("R_Node", ls, ls@weights, tm, ctrl)
stopifnot(isequal(node[[5]][[3]], 0.059))
