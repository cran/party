###################################################
### chunk number 1: setup
###################################################
options(width = 65, SweaveHooks = list(leftpar = 
    function() par(mai = par("mai") * c(1, 1.1, 1, 1))))
require("party")
require("coin")
set.seed(290875)


###################################################
### chunk number 2: bibtex eval=FALSE
###################################################
## bib <- HSAuR:::readBibtex("partyrefs.bib")
## subbib <- HSAuR:::extractBibtex("party.Rnw", bib)
## sink("party.bib")
## HSAuR:::toBibtex.txtBibtex(subbib)
## sink()


###################################################
### chunk number 3: party-data
###################################################
ls <- data.frame(y = gl(3, 50, labels = c("A", "B", "C")), x1 = rnorm(150) + rep(c(1, 0, 0), c(50, 50,
50)), x2 = runif(150))


###################################################
### chunk number 4: party-formula
###################################################
library("party")
ctree(y ~ x1 + x2, data = ls)


###################################################
### chunk number 5: party-xtrafo
###################################################
ctree(y ~ x1 + x2, data = ls, xtrafo = function(data) trafo(data,
numeric_trafo = rank))


###################################################
### chunk number 6: party-Bonf
###################################################
ctree_control(testtype = "Bonferroni")


###################################################
### chunk number 7: party-MC
###################################################
ctree_control(testtype = "MonteCarlo")           


###################################################
### chunk number 8: party-ss
###################################################
ctree_control(savesplitstats = TRUE)


###################################################
### chunk number 9: party-minsplit
###################################################
ctree_control(minsplit = 20)


###################################################
### chunk number 10: party-maxsurr
###################################################
ctree_control(maxsurrogate = 3)


###################################################
### chunk number 11: party-fitted
###################################################
ct <- ctree(y ~ x1 + x2, data = ls)


###################################################
### chunk number 12: party-print
###################################################
ct


###################################################
### chunk number 13: party-plot
###################################################
plot(ct)


###################################################
### chunk number 14: party-nodes
###################################################
nodes(ct, 1)


###################################################
### chunk number 15: party-nodelist
###################################################
names(nodes(ct, 1)[[1]])


###################################################
### chunk number 16: party-Predict
###################################################
Predict(ct, newdata = ls)


###################################################
### chunk number 17: party-treeresponse
###################################################
treeresponse(ct, newdata = ls[c(1, 51, 101),])


###################################################
### chunk number 18: party-where
###################################################
where(ct, newdata = ls[c(1,51,101),])


###################################################
### chunk number 19: treepipit-ctree
###################################################
data("treepipit", package = "coin")
tptree <- ctree(counts ~ ., data = treepipit)


###################################################
### chunk number 20: treepipit-plot
###################################################
plot(tptree, terminal_panel = node_hist(tptree, breaks = 0:6-0.5, ymax = 65, 
horizontal = FALSE, freq = TRUE))


###################################################
### chunk number 21: treepipit-x
###################################################
x <- tptree@tree


###################################################
### chunk number 22: glaucoma-ctree
###################################################
data("GlaucomaM", package = "ipred")
gtree <- ctree(Class ~ ., data = GlaucomaM)


###################################################
### chunk number 23: glaucoma-x
###################################################
x <- gtree@tree


###################################################
### chunk number 24: glaucoma-plot
###################################################
plot(gtree)


###################################################
### chunk number 25: glaucoma-plot-inner
###################################################
plot(gtree, inner_panel = node_barplot, 
     edge_panel = function(ctreeobj, ...) { function(...) invisible() }, tnex = 1)


###################################################
### chunk number 26: glaucoma-split-plot
###################################################
cex <- 1.6
inner <- nodes(gtree, 1:3)
layout(matrix(1:length(inner), ncol = length(inner)))
out <- sapply(inner, function(i) {
    splitstat <- i$psplit$splitstatistic
    x <- GlaucomaM[[i$psplit$variableName]][splitstat > 0]
    plot(x, splitstat[splitstat > 0], main = paste("Node", i$nodeID),
         xlab = i$psplit$variableName, ylab = "Statistic", ylim = c(0, 10),
         cex.axis = cex, cex.lab = cex, cex.main = cex)
    abline(v = i$psplit$splitpoint, lty = 3)
})


###################################################
### chunk number 27: glaucoma-prediction
###################################################
table(Predict(gtree), GlaucomaM$Class)


###################################################
### chunk number 28: glaucoma-classprob
###################################################
prob <- sapply(treeresponse(gtree), function(x) x[1]) +
               runif(nrow(GlaucomaM), min = -0.01, max = 0.01)
splitvar <- nodes(gtree, 1)[[1]]$psplit$variableName
plot(GlaucomaM[[splitvar]], prob, 
     pch = as.numeric(GlaucomaM$Class), ylab = "Conditional Class Prob.",
     xlab = splitvar)
abline(v = nodes(gtree, 1)[[1]]$psplit$splitpoint, lty = 2)
legend(0.15, 0.7, pch = 1:2, legend = levels(GlaucomaM$Class), bty = "n")


###################################################
### chunk number 29: GBSGS-ctree
###################################################
data("GBSG2", package = "ipred")  
stree <- ctree(Surv(time, cens) ~ ., data = GBSG2)


###################################################
### chunk number 30: GBSG2-plot
###################################################
plot(stree)


###################################################
### chunk number 31: GBSG2-KM
###################################################
treeresponse(stree, newdata = GBSG2[1:2,])


###################################################
### chunk number 32: mammo-ctree
###################################################
data("mammoexp", package = "party")
mtree <- ctree(ME ~ ., data = mammoexp)


###################################################
### chunk number 33: mammo-plot
###################################################
plot(mtree)


###################################################
### chunk number 34: spider-ctree eval=FALSE
###################################################
## data("spider", package = "mvpart")   
## sptree <- ctree(arct.lute + pard.lugu + zora.spin + pard.nigr + pard.pull +
##            aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce + 
##            alop.fabr + arct.peri ~ herbs+reft+moss+sand+twigs+water, 
##            controls = ctree_control(teststat = "maxabs", testtype = "MonteCarlo", nresample = 9999, 
##            minsplit = 5, mincriterion = 0.9), data = spider)
## sptree@tree$criterion
## 
## library("coin")
## data("spider", package = "mvpart")
## it <- independence_test(arct.lute + pard.lugu + zora.spin + pard.nigr + pard.pull +
##                        aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce +   
##                        alop.fabr + arct.peri ~ herbs+reft+moss+sand+twigs+water, 
##                        data = spider, distribution = approximate(B = 19999))
## statistic(it, "standardized")
## pvalue(it)


###################################################
### chunk number 35: spider-plot eval=FALSE
###################################################
## plot(sptree, terminal_panel = node_terminal)


