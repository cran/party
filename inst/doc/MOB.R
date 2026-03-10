### R code from vignette source 'MOB.Rnw'

###################################################
### code chunk number 1: setup
###################################################
require("party")
options(useFancyQuotes = FALSE)


###################################################
### code chunk number 2: MOB.Rnw:229-230
###################################################
library("party")


###################################################
### code chunk number 3: MOB.Rnw:235-236
###################################################
data("BostonHousing", package = "mlbench")


###################################################
### code chunk number 4: MOB.Rnw:250-252
###################################################
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2


###################################################
### code chunk number 5: MOB.Rnw:267-269
###################################################
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)


###################################################
### code chunk number 6: MOB.Rnw:280-282
###################################################
ctrl <- mob_control(alpha = 0.05, bonferroni = TRUE, minsplit = 40,
  objfun = deviance, verbose = TRUE)


###################################################
### code chunk number 7: MOB.Rnw:295-297
###################################################
fmBH <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, control = ctrl, model = linearModel)


###################################################
### code chunk number 8: MOB.Rnw:306-307
###################################################
fmBH


###################################################
### code chunk number 9: MOB.Rnw:313-314 (eval = FALSE)
###################################################
## plot(fmBH)


###################################################
### code chunk number 10: BostonHousing-plot
###################################################
plot(fmBH)


###################################################
### code chunk number 11: MOB.Rnw:347-348
###################################################
coef(fmBH)


###################################################
### code chunk number 12: MOB.Rnw:354-355
###################################################
summary(fmBH, node = 7)


###################################################
### code chunk number 13: MOB.Rnw:362-363
###################################################
sctest(fmBH, node = 7)


###################################################
### code chunk number 14: MOB.Rnw:369-372
###################################################
mean(residuals(fmBH)^2)
logLik(fmBH)
AIC(fmBH)


###################################################
### code chunk number 15: MOB.Rnw:375-377
###################################################
nt <- NROW(coef(fmBH))
nk <- NCOL(coef(fmBH))


###################################################
### code chunk number 16: MOB.Rnw:393-395
###################################################
data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes <- na.omit(PimaIndiansDiabetes2[,-c(4, 5)])


###################################################
### code chunk number 17: MOB.Rnw:415-417
###################################################
fmPID <- mob(diabetes ~ glucose | pregnant + pressure + mass + pedigree + age,
  data = PimaIndiansDiabetes, model = glinearModel, family = binomial())


###################################################
### code chunk number 18: MOB.Rnw:422-423 (eval = FALSE)
###################################################
## plot(fmPID)


###################################################
### code chunk number 19: PimaIndiansDiabetes-plot
###################################################
plot(fmPID)


###################################################
### code chunk number 20: MOB.Rnw:457-459
###################################################
coef(fmPID)
exp(coef(fmPID)[,2])


###################################################
### code chunk number 21: MOB.Rnw:462-463
###################################################
risk <- round(100 * (exp(coef(fmPID)[,2])-1), digits = 1)


