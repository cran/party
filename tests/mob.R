library("party")

data("BostonHousing", package = "mlbench")
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)
fmBH <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  control = mob_control(minsplit = 40, verbose = TRUE),
  data = BostonHousing, model = linearModel)
fmBH
summary(fmBH)

data("PimaIndiansDiabetes", package = "mlbench")
fmPID <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age,
  control = mob_control(verbose = TRUE),
  data = PimaIndiansDiabetes, model = glinearModel, family = binomial())
fmPID
summary(fmPID)
