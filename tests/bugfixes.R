
set.seed(290875)
library("party")

### get rid of the NAMESPACE
attach(asNamespace("party"))

### check if doxygen documentation is there
stopifnot(nchar(system.file("documentation/html/index.html", package = "party")) > 29)

### check nominal level printing
set.seed(290875)
x <- gl(5, 50)
df <- data.frame(y = c(rnorm(50, 0), rnorm(50, 1), rnorm(50, 2), rnorm(50, 3), rnorm(50, 4)), 
                 x = x, z = rnorm(250))
ctree(y ~ x, data = df)

### check asymptotic vs. MonteCarlo, especially categorical splits after
### MonteCarlo resampling
a <- ctree(y ~ x + z, data = df, control = ctree_control(stump = TRUE))
b <- ctree(y ~ x + z, data = df,
           control =  ctree_control(testtype = "Monte", stump = TRUE))
stopifnot(isequal(a@tree$psplit, b@tree$psplit))  
stopifnot(isequal(a@tree$criterion$statistic, b@tree$criterion$statistic))

### we did not check the hyper parameters
try(cforest_control(minsplit = -1))
try(cforest_control(ntree = -1))
try(cforest_control(maxdepth = -1))
try(cforest_control(nresample = 10))

### NA handling for factors and in random forest
### more than one (ordinal) response variable
xo <- ordered(x)
x[sample(1:length(x), 10)] <- NA
cforest(y + xo ~ x + z, data = df, 
        control = cforest_control(ntree = 50))

### make sure minsplit is OK in the presence of missing values
### spotted by Han Lee <Han.Lee@GeodeCapital.com>
load("t1.RData")
tr <- try(ctree(p ~., data = t1))
stopifnot(!inherits(tr, "try-error"))

### make sure number of surrogate splits exceeds number of inputs by 1
### spotted by Henric Nilsson <henric.nilsson@phadia.com>
airq <- subset(airquality, !is.na(Ozone))
tr <- try(ctree(Ozone ~ Wind, data = airq,
          controls = ctree_control(maxsurrogate = 3)))
stopifnot(inherits(tr, "try-error"))

### ctree() used only the first of a multivariate response
### spotted by Henric Nilsson <henric.nilsson@phadia.com>
airq <- subset(airquality, complete.cases(Ozone, Solar.R))
airOzoSol1 <- ctree(Ozone + Solar.R ~ Wind + Temp + Month + Day,
                    data = airq)
airOzoSol2 <- ctree(Solar.R + Ozone ~ Wind + Temp + Month + Day,
                    data = airq)
stopifnot(isequal(airOzoSol1@where, airOzoSol2@where))

### one variable with all values missing
dat <- data.frame(y = rnorm(100), x1 = runif(100), x2 = rep(NA, 100))
ctree(y ~ x1 + x2, data = dat)

### one factor with only one level
dat$x2 <- factor(rep(0, 100))
try(ctree(y ~ x1 + x2, data = dat))

### weights for sampling without replacement for cforest
### spotted by Carolin Strobl <carolin.strol@stat.uni-muenchen.de>
airq <- subset(airquality, !is.na(Ozone))
cctrl <- cforest_control(replace = FALSE, fraction = 0.5)
n <- nrow(airq)
w <- double(n)


if (FALSE) {
### forest objects have weights remove in 0.9-13

### case weights
x <- runif(w)
w[x > 0.5] <- 1
w[x > 0.9] <- 2

rf <- cforest(Ozone ~ .,data = airq, weights = w, control = cctrl)
rfw <- sapply(rf@ensemble, function(x) x[[2]])
stopifnot(all(colSums(rfw) == ceiling(sum(w) / 2)))
stopifnot(max(abs(rfw[w == 0,])) == 0)

### real weights
w <- runif(n)
w[1:10] <- 0
rf <- cforest(Ozone ~ .,data = airq, weights = w, control = cctrl)
rfw <- sapply(rf@ensemble, function(x) x[[2]])
stopifnot(all(colSums(rfw) == ceiling(sum(w > 0) / 2)))
stopifnot(max(abs(rfw[w == 0,])) == 0)
}

### cforest with multivariate response
df <- data.frame(y1 = rnorm(100), y2 = rnorm(100), x1 = runif(100), x2 = runif(100))
df$y1[df$x1 < 0.5] <- df$y1[df$x1 < 0.5] + 1
cf <- cforest(y1 + y2 ~ x1 + x2, data = df)
pr <- predict(cf)
stopifnot(nrow(pr) == nrow(df) || ncol(pr) != 2)

### varimp with ordered response
### spotted by Max Kuhn <Max.Kuhn@pfizer.com>
test <- cforest(ME ~ ., data = mammoexp, control = cforest_control(ntree = 50))
stopifnot(sum(abs(varimp(test)[,1])) > 0)

### missing values in factors lead to segfaults on 64 bit systems
### spotted by Carolin Strobl <carolin.strobl@lme.de>
y <- rnorm(100)
x <- gl(2, 50)
z <- gl(2, 50)[sample(1:100)]
y <- y + (x == "1") * 3
xNA <- x
xNA[1:2] <- NA
ctree(y ~ xNA )


y <- rnorm(100)
x <- y + rnorm(100, sd = 0.1)

tmp <- data.frame(x, y)

x[sample(1:100)[1:10]] <- NA

ct1 <- ctree(y ~ x, data = tmp)
ct2 <- ctree(y ~ x, data = tmp[complete.cases(tmp),])
w <- as.double(complete.cases(tmp))
ct3 <- ctree(y ~ x, data = tmp, weights = w)

xx <- data.frame(x = rnorm(100))
t1 <- max(abs(predict(ct2, newdata = xx) - predict(ct3, newdata = xx))) == 0
t2 <- nterminal(ct1@tree) == nterminal(ct2@tree)
t3 <- nterminal(ct3@tree) == nterminal(ct1@tree)
t4 <- all.equal(ct2@tree$psplit, ct1@tree$psplit)
stopifnot(t1 && t2 && t3 && t4)

y <- rnorm(100)
x <- cut(y, c(-Inf, -1, 0, 1, Inf))

tmp <- data.frame(x, y)

x[sample(1:100)[1:10]] <- NA

ct1 <- ctree(y ~ x, data = tmp)
ct2 <- ctree(y ~ x, data = tmp[complete.cases(tmp),])
w <- as.double(complete.cases(tmp))
ct3 <- ctree(y ~ x, data = tmp, weights = w)

stopifnot(all.equal(ct2@tree$psplit, ct1@tree$psplit))
stopifnot(all.equal(ct2@tree$psplit, ct3@tree$psplit))
