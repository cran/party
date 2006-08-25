
set.seed(290875)
gctorture(on = FALSE)
library(party)
gctorture(on = GCtorture)

### get rid of the NAMESPACE
nsparty <- attach(NULL, name="ns-party")
.Internal(lib.fixup(asNamespace("party"), nsparty))

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
