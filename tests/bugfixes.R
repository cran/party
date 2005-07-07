
set.seed(290875)
gctorture(on = FALSE)
library(party)
gctorture(on = GCtorture)

### get rid of the NAMESPACE
load(file.path(.find.package("party"), "R", "all.rda"))

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
