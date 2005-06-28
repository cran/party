
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
                 x = x)
ctree(y ~ x, data = df)
