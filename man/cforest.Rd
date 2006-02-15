\name{cforest}
\alias{cforest}
\title{ Random Forest }
\description{
    An implementation of the random forest and bagging ensemble algorithms
    utilizing conditional inference trees as base learners.
}
\usage{
cforest(formula, data = list(), subset = NULL, weights = NULL, 
        controls = ctree_control(teststattype = "maxabs", 
            testtype = "Raw", mtry = 5, savesplitstats = FALSE), 
        xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL, B = 500)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ an data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{controls}{an object of class \code{\link{TreeControl}}, which can be
                  obtained using \code{\link{ctree_control}}.}
  \item{xtrafo}{ a function to be applied to all input variables.
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{ytrafo}{ a function to be applied to all response variables.
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{scores}{ an optional named list of scores to be attached to ordered
               factors.}
  \item{B}{number of trees to grown.}
}
\details{

  This implementation of the random forest (and bagging) algorithm differs
  from the reference implementation in \code{\link[randomForest]{randomForest}]}
  with respect to the base learner used and the aggregation scheme applied.
  
  Conditional inference trees, see \code{\link{ctree}}, are fitted to each
  of the \code{B} bootstrap samples of the learning sample. There are many
  hyper parameters that can be controlled, see \code{\link{ctree_control}}.
  You MUST NOT change anything you don't understand completely.

  The aggregation scheme works by averaging observation weights extracted
  from each of the \code{B} trees and NOT by averaging predictions directly.
  See Hothorn et al. (2004) for a description. 

  Ensembles of conditional inference trees have not yet been extensively
  tested, so this routine is meant for the expert user only and its current
  state is rather experimental. However, there are some things that can't be
  done with \code{\link[randomForest]{randomForest}}, for example fitting
  forests to censored response variables.

}
\value{
  An object of class \code{\link{RandomForest-class}}.
}
\references{ 

    Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5--32.

    Torsten Hothorn, Berthold Lausen, Axel Benner and Martin Radespiel-Troeger
    (2004). Bagging Survival Trees. \emph{Statistics in Medicine}, \bold{23}(1), 77--91.

}
\examples{

    ### honest (i.e., out-of-bag) cross-classification of
    ### true vs. predicted classes
    table(mammoexp$ME, predict(cforest(ME ~ ., data = mammoexp, B = 50),
                               OOB = TRUE))

    ### fit forest to censored response
    if (require("ipred")) {

        data("GBSG2", package = "ipred")
        bst <- cforest(Surv(time, cens) ~ ., data = GBSG2, B = 50)

        ### estimate conditional Kaplan-Meier curves
        treeresponse(bst, newdata = GBSG2[1:2,], OOB = TRUE)

        ### if you can't resist to look at individual trees ...
        party:::prettytree(bst@ensemble[[1]], names(bst@data@get("input")))
    }
}
\keyword{tree}
