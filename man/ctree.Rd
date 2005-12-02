\name{Conditional Trees}
\alias{ctree}
\alias{conditionalTree}
\title{ Conditional Trees }
\description{
  Recursive partitioning for continuous, censored, ordered, nominal and
  multivariate response variables in a conditional inference framework. 
}
\usage{
ctree(formula, data, subset = NULL, weights = NULL, 
      controls = ctree_control(), xtrafo = ptrafo, ytrafo = ptrafo, 
      scores = NULL)
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
}
\details{

  Conditional trees estimate a regression relationship by binary recursive
  partitioning in a conditional inference framework. Roughly, the algorithm
  works as follows: 1) Test the global null hypothesis of independence between
  any of the input variables and the response (which may be multivariate as well). 
  Stop if this hypothesis cannot be rejected. Otherwise select the input
  variable with strongest association to the resonse. This
  association is measured by a p-value corresponding to a test for the
  partial null hypothesis of a single input variable and the response.
  2) Implement a binary split in the selected input variable. 
  3) Recursively repeate steps 1) and 2). 

  The implementation utilizes a unified framework for conditional inference,
  or permutation tests, developed by Strasser and Weber (1999). The stop
  criterion in step 1) is either based on a p-value 
  (\code{teststattype = "Bonferroni"} 
  or \code{teststattype = "MonteCarlo"} in \code{\link{ctree_control}})
  or on the raw (standardized) test
  statistic (\code{teststattype = "Raw"}). In both cases, the
  criterion is maximized, i.e., 1 - p-value is used. A split is implemented 
  when the criterion exceeds the value given by \code{mincriterion} as
  specified in \code{\link{ctree_control}}. For example, when 
  \code{mincriterion = 0.95}, the p-value must be smaller than
  $0.05$ in order to split this node. This statistical approach ensures that
  the right sized tree is grown and no form of pruning or cross-validation
  or whatsoever is needed. The selection of the input variable to split in
  is based on the univariate p-values avoiding a variable selection bias
  towards input variables with many possible cutpoints.

  By default, the scores for each ordinal factor \code{x} are
  \code{1:length(x)}, this may be changed using \code{scores = list(x =
  c(1,5,6))}, for example.

  For a general description of the methodology see Hothorn, Hornik and
  Zeileis (2004).

}
\value{
  An object of class \code{BinaryTree}.
}
\references{ 

   Torsten Hothorn, Kurt Hornik and Achim Zeileis (2004). Unbiased Recursive
   Partitioning: A Conditional Inference Framework. Technical Report Nr. 8, 
   Research Report Series / Department of Statistics and Mathematics, WU Wien.
   \url{http://epub.wu-wien.ac.at/dyn/openURL?id=oai:epub.wu-wien.ac.at:epub-wu-01_756}

   Helmut Strasser and Christian Weber (1999). On the asymptotic theory of permutation
   statistics. \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.

}
\examples{

    ### regression
    data("airquality")
    airq <- subset(airquality, !is.na(Ozone))
    airct <- ctree(Ozone ~ ., data = airq, 
                   controls = ctree_control(maxsurrogate = 3))
    airct
    plot(airct)
    mean((airq$Ozone - predict(airct))^2)

    ### classification
    irisct <- ctree(Species ~ .,data = iris)
    irisct
    plot(irisct)
    table(predict(irisct), iris$Species)

    ### estimated class probabilities, a list
    tr <- treeresponse(irisct, newdata = iris[1:10,])

    ### ordinal regression
    data("mammoexp", package = "party") 
    mammoct <- ctree(ME ~ ., data = mammoexp) 
    plot(mammoct)

    ### estimated class probabilities
    treeresponse(mammoct, newdata = mammoexp[1:10,])

    ### survival analysis
    if (require("ipred")) {
        data("GBSG2", package = "ipred")
        GBSG2ct <- ctree(Surv(time, cens) ~ .,data = GBSG2)
        plot(GBSG2ct)
        treeresponse(GBSG2ct, newdata = GBSG2[1:2,])        
    }

}
\keyword{tree}
