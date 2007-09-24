\name{Control Forest Hyper Parameters}
\alias{cforest_control}
\alias{cforest_classical}
\alias{cforest_unbiased}
\title{ Control for Conditional Tree Forests }
\description{

  Various parameters that control aspects of the `cforest' fit via
  its `control' argument.

}
\usage{
cforest_unbiased(\dots)
cforest_classical(\dots)
cforest_control(teststat = "max",
                testtype = "Teststatistic",
                mincriterion = qnorm(0.9),
                savesplitstats = FALSE,
                ntree = 500, mtry = 5, replace = TRUE,
                fraction = 0.632, \dots)
}
\arguments{
  \item{teststat}{ a character specifying the type of the test statistic
                       to be applied. }
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. }
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to implement a split. }
  \item{mtry}{ number of input variables randomly sampled as candidates 
               at each node for random forest like algorithms. The default
               \code{mtry = 0} means that no random selection takes place
               and bagging is performed.}
  \item{savesplitstats}{ a logical determining if the process of standardized
                         two-sample statistics for split point estimate
                         is saved for each primary split.}
  \item{ntree}{ number of trees to grow in a forest.}
  \item{replace}{ a logical indicating whether sampling of observations is 
                 done with or without replacement.}
  \item{fraction}{ fraction of number of observations to draw without 
                   replacement (only relevant if \code{replace = TRUE}).}
  \item{\dots}{ additional arguments to be passed to 
                \code{\link{ctree_control}}.}
}
\details{

  All three functions return an object of class \code{\link{ForestControl-class}}
  defining hyper parameters to be specified via the \code{control} argument
  of \code{\link{cforest}}.

  The arguments \code{teststat}, \code{testtype} and \code{mincriterion}
  determine how the global null hypothesis of independence between all input
  variables and the response is tested (see \code{\link{ctree}}). The 
  argument \code{nresample} is the number of Monte-Carlo replications to be
  used when \code{testtype = "MonteCarlo"}.

  A split is established when the sum of the weights in both daugther nodes
  is larger than \code{minsplit}, this avoids pathological splits at the
  borders. When \code{stump = TRUE}, a tree with at most two terminal nodes
  is computed.

  The argument \code{mtry > 0} means that a random forest like `variable
  selection', i.e., a random selection of \code{mtry} input variables, is
  performed in each node.

  It might be informative to look at scatterplots of input variables against
  the standardized two-sample split statistics, those are available when
  \code{savesplitstats = TRUE}. Each node is then associated with a vector
  those length is determined by the number of observations in the learning
  sample and thus much more memory is required.

  Function \code{cforest_unbiased} returns the settings suggested 
  for the construction of unbiased random forests by Strobl et al. (2007)
  and is the default since version 0.9-90.
  Hyper parameter settings mimicing the behaviour of
  \code{\link[randomForest]{randomForest}} are available in
  \code{cforest_classical} which have been used as default up to
  version 0.9-14. 

}
\value{
  An object of class \code{\link{ForestControl-class}}.
}
\references{

    Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007).
    Bias in Random Forest Variable Importance Measures: Illustrations, Sources and  
    a Solution. \emph{BMC Bioinformatics}, \bold{8}, 25. 
    \url{http://www.BioMedCentral.com/1471-2105/8/25/}
}
\keyword{misc}
