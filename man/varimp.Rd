\name{varimp}
\alias{varimp}
\title{ Variable Importance }
\description{
    Standard and conditional variable importance for `cforest', following the permutation
    principle of the `mean decrease in accuracy' importance in `randomForest'.
}
\usage{
varimp(object, mincriterion = 0, conditional = FALSE, 
       threshold = 0.2, nperm = 1, OOB = TRUE)
}
\arguments{
  \item{object}{ an object as returned by \code{cforest}.}
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to include a split in the 
                       computation of the importance. The default \code{mincriterion = 0}
                       guarantees that all splits are included.}
  \item{conditional}{ a logical determining whether unconditional or conditional 
                      computation of the importance is performed. }
  \item{threshold}{ the value of the test statistic or 1 - p-value of the association 
                    between the variable of interest and a covariate that must be 
                    exceeded inorder to include the covariate in the conditioning 
                    scheme for the variable of interest (only relevant if 
                    \code{conditional = TRUE}). }
  \item{nperm}{ the number of permutations performed.}
  \item{OOB}{ a logical determining whether the importance is computed from the out-of-bag 
              sample or the learning sample (not suggested).}
}
\details{

  Function \code{varimp} can be used to compute variable importance measures
  similar to those computed by \code{\link[randomForest]{importance}}. Besides the
  standard version, a conditional version is available, that adjusts for correlations between
  predictor variables. 
  
  If \code{conditional = TRUE}, the importance of each variable is computed by permuting 
  within a grid defined by the covariates that are associated  (with 1 - p-value 
  greater than \code{threshold}) to the variable of interest.
  The resulting variable importance score is conditional in the sense of beta coefficients in   
  regression models, but represents the effect of a variable in both main effects and interactions.
  See Strobl et al. (2008) for details.

  Note, however, 

  1. that all random forest results are subject to random variation. Thus, before
  interpreting the importance ranking, check whether the same ranking is achieved with a
  different random seed -- or otherwise increase the number of trees \code{ntree} in 
  \code{\link{ctree_control}}!

  2. that \code{varimp} cannot (yet) handle missing values.
  

}
\value{
  A vector of `mean decrease in accuracy' importance scores.
}
\references{ 

    Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5--32.

    Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006b). Unbiased
    Recursive Partitioning: A Conditional Inference Framework.
    \emph{Journal of Computational and Graphical Statistics}, \bold{15} (3),
    651-674.  Preprint available from 
    \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

    Carolin Strobl, Anne-Laure Boulesteix, Thomas Kneib, Thomas Augustin and Achim Zeileis (2008).
    Conditional Variable Importance for Random Forests. \emph{BMC Bioinformatics}, \bold{9}, 307. 
    \url{http://www.biomedcentral.com/1471-2105/9/307}
}
\examples{
    
   set.seed(290875)
   readingSkills.cf <- cforest(score ~ ., data = readingSkills, 
       control = cforest_unbiased(mtry = 2, ntree = 50))

   # standard importance
   varimp(readingSkills.cf)

   # conditional importance, may take a while...
   varimp(readingSkills.cf, conditional = TRUE)

}
\keyword{tree}