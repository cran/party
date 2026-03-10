\name{varimp}
\alias{varimp}
\alias{varimpAUC}
\title{ Variable Importance }
\description{
    Standard and conditional variable importance for `cforest', following the permutation
    principle of the `mean decrease in accuracy' importance in `randomForest'.
}
\usage{
varimp(object, mincriterion = 0, conditional = FALSE, 
       threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
varimpAUC(...)
}
\arguments{
  \item{object}{ an object as returned by \code{cforest}.}
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to include a split in the 
                       computation of the importance. The default \code{mincriterion = 0}
                       guarantees that all splits are included.}
  \item{conditional}{ a logical determining whether unconditional or conditional 
                      computation of the importance is performed. }
  \item{threshold}{ the threshold value for (1 - p-value) of the association 
                    between the variable of interest and a covariate, which must be 
                     exceeded inorder to include the covariate in the conditioning 
                     scheme for the variable of interest (only relevant if 
                    \code{conditional = TRUE}). A threshold value of zero includes 
                    all covariates.}
  \item{nperm}{ the number of permutations performed.}
  \item{OOB}{ a logical determining whether the importance is computed from the out-of-bag 
              sample or the learning sample (not suggested).}
  \item{pre1.0_0}{ Prior to party version 1.0-0, the actual data values
                   were permuted according to the original permutation
                   importance suggested by \bibcitet{party::Breiman:2001}. Now the assignments
                   to child nodes of splits in the variable of interest
                   are permuted as described by
                   \bibcitet{party::Hapfelmeier_Hothorn_Ulm_Strobl_2014},
                   which allows for missing values in the explanatory
                   variables and is more efficient wrt memory consumption and 
                   computing time. This method does not apply to conditional
                   variable importances.}
  \item{\dots}{Arguments to \code{\link[varImp]{varImpAUC}}.}
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
  See \bibcitet{party::Strobl+Boulesteix+Kneib:2008} for details.

  Note, however, that all random forest results are subject to random variation. Thus, before
  interpreting the importance ranking, check whether the same ranking is achieved with a
  different random seed -- or otherwise increase the number of trees \code{ntree} in 
  \code{\link{ctree_control}}.

  Note that in the presence of missings in the predictor variables the procedure
  described in \bibcitet{party::Hapfelmeier_Hothorn_Ulm_Strobl_2014} is performed.

  Function \code{varimpAUC} is a wrapper for
  \code{\link[varImp]{varImpAUC}} which implements AUC-based variables importances as
  described by \bibcitet{party::Janitza+Strobl+Boulesteix:2013}.  Here, the area under the curve
  instead of the accuracy is used to calculate the importance of each variable. 
  This AUC-based variable importance measure is more robust towards class imbalance.

  For right-censored responses, \code{varimp} uses the integrated Brier score as a 
  risk measure for computing variable importances. This feature is extremely slow and
  experimental; use at your own risk.

}
\value{
  A vector of `mean decrease in accuracy' importance scores.
}

\references{\bibshow{*}}

\examples{
    
   set.seed(290875)
   readingSkills.cf <- cforest(score ~ ., data = readingSkills, 
       control = cforest_unbiased(mtry = 2, ntree = 50))

   # standard importance
   varimp(readingSkills.cf)
   # the same modulo random variation
   varimp(readingSkills.cf, pre1.0_0 = TRUE)

   # conditional importance, may take a while...
   varimp(readingSkills.cf, conditional = TRUE)

   \dontrun{
   data("GBSG2", package = "TH.data")
   ### add a random covariate for sanity check
   set.seed(29)
   GBSG2$rand <- runif(nrow(GBSG2))
   object <- cforest(Surv(time, cens) ~ ., data = GBSG2, 
                     control = cforest_unbiased(ntree = 20)) 
   vi <- varimp(object)
   ### compare variable importances and absolute z-statistics
   layout(matrix(1:2))
   barplot(vi)
   barplot(abs(summary(coxph(Surv(time, cens) ~ ., data = GBSG2))$coeff[,"z"]))
   ### looks more or less the same
   }
}
\keyword{tree}
