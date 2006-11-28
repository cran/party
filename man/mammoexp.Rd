\name{mammoexp}
\alias{mammoexp}
\docType{data}
\title{ Mammography Experience Study }
\description{
  Data from a questionaire on the benefits of mammography.
}
\usage{data(mammoexp)}
\format{
  A data frame with 412 observations on the following 6 variables.
  \describe{
    \item{ME}{Mammograph experience, an ordered factor with 
              levels \code{Never} < \code{Within a Year} < \code{Over a Year}}
    \item{SYMPT}{Agreement with the statement: 
                 `You do not need a mamogram unless you develop symptoms.' 
                 A factor with levels \code{Strongly Agree}, \code{Agree},
                 \code{Disagree} and \code{Strongly Disagree}}
    \item{PB}{Perceived benefit of mammography, the sum of five scaled
              responses, each on a four point scale.
              A low value is indicative of a woman with strong agreement 
              with the benefits of mammography.}
    \item{HIST}{ Mother or Sister with a history 
                 of breast cancer; a factor with levels \code{No} and
                 \code{Yes}.}
    \item{BSE}{Answers to the question: 
               `Has anyone taught you how to examine your own breasts?'
               A factor with levels \code{No} and \code{Yes}. }
    \item{DECT}{Answers to the question: 
                `How likely is it that a mammogram
                could find a new case of breast cancer?'
                An ordered factor with levels \code{Not likely} < 
                \code{Somewhat likely} < \code{Very likely}.}
  }
}
\source{

  Hosmer and Lemeshow (2000). \emph{Applied Logistic Regression}, 2nd
  edition. John Wiley & Sons Inc., New York. Section 8.1.2, page 264.

}
\examples{

  ### fit tree with attached scores (equal to the default values)
  mtree <- ctree(ME ~ .,data = mammoexp, 
        scores = list(ME = 1:3, SYMPT = 1:4, DECT = 1:3))
  mtree
  plot(mtree)
}
\keyword{datasets}
