\name{RandomForest-class}
\docType{class}
\alias{RandomForest-class}
\alias{treeresponse,RandomForest-method}
\alias{weights,RandomForest-method}

\title{Class "RandomForest"}
\description{A class for representing random forest ensembles. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RandomForest", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ensemble}:}{Object of class \code{"list"}, each element
                           being an object of class \code{\link{BinaryTree-class}}.}
    \item{\code{data}:}{ an object of class \code{\link[modeltools]{ModelEnv}}.}
    \item{\code{responses}:}{ an object of class \code{"VariableFrame"}
                              storing the values of the response variable(s). }
    \item{\code{cond_distr_response}:}{ a function computing the conditional
                                        distribution of the response. }
    \item{\code{predict_response}:}{ a function for computing predictions. }
    \item{\code{prediction_weights}:}{ a function for extracting weights from
                                     terminal nodes. }
  }
}
\section{Methods}{
  \describe{
    \item{treeresponse}{\code{signature(object = "RandomForest")}: ... }
    \item{weights}{\code{signature(object = "RandomForest")}: ... }
	 }
}
\examples{

    ### honest (i.e., out-of-bag) cross-classification of 
    ### true vs. predicted classes
    table(mammoexp$ME, predict(cforest(ME ~ ., data = mammoexp, 
                               control = cforest_control(ntree = 50)), 
                               OOB = TRUE))
}
\keyword{classes}
