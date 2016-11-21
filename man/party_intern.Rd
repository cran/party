\name{party_intern}
\alias{party_intern}
\title{
Call internal functions.
}
\description{
Call one of the internal party functions.
}
\usage{
party_intern(..., fun = c("R_TreeGrow", "R_get_nodeID", 
             "R_getpredictions", "initVariableFrame", 
             "ctreedpp", "newinputs"))
}
\arguments{
  \item{\dots}{Arguments to \code{fun}.}
  \item{fun}{The name on an internal party function.}
}
\details{
This function must not be called under any circumstances.
}
\keyword{internal}
