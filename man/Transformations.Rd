\name{Transformations}
\alias{ptrafo}
\title{ Function for Data Transformations }
\description{
  Transformations of Response or Input Variables
}
\usage{
ptrafo(data, numeric_trafo = id_trafo, 
    factor_trafo = function(x) model.matrix(~ x - 1), 
    surv_trafo = logrank_trafo, 
    var_trafo = NULL)
}
\arguments{
  \item{data}{an object of class \code{data.frame}.}
  \item{numeric_trafo}{a function to by applied to \code{numeric} 
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{factor_trafo}{a function to by applied to \code{factor}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually a dummy or contrast 
    matrix).}
  \item{surv_trafo}{a function to by applied to 
    elements of class \code{Surv} of \code{data} returning a 
    matrix with \code{nrow(data)} rows and an arbitrary number of columns.}
  \item{var_trafo}{an optional named list of functions to be applied to the
                   corresponding variables in \code{data}.}
}
\details{

  \code{trafo} applies its arguments to the elements of \code{data}
  according to the classes of the elements. See \code{\link[coin]{trafo}}
  for more documentation and examples.

}
\value{

  A named matrix with \code{nrow(data)} rows and 
  arbitrary number of columns.

}
\examples{

### rank a variable
ptrafo(data.frame(y = 1:20), numeric_trafo = rank)

}
\keyword{manip}
