\name{Plot BinaryTree}
\alias{plot.BinaryTree}
\title{ Visualization of Binary Regression Trees }
\description{
  This method plots objects of class \code{BinaryTree}.
}
\usage{
\method{plot}{BinaryTree}(x, main = NULL,
     type = c("extended", "simple"), terminal_panel = NULL,
     inner_panel = node_inner, edge_panel = edge_simple,
     drop_terminal = (type[1] == "extended"), 
     tnex = (type[1] == "extended") + 1, new = TRUE,  
     pop = FALSE, \dots)
}
\arguments{
  \item{x}{ an object of class \code{BinaryTree}.}
  \item{main}{ an optional title for the plot.}
  \item{type}{ a character specifying the complexity of the plot:
               \code{extended} tries to visualize the distribution of the
               response variable in each terminal node whereas \code{simple} 
               only gives some summary information.}
  \item{terminal_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the terminal nodes.
                         Alternatively, a panel generating function of the form
                         \code{function(ctreeobj, \dots)} can be supplied.}
  \item{inner_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the inner nodes.}
  \item{edge_panel}{ an optional panel function of the form 
                    \code{function(split, ordered = FALSE, left = TRUE)}
                    plotting the edges.}
  \item{drop_terminal}{ a logical indicating whether all terminal nodes
                        should be plotted at the bottom.}
  \item{tnex}{ a numeric value giving the terminal node extention. }
  \item{new}{ a logical indicating whether a new plot should be created. }
  \item{pop}{ a logical whether the viewport tree should be popped before
              return. }
  \item{\dots}{ additional arguments passed to callies.}
}
\details{

  This \code{plot} method for \code{BinaryTree} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. Panel functions for plotting
  inner nodes, edges and terminal nodes are available for the most important
  cases and can serve as the basis for user-supplied extensions, see
  \code{\link{node_inner}}.

}
\examples{

  data(airquality)
  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ### boxplots in each node
  plot(airct, terminal_panel = node_boxplot, drop_terminal = TRUE)

}
\keyword{hplot}
