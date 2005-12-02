\name{Panel Generating Functions}
\alias{node_inner}
\alias{node_terminal}
\alias{edge_simple}
\alias{node_surv}
\alias{node_barplot}
\alias{node_boxplot}
\alias{node_hist}
\alias{node_density}
\title{ Visualization of Binary Regression Trees }
\description{
    The plot method for \code{BinaryTree} objects is rather flexible and can
    be extended by panel functions. Some pre-defined panel generating 
    functions for the most important cases are documented here.
}
\usage{
node_inner(ctreeobj, digits = 3, abbreviate = FALSE, 
  fill = "white", id = TRUE)
node_terminal(ctreeobj, digits = 3, abbreviate = FALSE, 
  fill = c("lightgray", "white"), id = TRUE)
edge_simple(ctreeobj, digits = 3, abbreviate = FALSE)
node_surv(ctreeobj, ylines = 2, id = TRUE, \dots)
node_barplot(ctreeobj, col = "black", fill = "lightgray", ymax = NULL,
  ylines = 3, widths = 1, gap = NULL, id = TRUE)
node_boxplot(ctreeobj, col = "black", fill = "lightgray",  width = 0.5,
  yscale = NULL, ylines = 3, id = TRUE)
node_hist(ctreeobj, col = "black", fill = "lightgray", freq = FALSE,
  horizontal = TRUE, xscale = NULL, ymax = NULL, ylines = 3,
  id = TRUE, \dots)
node_density(ctreeobj, col = "black", rug = TRUE, horizontal = TRUE,
  xscale = NULL, yscale = NULL, ylines = 3, id = TRUE)
}
\arguments{
  \item{ctreeobj}{ an object of class \code{BinaryTree}.}
  \item{digits}{ integer, used for formating numbers. }
  \item{abbreviate}{ logical indicating whether strings should be 
                     abbreviated. }
  \item{col}{ a color for points and lines. }
  \item{fill}{ a color to filling rectangles. }
  \item{id}{ a logical indicating whether the node ID should be plotted.}
  \item{ylines}{ number of lines in y-direction. }
  \item{widths}{ widths in barplots. }
  \item{width}{ width in boxplots. }
  \item{gap}{ gap between bars in a barplot (\code{node_barplot}). }
  \item{yscale}{ limits in y-direction}
  \item{xscale}{ limits in x-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{horizontal}{ logical indicating if the plots should be horizontal. }
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a representation
              of frequencies. If \code{FALSE}, probabilities are plotted.}
  \item{rug}{logical indicating if a rug representation should be added. }
  \item{\dots}{ additional arguments passed to callies.}
}
\details{

  The \code{plot} method for \code{BinaryTree} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. For convenience, the panel
  generating functions \code{node_inner} and \code{edge_simple} 
  return panel functions to draw inner nodes and left and right edges. 
  For drawing terminal nodes, the functions returned by the other panel 
  functions can be used. The panel generating function \code{panel_terminal_node} 
  is a terse text-based representation of terminal nodes.

  Graphical representations of terminal nodes are available and depend on
  the measurement scale of the response variable the tree was fitted to.
  \code{node_surv} returns a functions that plots Kaplan-Meier curves in each 
  terminal node. \code{node_barplot}, \code{node_boxplot}, \code{node_hist} and
  \code{node_density} can be used to plot barplots, boxplots, histograms and 
  estimated densities into the terminal nodes.

}
\examples{

  data("airquality")
  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ### very simple; the mean is given in each terminal node
  plot(airct, type = "simple")

  ### density estimates
  plot(airct, terminal_panel = node_density)
    
  ### histograms 
  plot(airct, terminal_panel = node_hist)

  ### boxplots
  plot(airct, terminal_panel = node_boxplot, drop_terminal = TRUE)

}
\keyword{hplot}
