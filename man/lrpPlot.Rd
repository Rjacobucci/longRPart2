\name{lrpPlot}
\alias{lrpPlot}

\title{Line plot of the nodes of a longRPart object}
\description{
Plotting the results of a longRPart regression tree: each line in the plot represents a different node
}
\usage{
lrpPlot(model,smoothing="n",color=NULL,place="bottomright")
}
\arguments{
  \item{model}{Output from \code{longRPart} function}
  \item{smoothing}{Type of smoothing to apply to the lines (NOT IMPLEMENTED)}
  \item{color}{colors to label the lines, default is rainbow()}
  \item{place}{placement for the legend}
}
\details{
This function performs a simple line plot of the regression coefficients from the mixed-effects models of each of the nodes in the longitudinal regression tree.  Eventually smoothing will be applied, but for now only 'n', for none, is implemented.
}
\seealso{
\code{\link{longRPart}}
\code{\link{lrpTreePlot}}
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\keyword{data}
