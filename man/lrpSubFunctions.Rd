\name{lrpSubFunctions}
\alias{lrpSubFunctions}
\alias{rpartco}
\alias{tree.depth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Helper functions from the rpart library}
\description{
Helper functions used in the plotting process for \code{longRPartTreePlot}
}
\usage{
rpartco(tree, parms =  paste(".rpart.parms", dev.cur(), sep = "."))
tree.depth(nodes)
}
\arguments{
  \item{tree}{Model returned from rpart}
  \item{parms}{secondary data for the rpartco function}
  \item{nodes}{The nodes of a tree}
}
\details{
\code{rpartco} is used to find the (x,y) coordinates of a tree-plot.

\code{tree.depth} is used to find the depth of a set of nodes.
}
\value{
each return approriate vectors/lists
}
\seealso{
\code{\link{lrpTreePlot}}
\code{\link{plot.rpart}}
}
\keyword{data}
