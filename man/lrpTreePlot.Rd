\name{lrpTreePlot}
\alias{lrpTreePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standard rpart tree plot adapted to longitudinal data}
\description{
Using the standard \pkg{rpart} plotting function and subfunction \code{rpartco}, this function produces a tree plot with simple line plots at the nodes of the regression coefficients of a mixed effects model.
}
\usage{
lrpTreePlot(model,use.n=TRUE,colors=NULL,place='bottomright')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Model produced by \code{longRPart}}
  \item{use.n}{Boolean indicating whether to include the number of observations at a leaf node in the plot.  If \code{F}, then the slope of the line is printed instead}
  \item{colors}{The desired colors for each node.}
  \item{place}{THe placement of the legend}
}
\details{
This function uses the standard \code{plot} and \code{text} commands from \pkg{rpart} and \code{longRPart}, but includes lineplots of the mixed-effects models calculated at each node.
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\seealso{
\code{\link{longRPart}}
\code{\link{lrpPlot}}
}
\keyword{data}
