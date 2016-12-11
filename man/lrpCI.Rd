\name{lrpCI}
\alias{lrpCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate confidence intervals for longitudinal regression trees}
\description{
Using bootstrapping to calculate confidence intervals for the first split of a longitudinal regression tree, as outlined in the referenced paper.  The details are outlined in the linked pdf.
}
\usage{
lrpCI(model,B,alpha)
}
\arguments{
  \item{model}{model produced by \code{longRPart}}
  \item{B}{Number of bootstrap samples to run}
  \item{alpha}{The desired confidence level}
}
\details{
To calculate the confidence interval, B bootstrap samples of the patients are taken from the model data set, and a longitudinal regression tree is built, resulting B initial splitting values (one from each tree).  The 100(1-a)\% confidence interval is calculated using the quantiles of the B splitting values.
}
\value{
a vector of split values is returned, from which the confidence intervals can be obtained
}
\references{
Abdolell et al. Binary partitioning for continuous longitudinal data: categorizing a prognostic variable. Statistics in Medicine (2002) vol. 21 pp. 3395-3409 
}
\seealso{
\url{../doc/bootstrap.pdf}
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\keyword{data}
