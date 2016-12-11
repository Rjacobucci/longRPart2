\name{lrpPVal}
\alias{lrpPVal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate p-values for the longitudinal regression tree}
\description{
Using permutation tests, this function calculates the p-value for the first split of the tree; see the linked pdf for details.  It also calculates the p-value for the entire tree, though the statistical theory behind this p-value is not as clear.
}
\usage{
lrpPVal(model,J)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{model produced by \code{longRPart}}
  \item{J}{number of permutations to run}
}
\details{
The basic principle behind the test is to randomized the covariates to the subject and then recalculate the tree: the p-value is the proportion of trees that have a greater reduction in deviance than the original tree.
}
\note{J indicates that a longRPart tree will be built J times, and a reasonable value for J would be 100 (at least), which should take 100 times as long as the original longRPart function}
\references{
Abdolell et al. Binary partitioning for continuous longitudinal data: categorizing a prognostic variable. Statistics in Medicine (2002) vol. 21 pp. 3395-3409 
}
\seealso{
\url{../doc/permutation.pdf}
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\keyword{data}
