\name{lrpCV}
\alias{lrpCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs cross-validation on the longRPart model}
\description{
This functions implements v-fold cross validation on the longRPart trees.  IT IS CURRENTLY UNDER DEVELOPMENT AND SHOULD NOT BE USED FOR ANALYSIS.
}
\usage{
lrpCV(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{model produced by \code{longRPart}}
}
\details{
Using the xval in the the rpart control object, this function performs v-fold cross-validation.  The object returned is deviance between the subject and the predicted values at the node it was navigated to as a test case.  The mean of that vector would be the average deviances.
}
\references{
Abdolell et al. Binary partitioning for continuous longitudinal data: categorizing a prognostic variable. Statistics in Medicine (2002) vol. 21 pp. 3395-3409 
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\keyword{data}
