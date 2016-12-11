\name{longRPart}
\alias{longRPart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursive partitioning of longitudinal data}
\description{
This function uses the custom partitioning rules in \pkg{rpart} with the mixed effects models from \pkg{nlme} to construct regression trees out of longitudinal data.
}
\usage{
longRPart(lmeFormula,rPartFormula,randomFormula,data,weight=NULL,R=NULL,control = rpart.control())
}
\arguments{
  \item{lmeFormula}{The mixed-effects formula, as would be passed to the \code{lme} function}
  \item{rPartFormula}{a list of variables to be used as potential classifiers for the rpart function}
  \item{randomFormula}{The grouping equation, as would be passed to the \code{lme} function}
  \item{data}{The data set being analyzed}
  \item{weight}{A weighting of the observations to be passed to the rpart function.}
  \item{R}{The correlation matrix for the mixed effects model.  Default is exponential, which is a generalization of the corAR1 structure}
  \item{control}{Standard rpart control structure}
}
\details{
Applying CART principles to longitudinal data using mixed effects models is outlined in the paper by Abdolell et al. refereneced below.  The basic principle is the same as the more common partitioning processes, with the deviance function for longitudinal data being taken from the mixed effects model created at each node.

The arguments passed to the function can be classified as either rpart arguments or lme arguments: rPartFormula and control are used to build the rpart structure, while lmeFormula, randomFormula and R are used to build the lme models.

To pass data to the function, it must have one observation per row, with at least the following three columns: a response column, an id column (to connect observations on the same subject) and a time column to record the time of response.  Other columns will be required for the covariates.

The function passes many of the arguments directly to the \code{rpart} or \code{lme} functions, so before the tree is a built, a regular mixed-effects model should be built using the \code{lme} function.

\code{rPartFormula} is a list of covariates that will be used as potential splits in \code{rpart}.  The basic structure of the function is \code{~cov1+cov2+cov3+...}.  Note that \code{~.} will not work in this situation, because rpart will try to include all the variables (including ID and tme) and the resulting classification will not make sense.
}
\value{
The object returned is a standard rpart object, with some additional attributes to be used in the plotting functions
}
\references{
Abdolell et al. Binary partitioning for continuous longitudinal data: categorizing a prognostic variable. Statistics in Medicine (2002) vol. 21 pp. 3395-3409 
}
\author{Mohamed Abdolell <mo@dal.ca> and Sam Stewart}
\seealso{
\code{\link{rpart}}
\code{\link{lme}}
\code{\link{corClasses}}
\code{\link{lrpPlot}}
\code{\link{lrpTreePlot}}
\code{\link{lrpPVal}}
\code{\link{lrpCI}}
}
\examples{
# data from Abdolell paper
data(pbkphData)
pbkphData$Time=as.factor(pbkphData$Time)

model = longRPart(pbkph~Time,~age+gender,~1|Subject,pbkphData,R=corExp(form=~time),control=rpart.control(minbucket=80,xval=10))

}
\keyword{data}
