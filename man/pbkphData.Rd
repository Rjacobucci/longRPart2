\name{pbkphData}
\alias{pbkphData}
\title{A dataset on the improvement in hearing over time for children with coclear implants}
\description{
This is the original dataset from Mohamed Abdolell's paper.
}
\details{
Note the format of this dataset, as this is the required format for the data to work.  The data should be in LONG format, i.e., it should have one row per observation, which means that a patient should have as many rows as they have time points.  For each row there needs to be a patient identifier, a time value and a response value.  Also note that, for the variables that are being used as splits, they must be completely populated, and the same for each patient.  For this particular dataset that means that age and gender are completely populated, having missing values in them would cause those observations to be dropped.
}
\references{
Abdolell et al. Binary partitioning for continuous longitudinal data: categorizing a prognostic variable. Statistics in Medicine (2002) vol. 21 pp. 3395-3409 
}
\keyword{data}
