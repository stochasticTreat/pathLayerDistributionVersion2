\name{DataSummary}
\alias{DataSummary}
\title{Display summaries of loaded data.}
\usage{
DataSummary(study)
}
\arguments{
  \item{study}{A \code{Study} object.}
}
\description{
This function prints to the screen summary statistics for
data sets from each data arm.
}
\examples{
study=getTestStudyObject()
DataSummary(study=study)
}

