\name{runArm}
\alias{runArm}
\title{runArm}
\usage{
runArm(armDescription, study, fromDescription = T, interactive = NULL)
}
\arguments{
  \item{armDescription}{the desription slot from a DataArm
  object}

  \item{study}{A Study object.}

  \item{fromDescription}{Indicates if arm should be
  retrieved by its description or by its title.}

  \item{interactive}{Indicates if arm should be run in
  interactive mode.}
}
\value{
Returns a Study object with results from running the study
arm.
}
\description{
Main function for taking and running a study arm. Provides
needed settings and study object, implementing the study
arm interface.
}

