\name{runInteractivePathAnalysis}
\alias{runInteractivePathAnalysis}
\title{runInteractivePathAnalysis()}
\usage{
runInteractivePathAnalysis(additionalArms = NULL)
}
\arguments{
  \item{additionalArms}{Optional. Function loading
  additional data input arms. Function should take a Study
  object as an argument, use the loadDataArm() to load data
  arms, and return the Study object with arms added.}
}
\value{
The Study object for the created in
runInteractivePathAnalysis
}
\description{
Main function for interactive usage of package. Provides
menu-based access for nearly all program functionality.
}

