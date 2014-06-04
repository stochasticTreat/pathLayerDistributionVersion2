\name{allInteractiveMainFunction}
\alias{allInteractiveMainFunction}
\title{Run pathMod in interactive mode.}
\usage{
allInteractiveMainFunction(additionalArms = NULL)
}
\arguments{
  \item{additionalArms}{Optional. Function loading
  additional data input arms. Function should take a Study
  object as an argument, use the loadDataArm() to load data
  arms, and return the Study object with arms added.}
}
\value{
The Study object created and filled by
allInteractiveMainFunction.
}
\description{
Main function for interactive usage of package. Provides
menu-based access for nearly all program functionality.
}

