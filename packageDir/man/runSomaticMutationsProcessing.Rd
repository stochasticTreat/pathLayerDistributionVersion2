\name{runSomaticMutationsProcessing}
\alias{runSomaticMutationsProcessing}
\title{Main execution function for entry of somatic mutation data, and somatic mutations input arm of this program.}
\usage{
runSomaticMutationsProcessing(settings, study)
}
\arguments{
  \item{settings}{The settings list object for the somatic
  mutation interface, or a blank list.}

  \item{study}{A Study object (a defined in
  initiateDataStructures.R)}
}
\value{
Somatic mutations path results set
}
\description{
Implements armMain interface, providing main execution
function for entry of somatic mutation data, and somatic
mutations input arm of this program.
}
\note{
This function is passed as armMain() inside the runArm()
function
}

