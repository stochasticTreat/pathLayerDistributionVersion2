\name{addSequenceCaptureArm}
\alias{addSequenceCaptureArm}
\title{Add sequence capture data input arm.}
\usage{
addSequenceCaptureArm(STUDY)
}
\arguments{
  \item{STUDY}{A \code{Study} object.}
}
\value{
STUDY A \code{Study} object with the sequence capture data
input arm loaded.
}
\description{
Add arm onto main study for input of sequence capture data.
}
\examples{
STUDY = getTestStudyObject()
STUDY = addSequenceCaptureArm(STUDY=STUDY)
\dontrun{
STUDY = allInteractiveMainFunction(additionalArms=addSequenceCaptureArm)
}
}

