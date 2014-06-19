\name{compareSources}
\alias{compareSources}
\title{Allows comparrison of summary statistics for all loaded aberration data types.}
\usage{
compareSources(study)
}
\arguments{
  \item{study}{A \code{Study} object. More than one
  aberration data type needs to have been loaded before
  this function is run.}
}
\value{
A \code{Study} object with a slot named, \code{"Aberration
data type comparrison"} in the results list filled with a
table comparring impacts of each aberration data type for
each pathway.
}
\description{
The function examines and compares the pathway summaries
for all loaded aberration data types
}
\examples{
stud = getTestStudyObject()
stud = compareSources(study=stud)
resslot = slot(object=stud, name="results")
View(resslot$'Aberration data type comparrison')
}

