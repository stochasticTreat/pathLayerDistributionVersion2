\name{loadBasicArms}
\alias{loadBasicArms}
\title{Load the basic set of analysis arms for use in a study}
\usage{
loadBasicArms(STUDY)
}
\arguments{
  \item{STUDY}{A \code{Study} object.}
}
\value{
A \code{Study} object with study arms loaded.
}
\description{
This function connects the basic set of analysis arms (drug
screen analysis, somatic mutation, arbitrary gene input
data and overlap analysis) to a Study object.
}
\examples{
study = getStudyObject(study.name="testDataSets",
											 geneIdentifierType="HUGO")
study = loadBasicArms(STUDY=study)
sa = slot(object=study, name="arms")
print(names(sa))
print(slotNames(sa$functional_drug_screen_summary))
}

