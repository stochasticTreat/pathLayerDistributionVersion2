\name{addCohortToIndividuals}
\alias{addCohortToIndividuals}
\title{Adds path enrichment information from cohort to paths found to contain active genes in individual patients.}
\usage{
addCohortToIndividuals(res, columnsToAdd = grep(pattern = "hyper",
  colnames(res$pathsummary))[1]:ncol(res$pathsummary))
}
\arguments{
  \item{res}{The results list object containing slots named
  path_summary_each_patient and pathsummary.}

  \item{columnsToAdd}{The columns from mainPathSum to add
  to each individual patient path summary}
}
\value{
The list of path summary results for each patient.
}
\description{
Adds the columns in columnsToAdd from the whole cohort path
enrichment to the path enrichments for each patient.
}
\examples{
tstudy = getTestStudyObject()
all_results = slot(object=tstudy, name="results")
psep = addCohortToIndividuals(res=all_results$somatic_mutation_aberration_summary)
View(psep$p1$pathsummary)
}

