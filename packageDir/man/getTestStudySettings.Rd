\name{getTestStudySettings}
\alias{getTestStudySettings}
\title{Get settings for the three basic study arms: functional, aberration and function-aberration overlap.}
\usage{
getTestStudySettings()
}
\value{
\code{list} object with slots named "overlap_analysis",
"somatic_mutation_aberration_summary" and
"functional_drug_screen_summary" containing settings lists
for their respective study arms.
}
\description{
Returns the settings list for the analysis of the test data
sets: functional drug screen data, somatic mutation data
and overlap analysis, of the functional and somatic
mutation (aberration) data.
}

