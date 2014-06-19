\name{abDrugOverlapAnalysis}
\alias{abDrugOverlapAnalysis}
\title{abDrugOverlapAnalysis}
\usage{
abDrugOverlapAnalysis(study, thresh = 0.05, enrich_col = "hyperg_p_w_FDR",
  verbose = T, settings = NULL)
}
\arguments{
  \item{study}{A study object (must have results loaded for
  aberration and functional analyses)}

  \item{thresh}{Numeric, the threshold used to establish
  which pathways will be considered to be drug sensitive
  and aberrational.}

  \item{enrich_col}{String, the name of the column in the
  path summary used to determine significance of a pathway.
  This is the column to which the threshold argument is to
  be applied.}

  \item{verbose}{Logical, a flag indicating if the program
  is to provide extra diagnostic output.}

  \item{settings}{List obejct. This contains the settings
  used in the overlap analysis. Use
  getDefaultOverlapSettings()}
}
\value{
Study object with overlap analysis in results list added or
replaced
}
\description{
main function for running analysis of functional and
genomic aberration arms
}
\examples{
#get a study object pre-loaded with genomic aberration and functional genomic data.
stud = getTestStudyObject()

#get default overlap analysis settings
defSettings = getTestStudySettings()
overlapSettings = defSettings$overlap_analysis

#run the overlap analysis

overlapRes1 = abDrugOverlapAnalysis(study=stud, settings=overlapSettings)

\dontrun{
#If settings for the overlap analysis cannot be found, the overlap
#analysis will run in interactive mode, requiring user input.
#This line of code will run the overlap anaysis in interactive mode,
#requesting user input to establish settings:

abDrugOverlapAnalysis(study=stud)
}
}

