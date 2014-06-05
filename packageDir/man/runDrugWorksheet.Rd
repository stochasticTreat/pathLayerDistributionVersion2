\name{runDrugWorksheet}
\alias{runDrugWorksheet}
\title{Run the drug selection worksheet.}
\usage{
runDrugWorksheet(STUDY = STUDY)
}
\arguments{
  \item{STUDY}{A \code{Study} object into which
  experimental data has been loaded.}
}
\value{
The \code{data.frame} displayed in the drug selection
worksheet, relating drug targets, drugs, clinical trial
phase, and other data.
}
\description{
Running this function will open up an HTML page allowing
interactive selection and exploration of targeted drugs.
Information is supplied to link genes to drugs and filter
by clinical trial phase, gene mutation in a cohort and
cellular pathways.
}

