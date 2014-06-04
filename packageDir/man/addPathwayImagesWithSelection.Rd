\name{addPathwayImagesWithSelection}
\alias{addPathwayImagesWithSelection}
\title{Add network diagrams of affected pathways.}
\usage{
addPathwayImagesWithSelection(study, limitCol = "hyperg_p_w_FDR",
  limitVal = 0.05)
}
\arguments{
  \item{study}{A \code{Study} object.}

  \item{limitCol}{The column from the results that is used
  to limit the number of path diagrams displayed.}

  \item{limitVal}{The value in the limitCol used to limit
  the number of pathways. Net work diagrams will only be
  produced for pathways with values in limitCol which are
  smaller than limitVal.}
}
\value{
A \code{Study} object will be returned with file names of
network diagrams added to the imageSlots slot of the
appropriate results sets.
}
\description{
Outputs to cytoscape diagrams of affected pathways .
}

