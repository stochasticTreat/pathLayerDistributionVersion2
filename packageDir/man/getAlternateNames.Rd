\name{getAlternateNames}
\alias{getAlternateNames}
\title{getAlternateNames}
\usage{
getAlternateNames(nodeIds, df)
}
\arguments{
  \item{nodeIds}{The ids for the nodes whose alternate
  names are being searched for.}

  \item{df}{The data frame containing the biopax data, as
  given by the rBiopaxParser package function readBiopax}
}
\value{
data.frame with columns nodeID and altName.
}
\description{
looks through the data frame assocaited with biopax file
for alternate names associated with a node/entity
}

