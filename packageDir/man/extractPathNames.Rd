\name{extractPathNames}
\alias{extractPathNames}
\title{Get names of significant pathways.}
\usage{
extractPathNames(resSetName, study, limitCol, limitVal)
}
\arguments{
  \item{resSetName}{the name of a summary table results
  set}

  \item{study}{a study object}

  \item{limitCol}{the column the returned pathways should
  be limited by}

  \item{limitVal}{the value in the limitCol the returned
  pathways should be limited to}
}
\value{
vector, a list of pathways from study meeting the
parameters described in the function argument
}
\description{
Extracts path names given a results set name, study and
limitCol
}

