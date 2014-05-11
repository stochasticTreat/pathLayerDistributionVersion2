\name{initiateStudy}
\alias{initiateStudy}
\title{initiateStudy}
\usage{
initiateStudy(studyFolderName = NULL, path_detail = NULL,
  root = "./output", study.name = NULL)
}
\arguments{
  \item{studyFolderName}{The name of the study folder.}

  \item{path_detail}{A Path_Detail object.}

  \item{root}{The path to the folder containing a study to
  be loaded}

  \item{study.name}{The name of the study to be initiated}
}
\value{
The initiated Study object.
}
\description{
Orchestrates establishment of a new study, or loading of an
old study
}

