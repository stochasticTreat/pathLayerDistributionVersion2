\name{getStudyObject}
\alias{getStudyObject}
\title{getStudyObject}
\usage{
getStudyObject(study.name = "", settings = list(), resf = list(),
  GeneIdentifierLookup = data.frame(), path_detail = Path_Detail$new(),
  arms = list(), geneIdentifierType = "HUGO", rootFolder = "./output/")
}
\arguments{
  \item{study.name}{The name of the study (this will be
  used as a folder name if study is saved, so this should
  be valid as a file name.)}

  \item{settings}{A settings list object.}

  \item{resf}{a results list object.}

  \item{GeneIdentifierLookup}{Depricated. A table of
  approved gene identifiers.}

  \item{path_detail}{A Path_Detail object.}

  \item{arms}{The list of arm objects.}

  \item{geneIdentifierType}{String, a one word name of the
  type of gene identifiers used (ex : HUGO, Uniprot)}

  \item{rootFolder}{character string. The file path where
  the study should be saved.}
}
\description{
Function to initilize a new study object.
}
\examples{
study = getStudyObject() #this will give a study object with no pathways loaded (pathways must first be loaded for analyses to take place)
study2 = getStudyObject(study.name="testStudy", path_detail=getDefaultPaths())
}

