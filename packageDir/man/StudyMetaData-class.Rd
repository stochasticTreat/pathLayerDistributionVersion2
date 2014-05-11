\docType{class}
\name{StudyMetaData-class}
\alias{StudyMetaData-class}
\title{StudyMetaData}
\description{
Object containing crucial study meta data.
}
\section{Slots}{

\describe{
\item{\code{paths}}{Path_Detail reference class object}

\item{\code{settings}}{A list of settings list objects. One slot in list is generally provided for each arm.}

\item{\code{studyName}}{String. The name of the study and name of the file folder used to save the study.}

\item{\code{geneIdentifierType}}{String. The name of the type of gene identifiers used in the study.}

\item{\code{GeneIdentifierLookup}}{depricated. A slot for a set of official, approved gene identifiers and mappings for unofficial gene identifiers.}

\item{\code{RootFile}}{character data. The path to folder where the study is saved.}
}}

