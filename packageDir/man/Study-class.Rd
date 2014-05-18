\docType{class}
\name{Study-class}
\alias{Study-class}
\title{Study}
\description{
Contains all data related to a single study, including
results for each data input arm and the overlap analysis,
the main functions for each data input arm, and the study
meta data object.
}
\section{Slots}{

\describe{
\item{\code{results}}{list. Names provided as data input arm names. Each slot contains the results set from its respective data input arm.}

\item{\code{arms}}{A list with each slot containing a DataArm object, slot names given as the names of each data input arm}

\item{\code{studyMetaData}}{the StudyMetaData object}
}}

