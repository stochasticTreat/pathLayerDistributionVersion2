\name{autoRunFromSettings}
\alias{autoRunFromSettings}
\title{Automatically run set of analyses from a set of loaded settings.}
\usage{
autoRunFromSettings(study, verbose = T)
}
\arguments{
  \item{study}{A \code{Study} object with settings
  available.}

  \item{verbose}{A flag indicating if analyses should be
  run in interactive/verbose mode, prompting user to check
  intermediate settings.}
}
\value{
A \code{Study} object with analysis described by settings
complete.
}
\description{
Function searches for settings objects in the
STUDY-studyMetaData-settings slot, then automatically
executes the analysis for each arm described by the
settings.
}

