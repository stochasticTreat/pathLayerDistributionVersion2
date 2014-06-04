\name{selectAndLoadSettings}
\alias{selectAndLoadSettings}
\title{Load settings from a file folder tree into a Study object.}
\usage{
selectAndLoadSettings(study, fname = NULL, root = "./output/")
}
\arguments{
  \item{study}{A \code{Study} object.}

  \item{fname}{Optional, a \code{string} giving the
  absolute or relative paths of the settings to be loaded.}

  \item{root}{Optional, a \code{string} giving a path to be
  checked for \code{Study} objects from which to retreive
  settings.}
}
\value{
A \code{Study} object with settings loaded.
}
\description{
Allows loading of settigns into a \code{Study} object. If a
directory path name is provided in the fname argument,
settings will be loaded from that file. If no directory is
provided, the user will be prompted to interactively select
a settings folder.
}

