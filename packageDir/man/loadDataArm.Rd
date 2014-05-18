\name{loadDataArm}
\alias{loadDataArm}
\title{loadDataArm}
\usage{
loadDataArm(mainFunction, arms, title, description,
  scriptFile = "no script file provided, loading data arm main function from local environment")
}
\arguments{
  \item{mainFunction}{The main execution function for the
  arm. If the script provided, the function does not need
  to be loaded into memory when the call is made to
  loadDataArm -- the script file will be run then the
  function will be called upon to be placed into the
  mainFunction slot of a DataArm object.}

  \item{arms}{The list containing all DataArm objects for
  the study.}

  \item{title}{character string. The one word description
  of the data arm. (ex: functional_drug_screen_summary)
  This will be used internally to refer to any results
  associated with the data arm, and externally to name
  folders when saving data.}

  \item{description}{The description of the data input arm
  to be used as a main menu option. (ex: "Process drug
  screen data")}

  \item{scriptFile}{A file providing the script where the
  functions implementing the data import arm can be found.
  Note: This allows that a script can be provied which
  contains the definition of the function provided for the
  mainFunction argument. Thus, this script will be executed
  before the mainFunction is called or referened.}
}
\value{
The arms list with the additional arm added.
}
\description{
Connects data input arm to program, providing essential
information and the main execution function.
}

