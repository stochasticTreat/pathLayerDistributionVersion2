\docType{class}
\name{DataArm-class}
\alias{DataArm-class}
\title{DataArm}
\description{
Contains the information necessarry to connect a data input
arm to the interactive version of this program.
}
\section{Slots}{

\describe{
\item{\code{description}}{character. The description of the data input arm provided in the main menu.}

\item{\code{title}}{character. The title of the path analysis arm. This is the title by which the results from the analysis arm will be saved in the results list.}

\item{\code{scriptFile}}{A string providing the full or relative file path of a script file containing all functions pertinent to running the data arm. This can be omitted if all functions are already loaded into the namespace, and the mainFunction slot is filled.}

\item{\code{mainFunction}}{function. The name of the main function for the arm in question.}
}}

