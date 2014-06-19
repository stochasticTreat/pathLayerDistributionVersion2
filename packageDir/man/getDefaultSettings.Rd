\name{getDefaultSettings}
\alias{getDefaultSettings}
\title{getDefaultSettings}
\usage{
getDefaultSettings()
}
\value{
A list with a slot named "defaultSummaryTable", containing
the default settings for the summaryTable function.
}
\description{
Returns the default settings object for the summaryTable
path analysis.
}
\examples{
defSettings = getDefaultSettings()
print(names(defSettings$defaultSummaryTable))
}

