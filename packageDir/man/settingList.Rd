\name{settingList}
\alias{settingList}
\title{settingList}
\usage{
settingList(s, prompt, set)
}
\arguments{
  \item{s}{the settings list}

  \item{prompt}{The message to be prompted to the user.}

  \item{set}{data frame matrix or vector: the options the
  user has to select amongst}
}
\value{
The settings list object is returned with the list of user
selections in the prompt slot and the full list of options
in a slot named prompt, '<prompt>(Full list of options.)'
}
\description{
Allows a user to select from a list of options.
}

