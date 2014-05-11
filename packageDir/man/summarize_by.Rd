\name{summarize_by}
\alias{summarize_by}
\title{summarize_by}
\usage{
summarize_by(col, display = F, barPlotTitle = "Counts across types",
  left_margin_factor = 1)
}
\arguments{
  \item{col}{Vector, the categorical values to be
  summarized}

  \item{display}{A flag indicating if the data should be
  displayed in a bar plot.}

  \item{barPlotTitle}{The title that should be given to the
  bar plot}

  \item{left_margin_factor}{A numerical factor expanding
  the left margin of the bar plot. Used to adjust in
  situations where there is insufficient room to display
  the text of each bar's label}
}
\value{
A two column data.frame with columns "types", the different
levels and "counts", the different numer of times each
level occurs
}
\description{
summarizes the data levels in a column of categorical data
and optionally provides a barplot output
}

