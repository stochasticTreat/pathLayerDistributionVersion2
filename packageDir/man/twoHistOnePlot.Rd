\name{twoHistOnePlot}
\alias{twoHistOnePlot}
\title{Puts two histograms on the same plot.}
\usage{
twoHistOnePlot(dataset1, dataset2, frequency_not_density = T,
  main_title = "main_title not set", x_label = "x_label not set",
  y_label = "Frequency", legend_titles = c("dataset 1", "dataset 2"),
  breaks = c(30, 30))
}
\arguments{
  \item{dataset1}{A numeric data set; \code{data.frame or
  matrix} with one column.}

  \item{dataset2}{A numeric data set; \code{data.frame or
  matrix} with one column.}

  \item{frequency_not_density}{Depricated.}

  \item{main_title}{The main title for the plot.}

  \item{x_label}{The label to be put along the x axis.}

  \item{y_label}{The label to be put along the y axis.}

  \item{legend_titles}{Character vector, length 2, giving
  the names of the data in dataset1 and dataset2.}

  \item{breaks}{Depricated.}
}
\description{
Puts two histograms on the same plot, shows each with a
different color, and overlap with a third color.
}

