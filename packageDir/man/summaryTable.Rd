\name{summaryTable}
\alias{summaryTable}
\title{summaryTable}
\usage{
summaryTable(study, activeGeneDescription = "active_gene",
  dataSetDescription = "analysis of undefined active genes", pgm = NULL,
  originalDataMatrix = NULL,
  settings = getDefaultSettings()$defaultSummaryTable, coverage = NULL,
  coverageGeneDescription = "covered_gene",
  coverageDataSetDescription = "Coverage analysis",
  individualEnrichment = NULL)
}
\arguments{
  \item{study}{A Study object.}

  \item{activeGeneDescription}{One word description of how
  the 'active' genes should be described (ex: mutated,
  amplified, drug_sensitive)}

  \item{dataSetDescription}{Description of the data set
  provided for the pathway analysis. (ex: "Somatically
  mutated genes")}

  \item{pgm}{The patient gene matrix for the patient data.
  Bipartate graph with row names set as gene names, column
  names set as patient IDs and values logical indicating if
  gene is considered "active" in patient.}

  \item{originalDataMatrix}{The original data matrix. Makes
  this available for novel types of pathway analysis, not
  treating genes as discretely "on" or "off"}

  \item{settings}{A settings list object controling how the
  path analysis will be conducted. Default given by call to
  getDefaultSettings()$defaultSummaryTable .}

  \item{coverage}{character vector of gene names. Should be
  provided if coverage provided by the analysis platform
  used to proved gene data is not considered to have full
  genome coverage.}

  \item{coverageGeneDescription}{The description of genes
  that are covered be the analysis platform (ex:
  drug_targeted or sequenced)}

  \item{coverageDataSetDescription}{A longer description of
  the set of genes provided in the coverage analysis}

  \item{individualEnrichment}{logical. A flag indicating if
  pathway analysis should be run for individual patients.}
}
\description{
The main function for conducting and pathway analysis for a
set of gene data.
}

