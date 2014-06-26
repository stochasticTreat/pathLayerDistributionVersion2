\docType{class}
\name{PathSummaryRunner-class}
\alias{PathSummaryRunner-class}
\title{PathSummaryRunner}
\description{
This reference class object is used in the pathway
analysis. It provides all available data for internal
functions and user-defined functions to conduct pathway
analysis.
}
\section{Fields}{
  \describe{ \item{\code{patientGeneMatrix}:}{ \code{Object
  of class "Matrix"} Bipartate graph with patient ids as
  columns, genes as rows and values, logical indicating if
  gene is 'active' (drug sensitive, aberrational, etc.) in
  patient.  } \item{\code{coverage_summary}:}{ \code{Object
  of class "list"} Same data as is found in [[result]] set,
  only filled if analysis platform in question is
  coverage-limited, ie, coverage slot of PathSummaryRunner
  is filled.  } \item{\code{coverage}:}{ \code{Object of
  class "Character"} A character vector containing symbols
  for the set of species (ie genes) covered by the analysis
  platform. If this is provided, pathway analysis will be
  limited to included only genes in the coverage set. This
  should be provided if an analysis platform is not
  considered to have full-genome coverage (ex: sequence
  capture data for 3000 unique genes.).  }
  \item{\code{indPatientSummarySettings}:}{ \code{Object of
  class "list"} Settings object for analysis of individual
  patients (as opposed to a full cohort).  }
  \item{\code{original_data_matrix}:}{ \code{Object of
  class "matrix"} The original, untransformed data set.  }
  \item{\code{min_gene_frequency}:}{ \code{Object of class
  "numeric"} The minimum frequency across the cohort for a
  gene to be considered as active (aberrational, drug
  sensitive, etc..) in the pathway analysis. (default 0
  imposes no threshold) } \item{\code{min_gene_count}:}{
  \code{Object of class "numeric"} Analogous to
  min_gene_frequency -- the minimum number of patients with
  gene in active state for it to be considered active in
  the pathway analysis. (default 1 imposes no threshold) }
  \item{\code{individualEnrichment}:}{ \code{Object of
  class "logical"} A flag indicating if path analysis
  should be run for individual patients. (note, if
  individual patient analysis is run, min_gene_frequency
  and min_gene_count will automatically be set to 0 and 1,
  respectively, in the analysis of individual patients
  (though will be maintained in the analysis of any full
  cohort)) } \item{\code{.verbose}:}{ \code{Object of class
  "logical"} Flag indicating if extra output should be
  shown.  } \item{\code{.significanceTests}:}{ \code{Object
  of class "list"} Each slot contains a function. Each
  function implements the path test interface, accepting a
  PathSummaryRunner as an input, and returning a data frame
  with path names or identifiers as row names and values
  indicating pathway's significance (ex: p-values).  }
  \item{\code{settings}:}{ \code{Object of class "list"}
  The settings list object.  } \item{\code{targetname}:}{
  \code{Object of class "character"} One word description
  of the genes affected (ex: mutated) }
  \item{\code{dataSetName}:}{ \code{Object of class
  "character"} Description of the data set analyzed by the
  current arm.  } \item{\code{study}:}{ \code{Object of
  class "Study"} The study object for the current analysis.
  } \item{\code{.targetMatrix}:}{ \code{Object of class
  "matrix"} Bipartate graph matrix. This is a version of
  the pathway bipartate graph reduced to only contain genes
  considered "active" (aberrational or functionally
  sensitive). This is the primary matrix used in computing
  pathway significance.  } \item{\code{patientsum}:}{
  \code{Object of class "matrix"} Numeric values indicating
  the number of genes affected in each patient. Patient IDs
  are given as row names, values are the number of genes
  affected in each patient.  }
  \item{\code{path_summary_each_patient}:}{ \code{Object of
  class "list"} The results sets for each patient. List
  names are set as patient identifiers.  }
  \item{\code{gene_count_matrix}:}{ \code{Object of class
  "matrix"} Values are \code{"numeric"} Row names are gene
  identifiers. Values are the number of patients with
  affected gene.  } \item{\code{.gene_frequency_matrix}:}{
  \code{Object of class "matrix"} The frequency across
  cohort that each gene is active.  }
  \item{\code{.gene_vector}:}{ \code{Object of class
  "matrix"} Gene identifiers are given as row names. Values
  are logical, indicating if genes are considered 'active'
  in current cohort (or patient, for individual patient
  analysis).  } \item{\code{genomicnotpw}:}{ \code{Object
  of class "matrix"} Values are \code{"character"}
  identifiers for affected genes not found in pathways.  }
  \item{\code{active_genes_ea_path}:}{ \code{Object of
  class "matrix"} Values are pasted together
  \code{"character"} vectors of identifiers for genes found
  active in each pathway. Row names are given as pathway
  identifiers.  } \item{\code{pathsummary}:}{ \code{Object
  of class "data.frame"} containing main coverage and
  enrichment summary data for each affected pathway.  } }
}

