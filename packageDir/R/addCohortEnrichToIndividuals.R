


#psro = psr
#pdo
#function to allow results from cohort to be added to per patient results durring summaryTable execution, durring path significance tests.
addCohortToEachIndividual<-function(psro,pdo){
	if(!is.null(psro$path_summary_each_patient)){
		tmp = list(path_summary_each_patient = psro$path_summary_each_patient, pathsummary=psro$pathsummary)
	}
	emptydf = data.frame(matrix(ncol=0,nrow=nrow(psro$pathsummary)))
	return(emptydf)
}


#'@title Adds path enrichment information from cohort to paths found to contain active genes in individual patients.
#'@description Adds the columns in columnsToAdd from the whole cohort path enrichment to the path enrichments for each patient.
#'@param res The results list object containing slots named path_summary_each_patient and pathsummary.
#'@param columnsToAdd The columns from mainPathSum to add to each individual patient path summary
#'@param returnFlag Flag indicating if the list of patient results should be returned. 
#'@return The list of path summary results for each patient.
#'@examples
#' tstudy = getTestStudyObject()
#' all_results = slot(object=tstudy, name="results")
#' psep = addCohortToIndividuals(res=all_results$somatic_mutation_aberration_summary)
#' View(psep$p1$pathsummary)
addCohortToIndividuals<-function(res, columnsToAdd=grep(pattern="hyper",colnames(res$pathsummary))[1]:ncol(res$pathsummary) ){
	psep = res$path_summary_each_patient
	mainPathSum=res$pathsummary
	cat("Patients: ")
	for(pat in names(psep)){
		cat(pat,".. ")
		#for each patient, get the current
		#gets the set of paths affected in the patient
		cur = psep[[pat]]$pathsummary
		if(!nullOrNoRows(cur)){
			#retreive the paths in the cohort that are affected in the patient
			colsFromMain = mainPathSum[ppaths$path_id,columnsToAdd]
			colnames(colsFromMain) = paste0("cohort_",colnames(colsFromMain))
			psep[[pat]]$pathsummary = cbind.data.frame(ppaths, colsFromMain, stringsAsFactors=F)	
		}
	}
	cat("\n")
	return(psep)
}

nullOrNoRows<-function(dat){
	if(is.null(dat)) return(TRUE)
	if(!nrow(dat)) return(TRUE)
	return(FALSE)
}

# pathAnalysisFunctions = list(hypergeometric=hypergeometricPathEnrichment, addCohortToIndividual=addCohortToEachIndividual)