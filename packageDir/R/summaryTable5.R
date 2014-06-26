
#save(stemp,file=testStudyObjectFileName, compress=T)
# 	somres = stemp@results$somatic_mutation_aberration_summary
# 	save(somres,
# 			 file=testStudyObjectFileNamesomatic, 
# 			 compress=T)

saveDefaultSettings<-function(psr){
	toSave =psr$settings
	fname="./reference_data/defaultSettings/defaultSummaryTableSettings.txt"
	dflist = saveSettings(set=toSave)
	write.table(x=dflist, file=fname, sep="\t")

}


#'@title Main pathway analysis function. 
#'@description The main function for conducting a pathway analysis on a set of gene data.
#'@param study A Study object. 
#'@param activeGeneDescription One word description of how the 'active' genes should be described (ex: mutated, amplified, drug_sensitive)
#'@param dataSetDescription Description of the data set provided for the pathway analysis. (ex: "Somatically mutated genes")
#'@param pgm The patient gene matrix for the patient data. Bipartate graph with row names set as gene names, column names set as patient IDs and values logical indicating if gene is considered "active" in patient.
#'@param originalDataMatrix The original data matrix. Makes this available for novel types of pathway analysis, not treating genes as discretely "on" or "off"
#'@param settings A settings list object controling how the path analysis will be conducted. Default given by call to getDefaultSettings()$defaultSummaryTable . 
#'@param coverage \code{character vector} of gene names. Should be provided if coverage provided by the analysis platform used to produce gene data entered into this function is not considered to have full genome coverage. 
#'@param coverageGeneDescription The description of genes that are covered be the analysis platform (ex: drug_targeted or sequenced)
#'@param coverageDataSetDescription A longer description of the set of genes provided in the coverage analysis
#'@param individualEnrichment logical. A flag indicating if pathway analysis should be run for individual patients. 
#'@export
#'@examples
#'studyObject = getStudyObject(study.name="testStudy", path_detail=getDefaultPaths())
#'
#'#conduct an analysis of pathway coverage only. 
#'targetlist  = c("NT5C2", "GUK1",  "ADAL",  "PCK1")
#'coverageAnalysis = summaryTable(study=studyObject, 
#'					coverageGeneDescription="testCoverageGene", 
#'					coverageDataSetDescription="test coverage analysis", 
#'					coverage=targetlist)
#'
#'#conduct full pathway analysis
#'pathwayAnalysis = summaryTable(study=studyObject, 
#'																pgm=getTestPGM(),
#'																dataSetDescription="test analysis of gene data",
#'																activeGeneDescription="analyzed_gene")
#'View(pathwayAnalysis$patientsums)
summaryTable<-function(study, 
											 activeGeneDescription="active_gene", 
											 dataSetDescription="analysis of undefined active genes",
											 pgm=NULL, 
											 originalDataMatrix=NULL,
											 settings=getDefaultSettings()$defaultSummaryTable, 
											 coverage=NULL,
											 coverageGeneDescription="covered_gene", 
											 coverageDataSetDescription="Coverage analysis", 
											 individualEnrichment=NULL){
	
	print("Inside summaryTable()")
	paths_detail = getPaths(path_file=study)
	#initilize the PathSummaryRunner, a reference object!!!
	
	psr = PathSummaryRunner$new(.verbose=F, path_summary_each_patient=list())
	
	psr$individualEnrichment = individualEnrichment
	psr$targetname = activeGeneDescription
	psr$dataSetName = dataSetDescription
	psr$study = study
	psr$coverage = coverage
	#prep everything for the pathway analysis
	covres = checkCoverage(psr=psr, 
												 paths_detail=paths_detail, 
												 geneDescription=coverageGeneDescription, 
												 dataSetDescription=coverageDataSetDescription)
	if(is.null(pgm)){
		cat("\nNo patient data found (ie, no pgm argument was passed to summaryTable).\nAssuming only platform coverage analysis was needed.\n")
		return(list(coverage_summary=covres$psr$coverage_summary))
	}
	psr=covres$psr
	paths_detail = covres$paths_detail
	cat("\nFound data for",paste(paste(dim(pgm),c("genes in","patients")), collapse=" "),"\n")
	##items to be set after the psr is used for coverage analysis
	psr$patientGeneMatrix = pgm
	psr = pathAnalysisSettings(study=study, 
														 s=settings, 
														 psr=psr, 
														 interactive=F)	#establish settings
	print("Path analysis settings established. . . ")
	print(dim(psr$patientGeneMatrix))
	

	if(!is.null(individualEnrichment)) psr$individualEnrichment = individualEnrichment
	
	psr$original_data_matrix = originalDataMatrix
	
	#assure correct formatting and data were supplied to summary table
	psr$patientGeneMatrix = psr$patientGeneMatrix==T

	#several preliminary computations
	print("Calculating patient sums")
	psr$patientsum = t((rep(x=1,times=nrow(psr$patientGeneMatrix))%*%psr$patientGeneMatrix)) #patientsum: matrix with the number of active genes found in each patient
	colnames(	psr$patientsum )<-"count_per_patient"
	
	#remove non-active genes from pgm
	print("Removing non-active genes from PGM")
	psr$gene_count_matrix = psr$patientGeneMatrix%*%rep(T, times=ncol(psr$patientGeneMatrix))
	psr$patientGeneMatrix = psr$patientGeneMatrix[psr$gene_count_matrix>0,,drop=F]
	
	checkSummaryTableInput(psr=psr, paths_detail=paths_detail)
	
	pathAnalysisResults = summaryTableInner(psr=psr, 
																					paths_detail=paths_detail)
	
	if(!is.null(coverage)) 	pathAnalysisResults$coverage_summary = psr$coverage_summary
	pathAnalysisResults$original_data_matrix = psr$original_data_matrix
	pathAnalysisResults$settings = psr$settings
	return(pathAnalysisResults)
}

individualPatientSummary<-function(psr, paths_detail){
	cat("\nRunny path summary for each individual patient. ")
	psr$path_summary_each_patient = list()
	psrTmp = psr$copy()
	if(!is.uninitilizedNull(psr$indPatientSummarySettings)) psrTmp$settings = psr$indPatientSummarySettings
	
	for(i in 1:ncol(psr$patientGeneMatrix)){#summarize by patient if multiple patients are passed
		pname = colnames(psr$patientGeneMatrix)[i]
		cat("\nPatient number",i,":", pname,"\n")
		psrTmp$patientGeneMatrix = psr$patientGeneMatrix[,i,drop=F]
	
		psrTmp$patientsum = t((rep(x=1,times=nrow(psrTmp$patientGeneMatrix))%*%psrTmp$patientGeneMatrix)) #patientsum: matrix with the number of active genes found in each patient
		colnames(	psrTmp$patientsum )<-"count_per_patient"
		#compress the patient gene matrix
		psrTmp$patientGeneMatrix = psrTmp$patientGeneMatrix[psrTmp$patientGeneMatrix%*%rep(T, ncol(psrTmp$patientGeneMatrix))>0,,drop=F]
		#re-make the target matrix
		psrTmp$.targetMatrix = getTargetMatrix(tgenes=rownames(psrTmp$patientGeneMatrix), paths=paths_detail$paths)
		#adjust basic settings
		psrTmp$individualEnrichment = FALSE
		psrTmp$min_gene_frequency=0
		psrTmp$min_gene_count=1
		psrTmp$dataSetName = paste("Patient",colnames(psrTmp$patientGeneMatrix),";",psr$dataSetName, sep=" ", collapse=" ")
		#add the results of enrichment analysis to the appropriate slot in the original pathSummaryRunner
		psr$path_summary_each_patient[[pname]]=summaryTableInner(psr=psrTmp, paths_detail=paths_detail)#individualPatientSummary(psr=psrTmp, paths_detail=paths_detail)
	}
	
	return(psr)
}

summaryTableInner<-function(psr, paths_detail){
	
	cat("\n----------------------------Entering path summary; patient records in current examination:",
			ncol(psr$patientGeneMatrix),"\n")
	cat("----------------------------genes in current examination:",
			nrow(psr$patientGeneMatrix),"\n")
	
	#Allows single patient record or set of patient records to be summarized, active gene 
	psr$.gene_vector = NULL#this will be a logical vector containing the genes considered "on" or active in the enrichment, names are the gene names
	psr$gene_count_matrix = NULL #For each gene, the numbers of samples in cohort that have that gene active
	psr$.gene_frequency_matrix = NULL #the frequency across cohort that each gene is active
	
	print( psr$individualEnrichment )
	print( class(psr$patientGeneMatrix) )
	
	if(ncol(psr$patientGeneMatrix) > 1 & psr$individualEnrichment){#if there are multiple patient samples, do the individual enrichments
		psr = individualPatientSummary(psr=psr, paths_detail=paths_detail)
	}
	
	#multiple columns of patient data passed as inputs; collapse them down into single column
	psr$gene_count_matrix = psr$patientGeneMatrix%*%rep(T, ncol(psr$patientGeneMatrix))
	colnames(	psr$gene_count_matrix )<-c("num_patients_with_affected_gene")
	psr$.gene_frequency_matrix = psr$gene_count_matrix/ncol(psr$patientGeneMatrix)
	
	if(psr$min_gene_frequency==0){
		psr$.gene_vector = psr$gene_count_matrix >= psr$min_gene_count
	}else{
		psr$.gene_vector = psr$.gene_frequency_matrix > min_gene_frequency
	}
	
	############### check that some genes are affected
	if(sum(psr$patientGeneMatrix)==0){ #true if there are no active genes in the current patients
		cpsdat = rep(1,times=nrow(psr$patientGeneMatrix))%*%psr$patientGeneMatrix
		curpatsum = matrix(data=cpsdat, ncol=1, dimnames=list(colnames(psr$patientGeneMatrix),"count_per_patient"))
		print("returning with no active genes in the current patient..")
		return(list(summarystats=generalSummary(psr=psr, paths_detail=paths_detail),
								patientGeneMatrix=psr$patientGeneMatrix, 
								patientsums=curpatsum, #return the number of active genes in the patient(s)
								genesummary=psr$gene_count_matrix))
	}#if no active genes

	if(psr$.verbose) cat(" .3.5 ")

	psr$genomicnotpw = setdiff(rownames(psr$patientGeneMatrix), 
												 colnames(paths_detail$paths))#genes not found in current patways
	psr$genomicnotpw = as.matrix(psr$genomicnotpw, ncol=1)

	if(psr$.verbose) cat(" st4.3.6 ")
	#find counts across pathways
	psr$.targetMatrix = getTargetMatrix(tgenes=rownames(psr$.gene_vector)[psr$.gene_vector], #gets only the genes that are 'active'
																		 paths=paths_detail$paths)#constricts pathways to set of targets found in patient gene matrix

	############## check that some paths are affected
	if(is.null(	psr$.targetMatrix) ) #if there are no targeted paths, then exit
	{
		cat("\n",psr$dataSetName,"was not found to have",psr$targetname,"genes annotated to the currently used set of pathways.\n")
		
		curpatsum = matrix(data=rep(1,times=nrow(psr$patientGeneMatrix))%*%psr$patientGeneMatrix,
											 ncol=1, 
											 dimnames=list(colnames(psr$patientGeneMatrix),"count_per_patient"))
		
		#if there are no targeted paths, the path enrichment matrix should be returned as a NULL value
		print("Returning with no targeted paths in current patient")
		return(list(summarystats=generalSummary(psr=psr, paths_detail=paths_detail),
								patientGeneMatrix=psr$patientGeneMatrix, 
								patientsums=curpatsum, #return the number of active genes in the patient(s)
								genesummary=psr$gene_count_matrix, 
					 			patientList = colnames(psr$patientGeneMatrix)))
	}

	################# get the pasted-together vectors of the active genes found in each path
	psr$active_genes_ea_path  =  getActiveGenesEachPath(psr=psr)
	if(psr$.verbose) cat(" st5.5.6: basic path info")
	
	################# build the pathSummary
	path_id = rownames(psr$.targetMatrix)
	tableOut = getBasicPathInformation(paths_detail=paths_detail, pathNames=path_id, psr=psr)
	
	if(psr$.verbose) cat(" st5.5.7: general paths affected info")
	#get general number of paths affected
	pathCounts = psr$.targetMatrix%*%rep(T, ncol(psr$.targetMatrix))
	colnames(pathCounts)<-paste(psr$targetname,"genes",sep="_")
	
	proportionAffected = pathCounts/tableOut$testable_path_length
	colnames(proportionAffected)<-paste0("proportion_",psr$targetname)
	
	if(psr$.verbose) cat(" st5.5.8: min max and freq info")
	minmaxfreq = getMinMaxAndFreq(psr=psr, paths_detail=paths_detail)
	
	psr$pathsummary = cbind.data.frame(tableOut, pathCounts, proportionAffected, minmaxfreq, stringsAsFactors=F)
	
	if(psr$.verbose) cat(" st5.5.9: sig tests..")
	if(!file.exists("./output/pathSummaryRunnerObject.rda")&(ncol(psr$patientGeneMatrix)>1)){
		dir.create(path="./output/", showWarnings=F, recursive=T)
		if(is.null(psr$path_summary_each_patient)) warning("path summaries for each patient is null... strange")
		warning("saving path summary runner for diagnostic purposes...")
		save(psr, file="./output/pathSummaryRunnerObject.rda")
	# 		path_summary_each_patient = psr$path_summary_each_patient
	# 		save(path_summary_each_patient, file="./output/pathSummaryRunnerObject_psep.rda")
	} 	
	
	for(tf in psr$.significanceTests){
		psr$pathsummary = cbind.data.frame(psr$pathsummary, tf(pathSigRunner=psr, paths_detail=paths_detail), stringsAsFactors=F)
		
	}
	psr$pathsummary = unfactorize(df=psr$pathsummary)
	
	psr$pathsummary = psr$pathsummary[order(psr$pathsummary[,5], decreasing=T), ,drop=F]
	
	print("Regular return from pathway analysis...")
	return(list(pathsummary=psr$pathsummary, #The path enrichment summary
							summarystats=generalSummary(psr=psr, paths_detail=paths_detail), #the overall summary stats
							patientGeneMatrix=psr$patientGeneMatrix, #the original patientGeneMatrix with rows removed if genes were not targeted
							patientsums=psr$patientsum, #counts of active genes for each patient
							genesummary=psr$gene_count_matrix, #count of patients who have activity in each gene
							active_genes_ea_path=psr$active_genes_ea_path, #list of genes that were found active in each path
							path_summary_each_patient=psr$path_summary_each_patient, #results from applying the summary table function
							active_genes_not_in_paths=psr$genomicnotpw, #list of active genes not in the pathways
							patientList = colnames(psr$patientGeneMatrix)))  #list of patient ids
}#summaryTable


