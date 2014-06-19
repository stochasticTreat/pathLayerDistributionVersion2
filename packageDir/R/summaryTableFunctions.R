# pathTestInterface<-function(){
# 	
# 	pathSigRunner = list(paths=paths,
# 										 original_data_matrix=original_data_matrix,
# 										 coverage=coverage,
# 										 gene_vector=gene_vector,
# 										 pathTargetMatrix=pathTargetMatrix)
# 	
# }


#ut1
#getActiveGenesEachPath()
#takes psr: the path summary runner object list
#returns: a 1 column matrix where each row is contains data for a single path and cell is a pasted together vector of gene names for its respective path
getActiveGenesEachPath <- function (psr) {

	tm = psr$.targetMatrix
	
	applyRes = apply(X=tm,MARGIN=1, FUN<-function(x){paste(colnames(tm)[x], collapse=" ")} )
	
	applyRes = matrix(applyRes, ncol=1, dimnames=list(names(applyRes)))
	
	return(applyRes)
	
}

#ut1
#used to get the first 3 columns of information going into the pathSummary
getBasicPathInformation <- function (paths_detail, pathNames, psr) {
	path_id = pathNames
	full_path_lengths = psr$study@studyMetaData@paths$full_path_length[path_id,,drop=FALSE]
	colnames(full_path_lengths)<-"full_path_length"
	testable_path_length = paths_detail$paths%*%rep(TRUE, times=ncol(paths_detail$paths))
	testable_path_length = testable_path_length[path_id,,drop=FALSE]
	return(cbind.data.frame(path_id, full_path_lengths, testable_path_length))
}

#ut1
generalSummary<-function(psr, paths_detail){
	cat("\n..Producing general summary of path coverage..")
	#summary set is the general summary table output in the list returned by this function
	summarySetLabels = c("Data set used:",
											 "Path analysis of genes that are:",
											 "Path source:", 
											 "Patient data sets:",
											 "Number of pathways examined:",
											 paste("Pathways containing",psr$targetname,"genes:"),
											 "Unique genes in all pathways:",
											 paste("Unique genes",psr$targetname,"in cohort"), 
											 paste("Unique,",psr$targetname,"genes from cohort, found in pathways:"), 
											 paste("Unique,",psr$targetname,"genes not in pathways:"), 
											 paste("Mean, median, min and max ",psr$targetname," genes in cohort's patients:"))
	defaultSummarySetValues = c(psr$dataSetName, #"Data set used:"  dataSetName,
															psr$targetname, #"Path analysis of genes that are:",  psr$targetname,
															paths_detail$file, #"Path source:",   	paths_detail$info,
															ncol(psr$patientGeneMatrix), #"Patient data sets:",  	ncol(patientGeneMatrix),
															nrow(paths_detail$paths), #"Number of pathways examined:", nrow(paths),
															ifelse(test=is.uninitilizedNull(psr$.targetMatrix),yes=0,no= nrow(psr$.targetMatrix)), # paste("Pathways containing",psr$targetname,"genes:"), 	nrow(pathTargetMatrix),
															ncol(paths_detail$paths),# 	"Unique genes in all pathways:", 				ncol(paths), 
															nrow(psr$patientGeneMatrix),# 	paste("Unique genes",psr$targetname,"in cohort"),  nrow(patientGeneMatrix),
															ifelse(test=is.uninitilizedNull(psr$.targetMatrix), yes=0, no=ncol(psr$.targetMatrix)),# 	paste("Unique,",psr$targetname,"genes from cohort, found in pathways:"), 
															ifelse(test=is.uninitilizedNull(psr$genomicnotpw), yes=0, no=nrow(psr$genomicnotpw)),# 	paste("Unique,",psr$targetname,"genes not in pathways:"), 		nrow(genomicnotpw),
															ifelse(test=is.uninitilizedNull(psr$patientsum), 
																		 yes=paste("Mean:", 0, "median:", 0,"min:", 0, "max:", 0),
																		 no=paste("Mean:", round(mean(psr$patientsum),digits=2),
																		 				 "median:", median(psr$patientsum),
																		 				 "min:", min(psr$patientsum),
																		 				 "max:", max(psr$patientsum))))
	
	covSumOut = cbind.data.frame(summarySetLabels, defaultSummarySetValues, stringsAsFactors=FALSE)
	rownames(covSumOut) <- NULL
	colnames(covSumOut) <- c("Data", "Value")
	cat("done..\n")
	return(covSumOut)
}

is.uninitilizedNull<-function(val){
	
	return(class(val)=="uninitializedField"|is.null(val))

}

# save(s, file=defaultSettingsFile)
#allow interactive determination of path analysis settings
pathAnalysisSettings <- function (psr, study, s=NULL, interactive=FALSE, 
																	defaultSettingsFile=NULL) {
	
	if(is.data.frame(s)) s = dfToList(df=s)
	
	print("inside pathAnalysisSettings")
	interactiveTmp = s$interactive

	if(is.null(s)){
		s = list()
	}else if(!is.list(s)){
		print(s)
	}

	if(interactiveTmp){
		interTest = s[["Use special path significance analysis settings for this data type? (y/n)"]]
		if(!is.null(interTest)){
			if(interTest=="y") s$interactive = TRUE
		}
	}else{
		s[["Use special path significance analysis settings for this data type? (y/n)"]] = FALSE
	}
	
	etPrompt = "To add additional path significance tests, please select the R script file containing the interfaces to these enrichment tests"
	s = setting(s, prompt="Use special path significance analysis settings for this data type? (y/n)")
	specialSettingsFile = s$.text == "y"
	
	s$interactive = s$.text == "y"
	if(!specialSettingsFile){
		print("Loading default path analysis settings..")
		#get default settings
		if(is.null(study@studyMetaData@settings$defaultSummaryTable)){
			stmp = loadSettings()
		}else{
			stmp = study@studyMetaData@settings$defaultSummaryTable
		}
		s = c(s, stmp)
		s = s[unique(names(s))]
	}

	
	#limit by proportion cohort with active gene
	s = setting(s, prompt="To limit by proportion of cohort with genes in examined state, enter minimum proportion here:\n(enter 0 for no limit): ")
	psr$min_gene_frequency=as.numeric(s$.text)
	#limit by count of patients 
	s = setting(s, prompt="To limit by count of patients with genes in examined state, enter count here:\n(Enter 1 to skip limiting by count)")
	psr$min_gene_count = as.numeric(s$.text)
	#
	s = setting(s, prompt="Analyze pathways for individual members of the cohort? (y or n) ")
	psr$individualEnrichment=s$.text=="y"
	
	s = setting(s, prompt="Would you like to run the path analysis in verbose mode? (y or n) ")
	psr$.verbose = s$.text=="y"
	
	if(specialSettingsFile){
		#select file of enrichment tests	
		s = setting(s, 
								requireInput=FALSE, 
								prompt=etPrompt)
	}

	if(is.null(s[[etPrompt]])){
		s[[etPrompt]] = ""
		s$.text = ""
	} 
	
	psr$.significanceTests = addSignificanceTests(fnames=s[[etPrompt]])
	
	s$interactive = interactiveTmp
	
	psr$settings = s
	

	return(psr)
}


#adds paths significance test from the file names provided in the argument, fnames
addSignificanceTests<-function(fnames, defaultTest = hypergeometricPathEnrichment){
	if(is.null(fnames)) return(c(defaultTest))
	if(!length(fnames)|fnames=="") return(c(defaultTest))
	ptests = c()
	for(n in fnames){
		source(n, local=TRUE)
		ptests = c(ptests, pathAnalysisFunctions)
	}
	return(ptests)
}

#'@title loadPathSigTests
#'@param sigTestFile An R script file containing one or more functions which implement the pathSigTest interface and a line defining a named list of variable name 'pathAnalysisFunctions' with names describing the path analysis functions and values set as the respective path analysis functions.
#'@return A named list containing the pathway analysis functions to be applied. 
loadPathSigTests<-function(sigTestFile='./hypergeometricPathAnalysis.R'){
	
	pathAnalysisFunctions = list()
	if(!is.null(sigTestFile)){
		source(sigTestFile, local=TRUE)
	}
	pathAnalysisFunctions[["hypergeometric"]] = hypergeometricPathEnrichment
	return(pathAnalysisFunctions)
}


checkSummaryTableInput<-function(psr,paths_detail){
	pass=TRUE
	if(!sum(psr$patientGeneMatrix)&ncol(psr$patientGeneMatrix)>1){
		readline(paste("Error, it appears there were no",targetname,
									 "genes in the cohort of data entered,\n though this may not be the case. Please check your data and your patient identifiers.\n Press enter to continue."))
		pass=FALSE
	}
	if(!is.matrix(psr$patientGeneMatrix)){
		tmp=readline(prompt="ERROR!! The patient gene matrix provided to the summarytable4 function was not of the R, matrix data type.\nPress any key to continue")
		pass=FALSE
	}
	pdnames = c("paths","source","date","info", "full_path_length")
	allgood=TRUE
	for(n in pdnames) if(!length(paths_detail[[n]])) allgood=FALSE
	if(!allgood){
		tmp=readline(paste("ERROR! The paths list is missing one or more of the following named elements: \n", 
											 paste(pdnames, collapse="; "),
											 "\nPress any key to continue.",
											 collapse=" "))
		pass=FALSE
	}
	print(ifelse(test=pass, yes="Path analysis inputs look correct..", no="Issues found with path analysis inputs.."))
}

adjustPathsForCoverage<-function(paths_detail, targetMatrix){
	
	pd = paths_detail$copy()
	pd$paths=targetMatrix
	pd[["gene_overlap_counts"]] = rep(TRUE,nrow(	pd$paths))%*%	pd$paths
	pd[["full_path_length"]] = 	pd$paths%*%rep(TRUE,ncol(	pd$paths))
	return(pd)	
}


#checkCoverage
#checks if a set of coverage genes was provided; if so, 
#coverage analysis is conducted and appended to psr
checkCoverage<-function(psr, paths_detail, 
												geneDescription=NULL, dataSetDescription=NULL){
	cat("\nChecking coverage ...\n")
	if(is.null(psr$coverage)){
		cat("..no coverage limitation established.\n")
		return(list(psr=psr, paths_detail=paths_detail))
	} 
	cat("\nAdding coverage analysis...\n")
	###set the target matrix
	psr$.targetMatrix = getTargetMatrix(tgenes=psr$coverage, 
																			paths=paths_detail$paths)
	###conduct a summaryTable analysis for the coverage
	#note: psr is a reference object, but in the summaryTable function, all the appropriate settings will be established after this
	psrCtmp = psr$copy()
	psrCtmp$min_gene_frequency = 0
	psrCtmp$min_gene_count = 1
	psrCtmp$individualEnrichment = FALSE
	psrCtmp$.significanceTests = c()
	psrCtmp$patientGeneMatrix = PGMFromVector(genevector=psr$coverage)
	#assure correct formatting and data were supplied to summary table
	psrCtmp$patientGeneMatrix = psrCtmp$patientGeneMatrix==TRUE
	
	#several preliminary computations
	psrCtmp$patientsum = t((rep(x=1,times=nrow(psrCtmp$patientGeneMatrix))%*%psrCtmp$patientGeneMatrix)) #patientsum: matrix with the number of active genes found in each patient
	colnames(	psrCtmp$patientsum )<-"count_per_patient"
	
	psrCtmp$dataSetName = ifelse(test=is.null(dataSetDescription), yes=paste("coverage for:",psr$dataSetName), no=dataSetDescription)
	psrCtmp$targetname = ifelse(test=is.null(geneDescription), yes="covered", no=geneDescription)
	
	psr$coverage_summary = summaryTableInner(psr=psrCtmp, paths_detail=paths_detail)
	psr$coverage_summary$pathsummary = psr$coverage_summary$pathsummary[,1:5]
	
	psr$coverage_summary$summarystats = psr$coverage_summary$summarystats[1:(nrow(psr$coverage_summary$summarystats)-1),] #remove the min median min and max row
	paths_detail = adjustPathsForCoverage(paths_detail=paths_detail, targetMatrix=psrCtmp$.targetMatrix)
	
	rm(psrCtmp)
	
	return(list(psr=psr, paths_detail=paths_detail))
}

getMinMaxAndFreq<-function(psr = psr, paths_detail = paths_detail){ # pathTargetMatrix,patientGeneMatrix,genesum, paths){
	#gets the minumumn, maximum and frequency of mutations in path genes, across the cohort
	#takes: 		pathTargetMatrix: bipartate graph, path matrix reduced to only the pathways and genes that are targeted by the patients in the cohort
	#						patientGeneMatrix: bipartate graph, indicates which patients have which genes active, columns=patient ids, rows=gene names, values=TRUE/FALSE
	#						paths: the bipartate graph of the full pathways
	#						genesum: data frame with three columns: types, the gene names
	#																										counts, the number of patients the genes are found active in
	#																										rownames(genesum) = genesum[,"types"]
	#returns: 	list with three elements: maxmuts: the maximum number of patients inwhich any one of the pathway's genes is found active
	#																			minmuts: the minimum number of patients inwhich any one of the pathway's genes is found active
	#																			pfreq: the frequency of active genes in the pathway, across the cohort
	#

	pathTargetMatrix = psr$.targetMatrix
	patientGeneMatrix = psr$patientGeneMatrix
	genesum = data.frame(types=rownames(psr$gene_count_matrix), 
											 counts=psr$gene_count_matrix, 
											 stringsAsFactors=FALSE)
	colnames(genesum)<-c("types", "counts")

	paths = paths_detail$paths
	
	gsnames = c(rownames(genesum), 
							setdiff(colnames(paths),
											rownames(genesum))) #[!colnames(paths)%in%genesum[,"counts"]])
	
	tmpgenesum = rep(0,times=length(gsnames))
	names(tmpgenesum)<-gsnames
	tmpgenesum[genesum[,"types"]] = genesum[,"counts"]
	
	names(tmpgenesum)<-gsnames
	maxmuts = rep(0,nrow(pathTargetMatrix)) #will hold the max number of active genes across the cohort for any gene in pathway
	minmuts = rep(0,nrow(pathTargetMatrix)) #will hold the minimum number of active genes across the cohort for any gene in pathway
	pfreq = rep(0,nrow(pathTargetMatrix)) #sum the number of patients across the cohort with active genes in any of the path genes
	p_act_count = rep(1,nrow(pathTargetMatrix))

	for(i in 1:nrow(pathTargetMatrix)){#for each path in the target matrix
		pname = rownames(pathTargetMatrix)[i]
		p = paths[pname,,drop=FALSE]#pull out one pathway
		mems = colnames(p)[p]#pull out the set of path memebers
		xm = patientGeneMatrix[intersect(x=mems,y=colnames(pathTargetMatrix)),,drop=FALSE]#extract just the parts of the patientGeneMatrix corresponding to the path
		p_act_count[i] = sum(xm)
		pathpatients = (rep(TRUE, nrow(xm))%*%xm) > 0 #find those patients with active genes in any of the genes in the pathway
		pfreq[i] = sum(pathpatients) #sum the number of patients across the cohort with active genes in any of the path genes
		maxmuts[i] = max(tmpgenesum[mems])#find the maximum number of active genes of any of the genes in the pathway
		minmuts[i] = min(tmpgenesum[mems])#find the minimum number of active genes of any of the genes in the pathway
	}
	
	pnames = rownames(pathTargetMatrix)
	
	plens = paths_detail$full_path_length[pnames,]
	
	datout = cbind.data.frame(pfreq/nrow(psr$patientsum),
														pfreq, 
														p_act_count, 
														maxmuts, 
														minmuts, 
														stringsAsFactors=FALSE)
	
	colnames(datout)<-c(paste0("proportion_of_cohort_w_",psr$targetname,"_gene_in_path"),
											paste0("Num_patient_with_",psr$targetname,"_gene_.s._in_path"), 
											paste0("count_of_",psr$targetname,"_genes_in_path_cross_cohort"),
											"max_in_one_gene", 
											"min_in_one_gene")
	
	rownames(datout)<-pnames
	
	return(datout)

}


