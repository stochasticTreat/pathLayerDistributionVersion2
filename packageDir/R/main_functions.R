#main functions

clearAllSettings<-function(study){
	study@studyMetaData@settings = list()
	for(n in names(study@results)){
		study@results[[n]]$settings = list()
	}
	return(study)
}

changeStudyNamePrompt<-function(study){
	name = ""
	while(T){
		nname = readline("Please enter a new study name: ")
		if(nname != "") break
	}
	study = changeStudyName(study=study, newName = nname)
	return(study)
}

changeStudyName<-function(study, newName){
	study@studyMetaData@studyName = newName
	if(grepl(pattern="^test", x=newName)){
		study@studyMetaData@RootFile = paste("./output/", newName, "/", sep="", collapse="")
	}else{
		study@studyMetaData@RootFile = paste("./output/study_", newName, "/", sep="", collapse="")
	}
	return(study)
}

#autoRunFromSettings
#runs all analyses described in settings found in meta data
#to make the data file, run the code below in examples, 
#make sure the working directory is set to 
#setwd("../packageDir/")
#then run this
# abacavirSettings = sres@studyMetaData@settings
# save(abacavirSettings, file="./inst/extdata/abacavirSettings.rda")
#'@title Automatically run set of analyses from a set of loaded settings. 
#'@description Function searches for settings objects in the STUDY-studyMetaData-settings slot, then automatically executes the analysis for each arm described by the settings. 
#'@param study A \code{Study} object with settings available. 
#'@param verbose A flag indicating if analyses should be run in interactive/verbose mode, prompting user to check intermediate settings. 
#'@return A \code{Study} object with analysis described by settings complete. 
#'@export
#'@examples
#'#load a set of settings
#'load(system.file("extdata/abacavirSettings.rda",package="packageDir"), verbose=TRUE)
#'#initialize the study
#'s1 = getStudyObject(path_detail=getDefaultPaths(), 
#'										settings=abacavirSettings,
#'										study.name="integrationTestAvacavirMetabolism")
#'s1 = loadBasicArms(STUDY=s1)
#'sres = autoRunFromSettings(study=s1)
autoRunFromSettings<-function(study, verbose=T){
	#pull each analysis type out
	allanalyses = names(study@studyMetaData@settings)	
	analyses = allanalyses[!grepl(pattern="^default|^overlap_analysis", x=allanalyses, ignore.case=T)]#remove any settings starting with the string "default"
	
	if(length(analyses)){
		cat("Settings being run for these intput arms:\n")
		cat(analyses, "\n")
		for(a in analyses){
			study = runArm(armDescription=a, study=study, 
										 fromDescription=F, interactive=F)
		}
		if("overlap_analysis"%in%allanalyses){
			study = runArm(armDescription="overlap_analysis", study=study, fromDescription=F, interactive=F)
		}
	}else if(verbose){
		tmp = readline("Sorry, no settings could be found.\nSelect 2 from main menu to load settings from another study.\nPress enter to continue. ")
	}
	return(study)
}#autoRunFromSettings

checkForceRowNames<-function(tab){
	if(sum(rownames(tab) != as.character(1:nrow(tab)))&!is.null(rownames(tab))) return(tab)
	rn = tab[,1]
	tab = tab[,2:ncol(tab),drop=F]
	rownames(tab)<-rn
	return(tab)	
}



#runArm
#function that executes the armsMain interface
#if NULL is returned by the arm, NULL will be assigned to the slot in the study and settings
#takes: 
#		study: the study object
#		arm description: the arm description string or the arm name (if arm name, the fromDescription flag must be set to F)
#		interactive: flag for interactive mode
#returns: study
#'@title runArm
#'@description Main function for taking and running a study arm. Provides needed settings and study object, implementing the study arm interface. 
#'@param armDescription the desription slot from a DataArm object
#'@param study A Study object.
#'@param fromDescription Indicates if arm should be retrieved by its description or by its title. 
#'@param interactive Indicates if arm should be run in interactive mode. 
#'@return Returns a Study object with results from running the study arm. 
#'@export
runArm<-function(armDescription, 
								 study, 
								 fromDescription=T, 
								 interactive=NULL){

	#pull out settings and arms
	arms = study@arms
	settings = study@studyMetaData@settings
	
	#establish which arm's settings to extract
	if(fromDescription){
		armName = arms$dictionary[armDescription]
	}else{
		armName = armDescription
	} 
	
	curSettings = settings[[armName]]
	if(is.data.frame(curSettings)) curSettings = dfToList(df=curSettings)
	
	if(is.null(interactive)&!is.null(curSettings)){
		inter_option = readline("Would you like to process the data interactively so that\nnew settings or files can be selected? (enter y or n)\n")
		interactive = inter_option=="y"
		if(interactive){
			changeFiles = readline("Would you like to only change the files and leave all the settings the same? (y/n)")
			if(changeFiles=="y"){
				interactive=F
				curSettings$changeFiles = T
			}
		}
	}
	
	if(is.null(curSettings)){
		curSettings = list()
		curSettings$interactive=T
		interactive = T
	}else{
		if(is.null(interactive)) interactive = T
		curSettings$interactive = interactive 
	}
	
	#pull out arm's main function
	armMain = arms[[armName]]@mainFunction
	
	#run the arm's main function
	res = armMain(settings=curSettings, study=study)
	
	if(!is.null(res)){
		if(is.null(res$resTypeName)){
			study@results[[armName]] = res
			study@studyMetaData@settings[[armName]] = res$settings
		}else{#if the arm is generic gene data input arm, there will be something in the resTypeName slot
			study@results[[res$resTypeName]] = res$summary
			study@studyMetaData@settings[[res$resTypeName]] = res$settings
		}
	}else{
		cat("\nNothing returned from analysis. Returning to main menu.\n")
	}
	return(study)
}

#printProgramState
#displays a summary of loaded data
#takes: results list, study_name string and path_detail list object
#returns: nothing
#'@title printProgramState
#'@description Displays to standard out general information about the state of the program, including pathways, settings and data currently loaed. 
#'@param stud A Study object from which the program state is to be determined. 
#'@export
printProgramState<-function(stud){#results, study_name, path_detail){
	results = stud@results
	study_name = studyName(stud)
	path_detail = getPaths(stud)
	########################################################################
	##############   Current state of loaded data          #################
	########################################################################
	cat("\n#########################################################")
	cat("\nStudy name:",study_name,"\n")
	cat("#########################################################")
	if(is.null(path_detail)){
		cat("\n***No paths loaded***\n")
	}else{
		cat("\nCellular pathway set from: ", path_detail$info, "containing ", as.integer(ncol(path_detail$paths)), "genes across", as.integer(nrow(path_detail$paths)), "pathways.\n" )
	}
	cat("######### Currently loaded data processing settings:")
	cat("\n")
	if(length(stud@studyMetaData@settings)){
		cat(paste("-",names(stud@studyMetaData@settings), "\n", sep=""), sep="")
	}else{
		cat("\n***There are currently no data processing settings loaded***\n")
	}
	cat("######### Currently loaded data summaries:\n")
	###show loaded data
	if(length(results)){
		
	# 		summaries:
	# 		summis = grepl(pattern="summary", x=names(results))
	# 		sumnames = names(results)[summis]
	# #			processed data:
	# 		odnames = names(results)[!summis]
		sumnames = names(results)
		for(n in sumnames){
			
			curnamev = strsplit(n,split="_")[[1]]
			curname = paste(curnamev[1:(length(curnamev)-1)])
			cursum = results[[n]]
			if(length(grep(pattern="drug",x=curname,ignore.case=T))){
				cat("-",curname,"data loaded.\t", as.character(nrow(cursum$coverage_summary$pathsummary)),"paths analyzed\n")
			}else if(grepl(pattern="^overlap", x=n,ignore.case=T)){
				cat("-",curname,"data loaded.\t", as.character(nrow(cursum$"Aberrationally enriched, containing drug targets")),"paths functionally targeted and significantly aberrational\n")
			}else{
				cat("-",curname,"data loaded.\t", as.character(nrow(cursum$pathsummary)),"paths affected\n")
			}
		}
		
	# 		cat("######### Other loaded data:\n")
	# 		for(n in odnames) cat("-",gsub(pattern="_",replacement=" ",x=n),"\n")
	}else{
		cat("\n***There is currently no aberration or drug screen data loaded***\n")
	}
	cat("\n#########################################################\n")
	# 	cat("#########################################################")
	# 	###indicated if patients will be limited to some subset
	# 	if(is.null(aberration_patient_subset)){
	# 		cat("\nUsing full set of patients\n")
	# 	}else if(length(grep(pattern="TCGA",x=patient_set_name,ignore.case=T))){
	# 		cat("\nSingle patient loaded, pid:",patient_set_name,"\n")
	# 	}else{
	# 		cat("\nLoaded",patient_set_name,"cohort consisting of ",length(aberration_patient_subset)," patients out of ",length(all_patients)," patients.\n")
	# 	}
	# 	cat("#########################################################\n")
}#printProgramState

# darkPaths = results$overlap_analysis$'Aberration enriched, not drug targeted'$path_id
#determine list of genes found most often in pathways
BangForBuck<-function(darkPaths, path_detail){
	cat("\nIn bang for buck analysis...")
	if(!is.vector(darkPaths)){
		darkPaths = darkPaths[,1]
	}
	pths = path_detail$paths
	#1 extract the paths from path_detail
	xpaths = pths[darkPaths,,drop=F]
	genecounts = rep(1, nrow(xpaths))%*%xpaths
	genecounts = genecounts[1,]
	genecounts = genecounts[genecounts!=0] #after this length(genecounts)=440 
	genecounts = genecounts[order(genecounts, decreasing=T)]
	
	pathcross = c()
	for(gn in names(genecounts)){
		# 		print(gn)
		cur = paste(rownames(xpaths)[xpaths[,gn]], sep=" | ", collapse=" | ")
		pathcross = c(pathcross, cur)
	}
	
	pathsPerGene = cbind.data.frame(genecounts, pathcross, stringsAsFactors=F)
	colnames(pathsPerGene)<-c("Number of paths", "Path names")
	print("..done\n")
	return(pathsPerGene)
}

#next, report which of the pathways interact with 

#'@title Display summaries of loaded data. 
#'@description This function prints to the screen summary statistics for data sets from each data arm. 
#'@param study A \code{Study} object. 
#'@export
#'@examples
#'study=getTestStudyObject()
#'DataSummary(study=study)
DataSummary<-function(study){
	results=study@results
	cat("\nDisplaying summaries for all loaded data:\n")
	sumin = grep("summary", x=names(results),ignore.case=T)

	for(i in sumin){
		titleTxt = gsub(pattern="summary", replacement="", x=names(results)[i])
		titleTxt = gsub(pattern="_", replacement=" ", x=titleTxt)
		if(length(grep(pattern="coverage", names(results[[i]]),ignore.case=T))){
			cat("\nOverview of",titleTxt, "coverage:\n")
			print(results[[i]]$coverage_summary$summarystats)
		}
		cat("\nOverview of",titleTxt,"data:\n")
		print(results[[i]]$summarystats)
	}
}

#options(warn=-1)
bmerge<-function(x,y){
	if(is.null(x)){
		return(y)
	}
	#merge by rowname
	#make rownames rownames again
	#
	out=merge(x=x,y=y,by=0,all=T,incomparable=NaN)
	rownames(out)<-out[,1]
	out = out[,2:ncol(out)]
	#fill in NAs with zeros: 
	return(out)
}


#'@title Allows comparrison of summary statistics for all loaded aberration data types. 
#'@description The function examines and compares the pathway summaries for all loaded aberration data types
#'@param study A \code{Study} object. More than one aberration data type needs to have been loaded before this function is run. 
#'@return A \code{Study} object with a slot named, \code{"Aberration data type comparrison"} in the results list filled with a table comparring impacts of each aberration data type for each pathway.
#'@export
#'@examples
#' stud = getTestStudyObject()
#' stud = compareSources(study=stud)
#' resslot = slot(object=stud, name="results")
#' View(resslot$'Aberration data type comparrison')
compareSources<-function(study){
	#takes:		 the results list
	#returns:	 Two item list: list(outtable=outtable, results=results)
	#							outtable slot: 	table containing comparisons of aberration data
	#							results slot:		the original results list object
	if(is.null(study@results[["combined_aberrations_summary"]])){
		study@results[["combined_aberrations_summary"]] = combineAberrationTypes(study=study)
	}
	results = study@results
	cat("\nSelecting columns for output table..\n")
	#columns in output: <the overall analysis> <number of nodes contributed from each> <pvalue from each
	#put first 10 columns of the  onto the outtable
	outtable = results$combined_aberrations_summary$pathsummary[,c(1,2,3,4,12)]
	colnames(outtable)[4:ncol(outtable)] = paste(colnames(outtable)[4:ncol(outtable)], "combined", sep="_")
	
	#go through the aberration summaries and 
	#append on the number of aberrational genes and the hyperg p w/FDR
	abi = grep(pattern="_aberration_",x=names(results))
	for(i in abi){ # for each slot in results holding aberration data
		cur = results[[i]]
		extract = cur$pathsummary[,c(1,4,12)]
		#rename the extract columns
		cname = gsub(pattern="_aberration_summary", replacement="", x=names(results)[i])
		colnames(extract)[2:ncol(extract)] = paste(colnames(extract)[2:ncol(extract)], cname, sep="_")
		#now append
		outtable = merge(outtable,
										 extract,
										 by="path_id",
										 all=F, 
										 sort=F)
	}
	cat("\nReordering by p. aberrational...\n")
	#reorder outtable by the proportion aberrational: 
	outtable = outtable[order(outtable[,5],decreasing=F),]
	cat("\nReturning from compareSources()\n")
	
	cat("\nComparrison of aberration sources complete. \n")
	
	study@results = results
	study@results[["Aberration data type comparrison"]] = outtable
	
	return(study)
}#compareSources

combineAberrationTypes2<-function(study, s, results, overlapPatients=NULL){
	results=study@results
	if(VERBOSE) print("inside combineAberrationTypes2()....................................")
	
	if(is.null(s))	s = pickSettings(study=study, settingType="aberration")
	
	big_pgm = NULL
	abtypenames = "Summary of"
	if(is.null(s)) s=study@studyMetaData@settings$defaultSummaryTable
	
	#get the indexes for the aberration data
	ri = grep(pattern="aberration_",x=names(results)) 
	
	if(!length(ri)&length(grep(pattern="aberrations_", x=names(results)))&(length(overlapPatients)==1)){
		#if we're looking at a single patient and they already have a combined aberrations analysis done, return it
		return(results$combined_aberrations_summary)
	}
	
	if(VERBOSE){
		print(overlapPatients)
		print(names(results))
	}
	#next: find the overlap of patients in each ab type
	
	# 	#if overlap patients are not yet provided , check if there are overlap patients
	# 	if(is.null(overlapPatients)){
	# 		for(i in ri){
	# 			curres = results[[i]]
	# 			print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	# 			if(VERBOSE) print(names(results)[i])
	# 			if(VERBOSE) print(curres$patientsums)
	# 			if(is.null(overlapPatients)){
	# 				overlapPatients = rownames(curres$patientsums)
	# 				if(VERBOSE) print("overlap patients rout 1:")
	# 				if(VERBOSE) print(overlapPatients)
	# 			}else{
	# 				numover = length(overlapPatients)
	# 				overlapPatients = intersect(overlapPatients, rownames(curres$patientsums))
	# 				nlost = numover - length(overlapPatients)
	# 				if(VERBOSE) print("overlap patients rout 2:")
	# 				if(VERBOSE) print(overlapPatients)
	# 				if(nlost > 0){
	# 					cat("Note: ",nlost,"patients were found not to have been analyzed by all platforms.\n")	
	# 				}
	# 				
	# 			}
	# 		}
	# 	}#if(is.null(overlapPatients))
	
	if(VERBOSE) cat("Total patients analyzed using all platforms:",length(overlapPatients),"\n")
	if(!length(overlapPatients)){
		cat("No patients were found to be analyzed on all platforms.\n",
				"Re-run analysis using fewer aberration data types or data sets\n",
				"with matching patient IDs to examine aberration overlap or \n",
				"combination of aberration types. If an error does not occur,",
				"any results that follow may be spurious.\n")
		readline("Press any key to continue.\n")
	}
	
	coverage_summary=NULL
	coverageSet = NULL

	if(VERBOSE) print(ri)
	#Establish joint coverage of all aberration types
	#for each aberration type
	cat("\nNote, limiting coverage to inlude only those genes\nanalyzed on all aberration analysis platforms\n")
	for(i in ri){
		curname = names(results)[i]
		if(VERBOSE){
			print("for(i in ri){ iteration:")
			print(i)
			print(curname)
		}
		curres = results[[i]]
		if("coverage_summary"%in%names(curres)){
			if(VERBOSE) cat(" ..coverage summary found.. ")
			#there is a coverage summary, so extract the coverage set
			covsum=curres$coverage_summary
			if(is.null(coverageSet)){
				if(VERBOSE) print("coverageSet is null")
				coverageSet = rownames(covsum$genesummary)
			}else{
				if(VERBOSE) print("coverageSet is not null")
				coverageSet = intersect(x=coverageSet, y=rownames(covsum$genesummary))
			}
		}
		abtypenames = paste(abtypenames,curres$summarystats[1,2], sep=" ")
		if(VERBOSE) print("abtypenames:")
		if(VERBOSE) print(abtypenames)
		trimGeneMatrix=curres$patientGeneMatrix[,overlapPatients,drop=F]#make sure only patients who have representation in each data type are being represented here
		big_pgm = as.matrix(joinTables(t1=big_pgm, t2=trimGeneMatrix, name_prefix=curname, fill=0))
		cat(".")
	}#for
	
	if(VERBOSE) print("infunction 2")
	if(!is.null(coverageSet)){
		#1: trim the big_pgm so there are no un-covered genes
		goodPGMRows = rownames(big_pgm)%in%coverageSet
		big_pgm = big_pgm[goodPGMRows,,drop=F]
	}
	
	if(VERBOSE){
		print("summaryTable from combineAberrationTypes, regular")
		cat("dim pgm going into summaryTable:", dim(big_pgm), "\n")
		print(typeof(big_pgm))
		print(ncol(big_pgm))
		print(nrow(big_pgm))
		print(length(big_pgm))
	}
	
	if( !is.null(big_pgm)&!is.null(coverageSet) ){
		combined_aberrations_summary = summaryTable(study=study, individualEnrichment=TRUE,
																								activeGeneDescription="aberrational",
																								pgm=big_pgm, 
																								coverage=coverageSet,
																								settings=s,
																								coverageDataSetDescription="analysis of coverage by poor-coverage genomic aberration analysis platform(s)",
																								coverageGeneDescription="analysis of coverage by poor-coverage genomic aberration analysis platform(s)", 
																								dataSetDescription=abtypenames)
	}else if(!is.null(big_pgm)){
		combined_aberrations_summary = summaryTable(study=study, 
																								individualEnrichment=TRUE,
																								activeGeneDescription="aberrational",
																								pgm=big_pgm, 
																								settings=s,
																								dataSetDescription=abtypenames)
	}else{
		warning("Cannot find patient gene matrix combined for all patients...")
	}
	
	return(combined_aberrations_summary)	
}#combineAberrationTypes2()

combineAberrationTypes<-function(study, s=NULL, results=NULL, runEachPatient=F, overlapPatients=NULL){
	if(is.null(results)) results = study@results
	if(is.null(s)){
		s = pickSettings(study=study, settingType="aberration")
	}
	cat("\nCombining aberration types ....................................\n")
	big_pgm = NULL
	abtypenames = "Summary of"
	if(is.null(s)) s=study@studyMetaData@settings$defaultSummaryTable
	#get the indexes for the aberration data
	ri = grep(pattern="aberration_",x=names(results)) 
	
	if(!length(ri)&length(grep(pattern="aberrations_", x=names(results)))&(length(overlapPatients)==1)){
		#if we're looking at a single patient and they already have a combined aberrations analysis done, return it
		return(results$combined_aberrations_summary)
	}
	
	if(VERBOSE){
		print(overlapPatients)
		print(names(results))
	}
	#next: find the overlap of patients in each ab type

	#if overlap patients are not provided , check if there are overlap patients
	if(is.null(overlapPatients)){
		for(i in ri){
			curres = results[[i]]
			cat("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n")
			if(VERBOSE) print(names(results)[i])
			if(VERBOSE) print(curres$patientsums)
			if(is.null(overlapPatients)){
				overlapPatients = rownames(curres$patientsums)
				if(VERBOSE) print("overlap patients rout 1:")
				if(VERBOSE) print(overlapPatients)
			}else{
				numover = length(overlapPatients)
				overlapPatients = intersect(overlapPatients,rownames(curres$patientsums))
				nlost = numover - length(overlapPatients)
				if(VERBOSE) print("overlap patients rout 2:")
				if(VERBOSE) print(overlapPatients)
				if(nlost > 0){
					cat("Note: ",nlost,"patients were found not to have been analyzed by all platforms.\n")	
				}
				
			}
		}
	}#if(is.null(overlapPatients))
	
	if(VERBOSE) cat("Total patients analyzed using all platforms:",length(overlapPatients),"\n")
	if(!length(overlapPatients)){
		cat("No patients were found to be analyzed on all platforms.\n",
				"Re-run analysis using fewer aberration data types or data sets\n",
				"with matching patient IDs to examine aberration overlap or \n",
				"combination of aberration types. If an error does not occur,",
				"any results that follow may be spurious.\n")
		readline("Press any key to continue.\n")
	}
	
	coverage_summary=NULL
	coverageSet = NULL
	tm=NULL
	print(ri)
	#Establish joint coverage of all aberration types
	#for each aberration type
	cat("\nNote, limiting coverage to inlude only those genes\nanalyzed on all aberration analysis platforms\n")
	for(i in ri){
		curname = names(results)[i]
		if(VERBOSE){
			print("for(i in ri){ iteration:")
			print(i)
			print(curname)
			print("point x")
		}
		curres = results[[i]]
		if("coverage_summary"%in%names(curres)){
			if(VERBOSE) print("coverage summary found")
			#there is a coverage summary, so extract the coverage set
			covsum=curres$coverage_summary
			if(is.null(coverageSet)){
				if(VERBOSE) print("coverageSet is null")
				coverageSet = rownames(covsum$genesummary)
			}else{
				if(VERBOSE) print("coverageSet is not null")
				coverageSet = intersect(x=coverageSet, y=rownames(covsum$genesummary))
			}
		}
		abtypenames = paste(abtypenames,curres$summarystats[1,2], sep=" ")
		if(VERBOSE) print("abtypenames:")
		if(VERBOSE) print(abtypenames)
		trimGeneMatrix=curres$patientGeneMatrix[,overlapPatients,drop=F]#make sure only patients who have representation in each data type are being represented here
		big_pgm = as.matrix(joinTables(t1=big_pgm, t2=trimGeneMatrix, name_prefix=curname, fill=0))
		if(VERBOSE) print("point z")
	}#for
	
	if(VERBOSE) print("infunction 2")
	if(!is.null(coverageSet)&is.null(big_pgm)){
		#1: trim the big_pgm so there are no un-covered genes
		goodPGMRows = rownames(big_pgm)%in%coverageSet
		big_pgm = big_pgm[goodPGMRows,,drop=F]
	}
	
	if(VERBOSE){
		print("summaryTable from combineAberrationTypes, regular")
		cat("dim pgm going into summaryTable:", dim(big_pgm), "\n")
		print(typeof(big_pgm))
		print(ncol(big_pgm))
		print(nrow(big_pgm))
		print(length(big_pgm))
	}
	
	if(!(is.null(big_pgm)&is.null(coverageSet))){
		combined_aberrations_summary = summaryTable(study=study, individualEnrichment=runEachPatient,
																								activeGeneDescription="aberrational",
																								pgm=big_pgm, 
																								coverage=coverageSet,
																								settings=s,
																								coverageDataSetDescription="analysis of coverage by poor-coverage genomic aberration analysis platform(s)",
																								coverageGeneDescription="analyzed", 
																								dataSetDescription=abtypenames)
	}else{
		combined_aberrations_summary = list()
	}
	
	return(combined_aberrations_summary)	
}#combineAberrationTypes()

#results$overlap_analysis$combined_aberrations_summary$path_summary_each_patient$X08.00173$

getPatientSubset<-function(){
	
	############################################
	############################################
	###############          load patient subset
	prompt = paste("Enter f to load a subset of patients from a file\n",
								 "Enter s to load a subset of patients from a clinical data table\n",
								 "Enter r to reset and use all patients: ",sep="")
	pidsource = readline(prompt)
	if(pidsource=="r"){ #reset and use all patients
		aberration_patient_subset = NULL
	}
	if(pidsource=="f"){#load a subset of patients from a file
		cat("\nPlease select file of patient identifiers\n")
		fname=file.choose()
		pids = read.delim(file=fname,sep="\t", quote=F, stringsAsFactors=F)
		pids = unlist(lapply(pids, extract_pid))
		aberration_patient_subset=pids
	}else if(pidsource=="s"){#load a subset of patients from a clinical data table
		cat("\nPlease select the clinical data file. The first column must contain patient IDs.\n")
		fname=file.choose()
		clinical = read.delim(file=fname,sep="\t", stringsAsFactors=F)
		pidcol = 1
		all_patients = unlist(lapply(clinical[,pidcol], extract_pid))
		extract_col = NULL
		while(T){
			extract_col = readline("Please enter the name of the column over which to divide the patient IDs: ")	
			if(extract_col%in%colnames(clinical)){
				tmp_patient_set_name = extract_col
				break
			}
			cat("Sorry, that column name did not match any of the columns")
		}
		
		summed_col = summarize_by(col=clinical[,extract_col],display=F)
		cat("\nHere are the values found in that column:\n")
		print(summed_col)
		val = NULL
		while(T){ #### Find the column to sort on. 
			val = readline("Please enter the value in that column you would like to extract patient IDs by: ")	
			if(val%in%summed_col[,"types"]){
				tmp_patient_set_name = paste(tmp_patient_set_name, val, sep="_", collapse="_")
				break
			}#if
			cat("Sorry, that value name did not match any of the values")
		}##while
		###############################################################################################
		### with column in extract_col, and taking the pids with val in extract_col, pull the trigger, 
		### and put the set of pids into aberration_patient_subset
		ecol = clinical[,extract_col]
		ecol[is.na(ecol)] = paste(val, 0, sep="_", collapse="_")#this is a hack to deal with NAs
		aberration_patient_subset = all_patients[ecol==val]
		cat("The subset of patients selected is", tmp_patient_set_name, 
				". This subset contains", length(aberration_patient_subset), "patients.")
		if("y"==readline("Is this the correct cohort? (please enter y or n) ")){
			patient_set_name=tmp_patient_set_name
			break
		}##if
	}
	return(aberration_patient_subset)
}#getPatientSubset
# 
