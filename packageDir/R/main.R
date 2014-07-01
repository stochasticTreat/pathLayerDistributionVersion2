#main.R

#'@title Load the basic set of analysis arms for use in a study
#'@description This function connects the basic set of analysis arms (drug screen analysis, somatic mutation, arbitrary gene input data and overlap analysis) to a Study object. 
#'@param STUDY A \code{Study} object.
#'@return A \code{Study} object with study arms loaded.
#'@export
#'@examples
#'study = getStudyObject(study.name="testDataSets", 
#' 											 geneIdentifierType="HUGO")
#'study = loadBasicArms(STUDY=study) 
#'sa = slot(object=study, name="arms")
#'print(names(sa))
#'print(slotNames(sa$functional_drug_screen_summary))
loadBasicArms<-function(STUDY){
	
	arms = STUDY@arms
	cat("\nLoading basic arms...\n")
	arms = loadDataArm(description="Load drug screen data",
										 title="functional_drug_screen_summary", 
										 mainFunction=packageDir:::RunDrugScreen, 
										 arms=arms)
	
	arms = loadDataArm(description="Load somatic mutation data",
										 title="somatic_mutation_aberration_summary", 
										 mainFunction=packageDir:::runSomaticMutationsProcessing, 
										 arms=arms)
	
	arms = loadDataArm(description="Load abitrary set of genes for path enrichment",
										 title="arbitrary_gene_data_input", 
										 mainFunction=packageDir:::RunGenericEnrichment, 
										 arms=arms)
	
	arms = loadDataArm(description="Run overlap analysis",
										 title="overlap_analysis", 
										 mainFunction=packageDir:::RunOverlapAnalysis, 
										 arms=arms)
	
	STUDY@arms = arms
	
	return(STUDY)
}


#'@title Run pathMod in interactive mode. 
#'@description Main function for interactive usage of package. Provides menu-based access for nearly all program functionality.
#'@param additionalArms Optional. Function loading additional data input arms. Function should take a Study object as an argument, use the loadDataArm() to load data arms, and return the Study object with arms added.
#'@export
#'@return The Study object created and filled by allInteractiveMainFunction.
#'@examples
#'\dontrun{
#'STUDY=allInteractiveMainFunction()
#'} 
########################################################################
allInteractiveMainFunction<-function(additionalArms=NULL){
	
	#initiate variables
	prevCohortFile = NULL
	arms = list()
	
	if(!exists("path_detail")){
		if(exists("STUDY")){
			cat("\nStudy object found in accessible envirnment\nUsing this Study\n")
			
			path_detail=STUDY@studyMetaData@paths
		}else{
			cat("\nNo Study object found, initilizing pathways.\n")
			path_detail = NULL
		}
	}

	########################################################################
	######## main menu selection logic
	########################################################################
	while(T){
		######################################################################
		######## establish the STUDY
		######################################################################
		initResults = F
		if(!exists("STUDY")){
			initResults=T
		}else if(STUDY@studyMetaData@studyName==""){
			initResults=T
		}#if(!exists("results"))
		
		if(initResults){
			STUDY = initiateStudy(path_detail=path_detail)
			cat("\nStudy initiated.. \n")
			# 		STUDY = getStudyObject(results=results)
			results = STUDY@results
			path_detail = getPaths(STUDY)
			cat("\nPaths set.. \n")
		}

		#establish path_detail, paths list object
		if(is.null(path_detail)){
			path_detail = getPaths()
			STUDY@studyMetaData@paths = path_detail
		}
		cat("\nReady to load basic study arms...\n")
		#establish arms
		STUDY = loadBasicArms(STUDY=STUDY)
		if(class(additionalArms)=="function"){
			print("loading additional arms")
			STUDY = additionalArms(STUDY=STUDY)
		}
		#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
		printProgramState(STUDY)
		
		########################################################################
		##############   main selection structure              #################
		########################################################################
		InputOptions = c("Change pathways", 
										 "Load settings from previous study",
										 "Run analysis from loaded settings",
										 armDescriptionList(study=STUDY))			
		optionframe1 = as.data.frame(matrix(1:length(InputOptions),ncol=1,dimnames=list(InputOptions,"Option number")))
		
		ProcessingOptions = c("combine aberration data and summarize by pathway", 
													"View summary of loaded data",
													"Compare sources of aberration data",
													"Create network diagrams for affected pathways",
													"Save a data summary to HTML",
													"Save current study", 
													"Change study name",
													"Clear all loaded settings",
													"Clear current study and study data", 
													"Run drug selection worksheet")
		optionframe2 = as.data.frame(matrix((length(InputOptions)+1):(length(InputOptions)+length(ProcessingOptions)),
																				ncol=1,
																				dimnames=list(ProcessingOptions,"Option number")))
		optionframe = rbind(optionframe1, 
												optionframe2)
		
		
		#########################################	prompt #########################################
		cat("\nThese are the options available at this time\n")
		cat("\nData input options:\n")
		print(optionframe1)
		cat("\nData processing options:\n")
		print(optionframe2)
		cat("\t\t\tTo quit program enter\t q\n")
		sel_line = readline("Please enter the option number: ")
		#########################################	selection logic #########################################
		if(sel_line=="q"){
			#### quitter 
			break
		}
		sel="Unknown input"
		if(sel_line%in%c("q",1:nrow(optionframe))) sel = rownames(optionframe)[as.integer(sel_line)]
		
		if(sel=="Run overlap analysis" && (!"functional_drug_screen_summary"%in%names(STUDY@results)|2>length(grep(pattern="_summary",x=names(STUDY@results))))){
			readline("\n\n***Drug screen and aberration data must be loaded before overlap is examined.***\n\n")
			sel="quit"
		}
		if(sel%in%armDescriptionList(study=STUDY)){#if this, then run data input arm that was initilized in armLoader.R
			cat("\nRunning arm to", sel, "\n")
			#run the arm corresponding to the selection
			STUDY=runArm(armDescription=sel, study=STUDY)
			results = STUDY@results
		}else if(sel=="Load settings from previous study"){
			STUDY = selectAndLoadSettings(study=STUDY)
		}else if(sel == "Change pathways"){#load paths
			path_detail=NULL
		}else if(sel=="Run analysis from loaded settings"){
			STUDY = autoRunFromSettings(study=STUDY)
			results = STUDY@results
			
		}else if(sel=="combine aberration data and summarize by pathway"){
			STUDY = combineAberrationTypes(study=STUDY)
		}else if(sel=="View summary of loaded data"){#View summary of loaded data 
			DataSummary(STUDY)
			stemp = readline("Summaries of all loaded data can be seen above.\nPress any key to continue.")
		}else if(sel=="Save current study"){
			saveStudy(study=STUDY, path="./output")
		}else if(sel=="Clear current study and study data"){
			saveprompt = readline("Would you like to save the current study before you clear it and its data? (y/n)")
			if(saveprompt=="y"){
				saveStudy(study=STUDY, path="./output")
			}
			results = NULL
			study_name = NULL
			STUDY@results=list()
			STUDY@studyMetaData@settings = list()
			STUDY@studyMetaData@studyName = ""
		}else if(sel=="Save a data summary to HTML"){
			#list summaries available to save
			# 		checkAndSaveStudy(study=STUDY)
			SaveToHTML(study=STUDY)
		}else if(sel=="Compare sources of aberration data"){
			print("Comparring..")
			STUDY = compareSources(study=STUDY)
			if("y"==readline("Would you like to make an HTML summary of the data type comparrison? (y/n)")){
				toHTML(table_list=list(comparisons=STUDY@results[["Aberration data type comparrison"]]),
							 limit_col="hyperg_p_w_FDR",
							 reorderTables=T,plimit=.05,maxrows=5000,
							 fname=paste("./output/",STUDY@studyMetaData@studyName,"/aberration_data_comparison.html",sep=""),
							 path_detail=path_detail)
			}
		}else if(sel=="Create network diagrams for affected pathways"){
			STUDY = addPathwayImagesWithSelection(study=STUDY)
		}else if(sel=="Change study name"){
			STUDY = changeStudyNamePrompt(study=STUDY)
		}else if(sel=="Clear all loaded settings"){
			STUDY = clearAllSettings(study=STUDY)
		}else if(sel=="Unknown input"){
			tmp = readline(paste("\nSorry, your input, \"", sel_line,"\", was not recognized.\nPress any key to return to the main menu.", sep=""))
		}else if(sel=="Run drug selection worksheet"){
			cat("...")
			#sfile=system.file("shinyDrugSelect/runShiny.R",package = "packageDir")
			print("Runing shiny main")
			#source(sfile, local=T)
			STUDY@results[["drugSelectionWorksheet"]] = runDrugWorksheet(STUDY=STUDY)
			cat("drugSelectionWorksheet found:",!is.null(STUDY@results$drugSelectionWorksheet))
		}
	}
	return(STUDY)
}


