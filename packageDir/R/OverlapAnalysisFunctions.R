#overlapAnalaysisFunctions


RunOverlapAnalysis<-function(settings, study){
	
	results = abDrugOverlapAnalysis(study=study,  
																	settings=settings)
	return(results)
}


test.abDrugOverlapAnalysis<-function(){
	
	#stud = STUDY
	#test all instances with individuals and cohorts. 
	#test matched cohort
	
	#get default overlap analysis settings
	defSettings = getTestStudySettings()
	overlapSettings = defSettings$overlap_analysis
	#get a study object pre-loaded with data aberration and functional data. 
	stud = getTestStudyObject()
	#run the overlap analysis with functional genomic and genomic aberration data from the same cohort
	overlapRes1 = abDrugOverlapAnalysis(study=stud, settings=overlapSettings)

	#simulate the analysis with functional genomic and genomic aberration data from sepparate cohorts
	studySepparateCohort = stud
	patientSums = studySepparateCohort@results$functional_drug_screen_summary$patientsums
	rownames(patientSums) <- toupper(rownames(patientSums))
	
	studySepparateCohort@results$functional_drug_screen_summary$patientsums  = patientSums
	
	overlapRes2 = abDrugOverlapAnalysis(study=studySepparateCohort, 
															settings=overlapSettings)
	
	
}


#settings slots expected:
#	combinedAberrations
#'@title abDrugOverlapAnalysis
#'@description main function for running analysis of functional and genomic aberration arms
#'@param study A study object (must have results loaded for aberration and functional analyses)
#'@param thresh Numeric, the threshold used to establish which pathways will be considered to be drug sensitive and aberrational.
#'@param enrich_col String, the name of the column in the path summary used to determine significance of a pathway. This is the column to which the threshold argument is to be applied.
#'@param verbose Logical, a flag indicating if the program is to provide extra diagnostic output.
#'@param settings List obejct. This contains the settings used in the overlap analysis. Use getDefaultOverlapSettings()
#'@return Study object with overlap analysis in results list added or replaced 
#'@export
#'@examples
#'#get a study object pre-loaded with genomic aberration and functional genomic data.
#' stud = getTestStudyObject()
#' 
#'#get default overlap analysis settings
#'defSettings = getTestStudySettings()
#'overlapSettings = defSettings$overlap_analysis
#' 
#'#run the overlap analysis
#'
#'overlapRes1 = abDrugOverlapAnalysis(study=stud, settings=overlapSettings)
#' 
#'\dontrun{
#'#If settings for the overlap analysis cannot be found, the overlap 
#'#analysis will run in interactive mode, requiring user input.
#'#This line of code will run the overlap anaysis in interactive mode, 
#'#requesting user input to establish settings:
#'
#'abDrugOverlapAnalysis(study=stud)
#'}
#'@importFrom VennDiagram venn.diagram
#'@importFrom calibrate textxy
abDrugOverlapAnalysis<-function(study,
																thresh=0.05,
																enrich_col="hyperg_p_w_FDR",
																verbose=T,
																settings=NULL){
	if(!exists("VERBOSE")) VERBOSE=T
	cat("\n******************************************************************************************\n")
	cat("\n*************************************overlap analysis*************************************\n")
	cat("\n******************************************************************************************\n")
	
	sfolder = paste(studyFolder(s=study),"/results/", sep="")
	#set up ola object
	ola = list() #ola = overlap analysis
	ola$results = study@results
	ola$path_detail = getPaths(study)
	ola$threshold = thresh
	ola$enrich_col = enrich_col

	ola$outfname_path = sfolder
	ola$verbose = verbose

	if("interactive"%in%names(settings)){ 
		#if interactive is found as a flag in settings it is the first time overlap analysis has been run for the study
		#put it in the right place
		if( "basic_overlap"%in%names(settings) ){#overlap analysis has already been run
			settings_full = settings #transfer the settings list
			settings_full$interactive = NULL #erase the flag from the master list
			interactiveFlag = settings$interactive # put the interactive flag in the right place
			settings = settings_full$basic_overlap # in the 
			settings$interactive = interactiveFlag
		}else{#it's the first time the overlap was run on the current study
			settings_full = list()
			settings_full$basic_overlap = settings
		}
	}else{
		settings_full = settings
		settings = settings_full$basic_overlap
	}
		
	# 	settings$interactive = settings_full$interactive
	# 	settings_full$interactive = NULL

	settings = setting(s=settings, prompt="Would you like to run the function/aberration overlap analysis on each patient? (y/n)")
	runEachPatient=settings$.text == "y"
	ola$runEachPatient = runEachPatient

	while(T){
		settings = setting(s=settings, "Would you like to use default path analysis settings or use those associated with study arms? (please enter d or a)")
		if(settings$.text%in%c("a","d")) break
		print("Sorry, your input was not understood, please try again.")
	}
	if(settings$.text=="d"){
		settings_full$defaultSummaryTable = getDefaultSettings()$defaultSummaryTable
		settings$interactive = F
	}else{
		settings_full$defaultSummaryTable = NULL
		settings$interactive = T
	}
	
	settings_full$basic_overlap = settings
	
	ola$settings = settings_full
	
	ola = selectOverlapType(ola=ola, study=study)
	
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	ola$results$overlap_analysis$settings = ola$settings

	return(ola$results$overlap_analysis)
}


selectOverlapType<-function(ola, study){
	
	#determine if matched
	overlapPatients = CheckPatientOverlap(results=ola$results)
	
	#Assure only the sensitive patients are 
	overlapPatients = adjustOverlapToSensitive(overlap=overlapPatients, 
																						 functionalEnrichmentAnalysis=ola$results$functional_drug_screen_summary)
	
	matched = length(overlapPatients)>0 #if there are overlap patients, it's assumed to be a matched sample
	
	if(matched){
		#1: cohort, matched
		#limit cohort
		ola = runMatchedCohort(overlapPatients=overlapPatients, ola=ola, study=study)
		
	}else{
		cat("\nNote: no patient sample ids were found to be associated with\nboth sensitive gene targets and aberration analysis/analyses.\n")
		#2: cohort, unmatched
		#no limiting, just run it once the combination is done
		ola = runUnmatchedCohort(ola=ola, study=study)
		
	}
	
	if(ola$runEachPatient){#if the overlap analysis should be run for each patient
		
		ola$outfname_path = paste(ola$outfname_path, "overlap_analysis_each_patient/", sep="")
		
		if(matched){
			#3: individual patient, matched
			#pass drug screen summary for patient plus aberration summary for patient
			ola = runMatchedIndividuals(overlapPatients=overlapPatients, ola=ola, study=study)
		}else{
			#4: individual patient, unmatched
			#pass drug screen summary for individual patient, plus aberration summary for whole cohort
			ola = runUnmatchedIndividuals(ola=ola, study=study)
		}
		#indresults should be the list with each slot named with a patient identifier
	}
	
	return(ola)
}

runUnmatchedCohort<-function(ola, study, functional_analysis_name="functional_drug_screen_summary"){
	if(VERBOSE) cat("\nola$settings:\n") 
	if(VERBOSE) print(ola$settings)
	functionalEnrichmentAnalysis = ola$results[[functional_analysis_name[1]]]
	#case 1: use what was passed as settings
# 	combSettings = ola$settings$combinedAberrations
	if(is.null(ola$settings$defaultSummaryTable)){#if the default settings aren't provided, look for the settings in the arms

		ola$settings$combinedAberrations = NULL
		#case 2: use the default settings
		#case 3: use the settings from an aberration arm (taken care of in combineAberrationTypes)
		#					---> if the user said not to use the default, the combineAberrationTypes2 will make the user pick
	}else{
		ola$settings$combinedAberrations = ola$settings$defaultSummaryTable
	}
	
	if(is.null(ola$results$combinedAberrations)) ola$results$combinedAberrations = combineAberrationTypes(study=study, 
																																																			 s=ola$settings$combinedAberrations, 
																																																			 results=ola$results, 
																																																			 runEachPatient=F)
	ola$settings$combinedAberrations = ola$results$combinedAberrations$settings
	
	combinedAberrations = ola$results$combinedAberrations 
	
	ola$outfname_path = paste(ola$outfname_path, "overlap_analysis/", sep="")
	ola$results$overlap_analysis = coreOverlapAnalysis(functionalEnrichmentAnalysis=functionalEnrichmentAnalysis, 
																										 combinedAberrations=combinedAberrations, 
																										 ola=ola)
	
	return(ola)
}#runUnmatchedCohort

pickSettings<-function(study, settingType = "aberration"){
	si = grep(pattern=settingType, x=names(study@studyMetaData@settings), ignore.case=T)
	targetArmSettings = names(study@studyMetaData@settings)[si]
	if(length(si)>1){
		sres = settingList(s=list(interactive=T), 
								prompt="Please pick a study arm whose path analysis settings you would like to use.", 
								set=names(study@studyMetaData@settings)[si])
		targetArmSettings = sres$.text
	}
	return(study@studyMetaData@settings[[targetArmSettings]])
}

runMatchedCohort<-function(ola, overlapPatients, study){
	
	# 	if(is.null(ola$settings$combinedAberrations)) ola$settings$combinedAberrations = 
	#case 1: use the settings passed
	funSettings = ola$settings$functionalPathAnalysis
	if(is.null(ola$settings$defaultSummaryTable)){#if case 1 is null
		
		ola$settings$functionalPathAnalysis = NULL
		#case 2: use the default settings
		#case 3: use the settings from an aberration arm (taken care of in combineAberrationTypes)
		#					---> if the user said not to use the default, the combineAberrationTypes2 will make the user pick
	}else{
		ola$settings$functionalPathAnalysis = ola$settings$defaultSummaryTable
	}
	cat("\n\nAnalysis of functional drug screen for patients who have had both aberration and functional analysis..\n")
	functionalEnrichmentAnalysis = extractedDrugScreenAnalysis(psubset=overlapPatients, 
																														 study=study,
																														 s=ola$settings$functionalPathAnalysis,
																														 results=ola$results, 
																														 path_detail=ola$path_detail, 
																														 eachPatient=ola$runEachPatient) #holds the results for the functional enrichment analysis
	
	ola$results$overlap_analysis$functional_drug_screen_summary = functionalEnrichmentAnalysis
	
	ola$settings$functionalPathAnalysis = functionalEnrichmentAnalysis$settings
	
# 	overlapPatients = intersect(overlapPatients, colnames(functionalEnrichmentAnalysis$patientsums))
	
	if(is.null(ola$results$combinedAberrations)){
		
		#case 1: use the settings passed
		if(is.null(ola$settings$defaultSummaryTable)){#if the default settings aren't provided, look for the settings in the arms
			
			ola$settings$combinedAberrations = NULL
			#case 2: use the default settings
			#case 3: use the settings from an aberration arm (taken care of in combineAberrationTypes)
			#					---> if the user said not to use the default, the combineAberrationTypes2 will make the user pick
		}else{
			ola$settings$combinedAberrations = ola$settings$defaultSummaryTable
		}
		cat("\n\nAnalysis of combined aberration data for patients who have had both aberration and functional analysis..\n")
		
		ola$results$combinedAberrations = combineAberrationTypes2(study=study, 
																															s=ola$settings$combinedAberrations, 
																															results=ola$results, 
																															overlapPatients=overlapPatients)
		
	} 
	
	ola$settings$combinedAberrations = ola$results$combinedAberrations$settings
	
	
	ola$outfname_path = paste(ola$outfname_path, "overlap_analysis/", sep="")
	ola$results$overlap_analysis = coreOverlapAnalysis(functionalEnrichmentAnalysis=functionalEnrichmentAnalysis, 
																										 combinedAberrations=ola$results$combinedAberrations, 
																										 ola=ola)

	return(ola)
}

runMatchedIndividuals<-function(ola, overlapPatients, study){
	
	# 	functionalEnrichmentAnalysis = NULL #holds the results for the functional enrichment analysis
	# 	combinedAberrations = NULL #holds the results for the combined aberration enrichment analysis
	pol = list() #patient overlap analyses
	
	tmpAbs = ola$settings$combinedAberrations
	#pull out the appropriate summaries
	if(is.null(tmpAbs)){
		tmpAbs = ola$results$combinedAberrations$settings
		ola$settings$combinedAberrations = tmpAbs
	}
	tmpAbs$"To limit by count of patients with genes in examined state, enter count here:\n(Enter 1 to skip limiting by count)" = "1"
	tmpAbs$'To limit by proportion of cohort with genes in examined state, enter minimum proportion here:\n(enter 0 for no limit): '="0"
	
	tmpDs  = ola$settings$functional_drug_screen_summary
	if(is.null(tmpDs)){
		tmpDs = ola$settings$defaultSummaryTable
	# 		tmpDs = ola$results$functional_drug_screen_summary$settings
	# 		if(class(tmpDs)!="list") tmpDs = dfToList(df=tmpDs)
	# 		tmpDs$interactive = ola$settings$basic_overlap$interactive
		ola$settings$functional_drug_screen_summary = tmpDs
	}
	
	tmpDs$"To limit by count of patients with genes in examined state, enter count here:\n(Enter 1 to skip limiting by count)" = "1"
	tmpDs$'To limit by proportion of cohort with genes in examined state, enter minimum proportion here:\n(enter 0 for no limit): '="0"
	
	fpathtmp  = ola$outfname_path
	for(p in overlapPatients){
		print(p)
		ola$outfname_path=paste(fpathtmp, p,"/overlap_analysis/", sep="")
	
		combsum = combineAberrationTypes2(results=ola$results, 
																			study=study, 
																			s=tmpAbs,
																			overlapPatients=p)
		
		functionalSum = extractedDrugScreenAnalysis(psubset=p, 
																								eachPatient=F,
																								study=study, 
																								s=tmpDs,
																								results=ola$results, 
																								path_detail=ola$path_detail)
		
		#run the overlap analysis
		pres = coreOverlapAnalysis(functionalEnrichmentAnalysis=functionalSum, 
															 combinedAberrations=combsum, 
															 ola=ola, 
															 coverage=ola$results$functional_drug_screen_summary$coverage_summary)
		
		#put the results in the right place
		pol[[p]] = list()
		pol[[p]][["overlap_analysis"]] = pres
	}
	ola$results$overlap_analysis$overlap_analysis_each_patient = pol
	return(ola)
}#runMatchedIndividuals

runUnmatchedIndividuals<-function(ola, study){
	
	functionalEnrichmentAnalysis =  ola$results$overlap_analysis$functional_drug_screen_summary #holds the results for the functional enrichment analysis
	
	combinedAberrations = ola$results$overlap_analysis$combined_aberrations_summary #holds the results for the combined aberration enrichment analysis

	#get the set of individuals to be analyzed
	patientSet = rownames(ola$results$functional_drug_screen_summary$patientsums)
	ola$results$overlap_analysis$overlap_analysis_each_patient = list()
	
	pol=list()
	
	fpathtmp  = ola$outfname_path
	for(p in patientSet){
		print(p)
		
		ola$outfname_path=paste(fpathtmp, p,"/overlap_analysis/", sep="")
		#extract patient record, extract patient drug results	
		patientFunctionalAnalysis = functionalEnrichmentAnalysis$path_summary_each_patient[[p]]
		
		pres = coreOverlapAnalysis(functionalEnrichmentAnalysis=patientFunctionalAnalysis,
															 combinedAberrations=combinedAberrations,
															 ola=ola, 
															 coverage=ola$results$functional_drug_screen_summary$coverage_summary)
		pol[[p]] = list()
		pol[[p]][["overlap_analysis"]] = pres
	}
	ola$results$overlap_analysis$overlap_analysis_each_patient = pol
	return(ola)
}#runUnmatchedIndividuals

CheckPatientOverlap<-function(results){
	cat("..checking patients with functional and aberrational analysis..")
	psetOut = c()
	#check if analysis will be limited to patients who are in both cohorts, the drugs screen and the aberration sets
	funIndex = grep("functional", x=names(results))
	psetOut = rownames(results[[funIndex]]$patientsums)
	
	abIndex = grep("aberration_", x=names(results))
	for(i in abIndex){
		psetOut = intersect(x=psetOut, y=rownames(results[[i]]$patientsums))
	}
	cat(" ..checked...")
	return(psetOut)
}#CheckPatientOverlap

plotOverlapArea<-function(tabout_targAndAb, results, abandtargplotname){
	
	if(!nrow(tabout_targAndAb)){
		print("There are no pathways in the overlap between drug targeted/drug sensitive and aberrational")
		return()
	}
	#######################
	#### plot overlap area
	if(VERBOSE) print("Plotting overlap area.. ")
	if(VERBOSE) print(colnames(tabout_targAndAb))
	png(filename=abandtargplotname,width=800,height=500)
	
	maxpathlength = max(tabout_targAndAb[,3])
	cexvals = (tabout_targAndAb[,3]/maxpathlength)*10
	
	plot(cex=cexvals,
			 ylim=c(0,1.2*max(tabout_targAndAb[,"drug_targeted_genes"])),#)),
			 x=tabout_targAndAb[,"proportion_of_cohort_w_aberrational_gene_in_path"], 
			 y=tabout_targAndAb[,"drug_targeted_genes"],#"drug_targeted_genes"], 
			 main="Proportion affected vs number of drug targets.",
			 xlab="Proportion of cohort with aberrational gene in path", 
			 ylab="Number of drug targets")
	textxy(tabout_targAndAb[,"proportion_of_cohort_w_aberrational_gene_in_path"], 
				 tabout_targAndAb[,"drug_targeted_genes"], 
				 labs=1:nrow(tabout_targAndAb))#,
	#cx=1.3)
	legend(x="topright", 
				 legend=paste("Size = path length (",paste(range(tabout_targAndAb[,3]), collapse=" to "),"nodes)\nNumber = path number in the \'enriched and targeted\' table"), 
				 cex=.8, pch=1)
	dev.off()
	return()
}#plotOverlapArea

adjustOverlapToSensitive<-function(overlap, functionalEnrichmentAnalysis){
	#get the total sensitive genes per patient
	psums = functionalEnrichmentAnalysis$patientsums
	#extract the patient ids for those patients with active genes
	keepers = rownames(psums)[psums[,1]>0]
	#subset the overlap patients to only those patients with sensitive genes
	newOverlap = intersect(overlap, keepers)
	return(newOverlap)
}#adjustOverlapToSensitive

extractedDrugScreenAnalysis<-function(psubset, results, path_detail, study, s, eachPatient){
	#extractedDrugScreenAnalysis
	#runs path analysis on drug screen patient subset
	#returns drug screen results set with coverage analysis appended
	#first get the functional analysis out of results
	cat("...inside extractedDrugScreenAnalysis()...")
	findex=grep(pattern="functional", names(results),ignore.case=T)
	if(length(findex)>1) warning("More than one functional analysis was found in the input to the overlap analysis. (ex: more than one drug screen) This is currently a functionality of the program which has not been tested!")
	if(VERBOSE) print(names(results))
	if(VERBOSE) print(findex)
	fa = results[[findex]]
	pgm = as.matrix(fa$patientGeneMatrix[,psubset,drop=F]) #extract the subset of patients from the patient gene matrix
	covset = rownames(fa$coverage_summary$genesummary)
# 	tm = getTargetMatrix(tgenes=covset, paths=path_detail$paths)
	if(VERBOSE) print("summaryTable from extractedDrugScreenAnalysis")
		
	cat("\ndim pgm going into summaryTable:", dim(pgm), "\n")
	
	if(!is.null(pgm) & (nrow(pgm)>0)){
		
		if(is.null(s)) s = fa$settings
			
		summ = summaryTable(study=study, 
												coverageGeneDescription="targeted",
												coverageDataSetDescription="Coverage of functional analysis platform",
												individualEnrichment=eachPatient,
												coverage=covset,
												settings=s,
												pgm=pgm, 
												activeGeneDescription="functionally_affected", 
												dataSetDescription="Functional path enrichment in patient subset for function-aberration overlap analysis")
	}else{
		print("There were no sensitive genes in this patient")
		summ = list()
	}
	summ$coverage_summary = fa$coverage_summary
	return(summ)
}#extractedDrugScreenAnalysis

#'@title The central function used in an overlap analysis
#'@description Examines the actual overlap between functionally significant and aberrationally signficant pathways.
#'@param functionalEnrichmentAnalysis The functional enrichment analysis to be used. 
#'@param combinedAberrations The aberration analysis to be used
#'@param ola The overlap analysis runner object to be used. 
#'@param coverage If coverage is limited, this will contain the set of genes the coverage is limited to. 
#'@return An overlap analysis
#'@import VennDiagram
coreOverlapAnalysis<-function(functionalEnrichmentAnalysis, combinedAberrations, ola, coverage=NULL){
	
	if(is.null(coverage)) coverage=ola$results$functional_drug_screen_summary$coverage_summary
	
	dir.create(path=ola$outfname_path, showWarnings=F, recursive=T)
	outlist = list()
	outlist[["imageSlots"]] = list()
	drugcoverage = coverage$pathsummary[,"path_id"]
	if(nullOrEmpty( functionalEnrichmentAnalysis$pathsummary )){
		senspaths = NULL
		sensen=NULL
		print("There are no sensitive targets in this sample.")
	}else{
		print("Sensitive targets found.")
		senspaths = functionalEnrichmentAnalysis$pathsummary[,"path_id"]
		sensen  = functionalEnrichmentAnalysis$pathsummary[functionalEnrichmentAnalysis$pathsummary[,ola$enrich_col]<ola$threshold,"path_id"]
	}
	
	# 	if(is.null(senspaths)){
	# 		print("There are no sensitive targets in this sample.")
	# 		sensen=NULL
	# 	}else{
	# 		#T3			paths enriched for drug-sensitive genes
	# 		print("Sensitive targets found.")
	# 		sensen  = functionalEnrichmentAnalysis$pathsummary[functionalEnrichmentAnalysis$pathsummary[,ola$enrich_col]<ola$threshold,"path_id"]
	# 	}
	
	#vennList[["Paths enriched for drug-sensitve targets"]]
	#T4			aberration_coverage
	abpathsum = combinedAberrations$pathsummary
	if(VERBOSE) print(names(combinedAberrations))
	if(VERBOSE) print(is.null(abpathsum))
	aben = abpathsum[abpathsum[,ola$enrich_col]<ola$threshold,"path_id"] #aberration enriched
	
	if(VERBOSE) print("Creating venn diagram..")
	
	vennList = list()
	vennFileName= "Overlap_venn_diagram.png"
	fullvennFileName = paste(ola$outfname_path,vennFileName,sep="")
	outlist$imageSlots[[vennFileName]] = vennFileName
	
	vennList[["Drug-screen targeted paths"]] = drugcoverage
	vennList[["Aberration enriched paths"]] = aben
	# 	print("search():")
	# 	print(search())
	# 	print("ls()")
	# 	print(ls())
	require(VennDiagram)
	if(is.null(functionalEnrichmentAnalysis)){#ola$results$functional_drug_screen_summary$pathsummary)){#if this is null, then there is only a drug coverage analysis available
		print("venn option 1 -- functionalEnrichmentAnalysis is NULL")
		vennTmp = venn.diagram(vennList, 
													 filename=NULL,
													 main.cex=1.5,
													 cex=1.5,
													 cat.cex=c(1.3,1.3),
													 col="transparent",
													 cat.pos=c(340,20),
													 cat.dist=c(.001,.001),
													 cat.col=c("red","blue"),
													 fill=c("red","blue"),
													 main="Overlap of significantly aberrational and drug-targeted paths")
	}else{
		if(VERBOSE) print("venn option 2")
		vennList[["Paths containing drug-sensitive targets"]] = senspaths
		vennTmp = venn.diagram(vennList, filename=NULL,
													 main.cex=1.5,
													 cex=1.5,
													 cat.cex=c(1.3,1.3,1.3)[1:length(vennList)],
													 col="transparent",
													 cat.pos=c(340,20,0)[1:length(vennList)],
													 cat.dist=c(.001,.001,.001)[1:length(vennList)],
													 cat.col=c("red","darkorchid4","blue")[1:length(vennList)],
													 fill=c("red","darkorchid4","blue")[1:length(vennList)],
													 main="Overlap of significantly aberrational,\ndrug-targeted and drug sensitive paths")
	}
	
	png(fullvennFileName)
	grid.draw(vennTmp)
	dev.off()
	if(VERBOSE) print(fullvennFileName)
	
	if(VERBOSE) print(vennFileName)
	
	if(VERBOSE) print("Venn diagram complete.")
	# 	#this is the base name (no folder root) of the venn diagram; 
	# 	basename = strsplit(vennFileName,split="\\/")[[1]][length(strsplit(vennFileName,split="\\/")[[1]])]
	# 	outlist[[basename]] = basename #setting up a base name with a .png should compel toHTML to insert an image. .. lets see if it works
	# 	
	if(VERBOSE) print("Creating overlap tables...")
	#overlap tables	
	######1
	abNotDrug = setdiff(aben, drugcoverage)#aberration enriched, not drug targeted
	if(VERBOSE) print("tabout_abNotDrug")
	if(VERBOSE) print(length(abNotDrug))
	if(VERBOSE) print(dim(combinedAberrations$pathsummary))
	if(VERBOSE) print(is.null(combinedAberrations$pathsummary))
	tabout_abNotDrug = combinedAberrations$pathsummary[abNotDrug,,drop=F]
	if(VERBOSE) print("got combinedAberrations$pathsummary[abNotDrug,,drop=F]")
	outlist[["Aberration enriched, not drug targeted"]] = tabout_abNotDrug
	if(VERBOSE) print("set tabout_abNotDrug")
	if(length(tabout_abNotDrug)>0){
		outlist[["Pathway overlaps of genes in aberration enriched, not drug targeted paths"]]=BangForBuck(path_detail=ola$path_detail,darkPaths=tabout_abNotDrug)
	}else{
		outlist[["Pathway overlaps of genes in aberration enriched, not drug targeted paths"]]="No dark pathways were found (ie, no pathways were found to be aberrationally enriched but not drug targeted)"
	}
	if(VERBOSE){
		print("targAndAb")
		######2
		print(is.null(aben))
		print(is.null(drugcoverage))
	}
	
	targAndAb = intersect(drugcoverage,aben)#Aberrationally enriched, containing drug targets
	if(VERBOSE) print("tabout_targAndAb")
	if(VERBOSE) print(is.null(targAndAb))
	tabout_targAndAb = NULL
	if(!is.null(targAndAb)){
		tabout_targAndAb = merge(coverage$pathsummary[targAndAb,],
														 combinedAberrations$pathsummary[targAndAb,],
														 by="path_id",
														 all=F, 
														 sort=F)
		if(VERBOSE) print("Checking targeted and aberrational overlap..")
		if(nrow(tabout_targAndAb)){
			if(VERBOSE) print("Targeted and aberrational found..")
			plotfname = "aberrational_and_targeted.png"
			fullplotfname = paste(ola$outfname_path,plotfname, sep="")
			#show a bubble plot of the pathways that are targeted and aberrational
			plotOverlapArea(tabout_targAndAb=tabout_targAndAb,
											results=ola$results, 
											abandtargplotname=fullplotfname)
			outlist$imageSlots[[plotfname]] = plotfname
			tabout_targAndAb=cbind(1:nrow(tabout_targAndAb),tabout_targAndAb)
			colnames(tabout_targAndAb)[1]<-"path number"
			if(VERBOSE) print("Plot created.")
		}
	}
	outlist[["Aberrationally enriched, containing drug targets"]] = tabout_targAndAb
	
	if(VERBOSE) print("drugNotAb")
	
	######3
	drugNotAb = setdiff(drugcoverage, aben)#drug targeted, not aberrationally enriched
	if(VERBOSE) print("drugNotAb; merging to make tabout_drugNotAb")
	if(VERBOSE) print(dim(coverage$pathsummary[drugNotAb,]))
	if(VERBOSE) print(dim(combinedAberrations$pathsummary))
	if(is.null(dim(combinedAberrations$pathsummary))){
		if(VERBOSE) print("it's null")
		if(VERBOSE) print(combinedAberrations$pathsummary)
		tabout_drugNotAb = "Patient was not found to have aberrations in currently employed set of pathways."
		if(VERBOSE) print("Patient was not found to have aberrations in currently employed set of pathways.")
	}else{
		tabout_drugNotAb = merge(x=coverage$pathsummary[drugNotAb,],
														 y=combinedAberrations$pathsummary,
														 by="path_id",
														 all.x=T,
														 sort=F)#merges so that all rows with drug screen summary will be included
		if(VERBOSE) print(names(combinedAberrations))
		
		if(VERBOSE) print("merged; zeroing the na's")
		tabout_drugNotAb[is.na(tabout_drugNotAb[,11]),11]<- 0
		
		hyperg_p_value_indexes = grep(pattern="hyperg_p_value",x=colnames(tabout_drugNotAb),ignore.case=T)
		hyperg_p_value_index = hyperg_p_value_indexes[length(hyperg_p_value_indexes)]
		#sort targeted, not aberrationally enriched by the proportion aberrational
		tabout_drugNotAb = tabout_drugNotAb[order(tabout_drugNotAb[,hyperg_p_value_index],decreasing=F),]
	}
	outlist[["Drug targeted, not aberrationally enriched"]] = tabout_drugNotAb
	
	allSensPathsMergAb = "No pathways were found to contain drug-sensitive genes"
	
	if(!nullOrEmpty(functionalEnrichmentAnalysis$pathsummary)){
		allSensPathsMergAb = functionalEnrichmentAnalysis$pathsummary
		if(VERBOSE) print("allSensPathsMergAb; merging to make allSensPathsMergAb")
		if(is.null(dim(combinedAberrations$pathsummary))){
			print("Patient was not found to have aberrations in currently employed set of pathways.")
		}else{
			tabout_allSensPathsMergAb = merge(x=functionalEnrichmentAnalysis$pathsummary,
																				y=combinedAberrations$pathsummary,
																				by="path_id",
																				all.x=T,
																				sort=F)#merges so that all rows with drug screen summary will be included
			if(VERBOSE) print("merged; zeroing the na's")
			tabout_allSensPathsMergAb[is.na(tabout_allSensPathsMergAb[,11]),11]<- 0
			
			hyperg_p_value_indexes = grep(pattern="hyperg_p_value",x=colnames(tabout_allSensPathsMergAb),ignore.case=T)
			hyperg_p_value_index = hyperg_p_value_indexes[length(hyperg_p_value_indexes)]
			#sort targeted, not aberrationally enriched by the proportion aberrational
			tabout_allSensPathsMergAb = tabout_allSensPathsMergAb[order(tabout_allSensPathsMergAb[,hyperg_p_value_index],decreasing=F),]
			allSensPathsMergAb = tabout_allSensPathsMergAb
		}
	}
	outlist[["Paths containing drug-sensitive genes"]] = allSensPathsMergAb
	
	if(VERBOSE) print("sensenAndAb")
	if(is.null(aben)){
		outlist[["Enriched for aberration and enriched for sensitive drug targets"]] = "Patient was not found to have aberrations in currently employed set of pathways."
		outlist[["Aberration enriched, containing sensitive targets"]] = "Patient was not found to have aberrations in currently employed set of pathways."
	}else{
		sensenAndAb = intersect(sensen,aben)#enriched for sensitive targets and enriched for aberration
		outlist[["Enriched for aberration and enriched for sensitive drug targets"]] = as.matrix(sensenAndAb,ncol=1)
		if(VERBOSE) print("conSensAndAb")
		conSensAndAb = intersect(senspaths,aben)#paths containing sensitive genes, that are also aberration enriched
		outlist[["Aberration enriched, containing sensitive targets"]] = as.matrix(conSensAndAb, ncol=1)
	}
	
	# 	cat("\nWriting set of overlap-summary text files to:", paste("overlap_summary", as.character(Sys.time())), "\n")
	# 	outHTMLname = paste(outfname, ".html",sep="")
	#cat("\nWriting html-based summary to ", outHTMLname,"\n")
	ret = list()
	
	outlist[["combined_aberrations_summary"]] = combinedAberrations
	outlist[["functional_drug_screen_summary"]] = functionalEnrichmentAnalysis
	if(VERBOSE){
		print(names(combinedAberrations))
		print(names(combinedAberrations$path_summary_each_patient))
	}
	return(outlist)
}#coreOverlapAnalysis


#nullOrEmpty
#checks if an object, dat, is either null or empty, with zero rows or zero length
nullOrEmpty<-function(dat){
	
	if(is.null(dat)) return(TRUE)
	if(class(dat)=="matrix"|class(dat)=="data.frame"){
		if(!nrow(dat)) return(TRUE)
	}else{
		if(!length(dat)) return(TRUE)
	}
	return(FALSE)
}

# olatest = abDrugOverlapAnalysis(results=results, path_detail=path_detail)
