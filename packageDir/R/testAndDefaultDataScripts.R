#'@title checkFileCopyDefault
#'@description Checks if a file exists. If it doesn't, function checks for a default in the package directory's extdata/ folder and copies the default to the location described by the file name passed.
#'@param fname The path of the file to be checked. 
#'@return string, the file name.
checkFileCopyDefault<-function(fname){
	if(file.exists(fname)){
		return(fname)
	}else{
		defname = system.file(paste0("extdata/",basename(fname)), package = "packageDir")
		if(file.exists(defname)){
			message(paste("Moving default file,", defname, "\nto current working directory at ",fname))
			dir.create(path=dirname(fname), showWarnings=F, recursive=T)
			file.copy(from=defname, to=fname)
			return(fname)
		}else{
			warning(paste("A default file for",fname,"could not be found!!"))
			return(fname)
		}
	}#"extdata/uniqueColorPalette.txt"
}


#'@title Get an example patient gene matrix
#'@description Retreives from package data a patient gene matrix with corresponds with the Reactome pathway, 'Abacavir metabolism' being completely enriched in active (mutated, drug sensitive, etc.) genes.
#'@return \code{matrix} object with column names set as patient IDs, row names set as HUGO gene symbols, and \code{logical} cell values indicating if the corresponding gene is active in the corresponding patient.
#'@export
getTestPGM<-function(){
	#the pgm used is the abacavir metabolism pgm from the moc somatic mutation data
	#pgm=STUDY@results$somatic_mutation_aberration_summary$patientGeneMatrix
	#save(pgm, file="./testData/pgm.rda")
	defname = system.file("testData/pgm.rda", package = "packageDir")
	# 	load(file="./testData/pgm.rda", verbose=T)
	lres = load(file=defname, verbose=TRUE)
	pgmdat = get(lres[1])
	return(pgmdat)
}



#'@title Get a Study object containing moc data. 
#'@description Returns a \code{Study} object loaded with moc somatic mutation data which indicates the Reactome pathway, Abacavir metabolism, is enriched in somatic mutation. 
#'@param noSettings \code{logical} flag indicating if analysis setting should not be automatically added to \code{Study} object.
#'@param noResults \code{logical} flag indicating if analysis results should not be automatically added to \code{Study} object.
#'@param noPaths \code{logical} flag indicating if cellular pathway repository should not be automatically added to \code{Study} object.
#'@param noArms \code{logical} flag indicating if data processing arms should not be loaded with \code{Study} object.
#'@return \code{Study} object
#'@export
getTestStudyObject<-function(noSettings=FALSE, noResults=FALSE, noPaths=FALSE, noArms=FALSE){
	# 	exampleResults = STUDY@results
	# 	save(exampleResults, file="../packageDir/inst/testData/exampleResults.rda")
	# 	exampleSettings = STUDY@studyMetaData@settings
	# 	save(exampleSettings, file="../packageDir/inst/testData/exampleSettings.rda")
	
	exampleSettings=list()

	if(!noSettings){
		l2 = load(file=system.file("testData/exampleSettings.rda", package = "packageDir"), verbose=T)
		exampleSettings = get(l2[1])
	}

	exampleStudy = getStudyObject(study.name="testDataSets", 
																geneIdentifierType="HUGO",
																path_detail=getDefaultPaths(),
																settings=exampleSettings)
	
	if(!noArms) exampleStudy = loadBasicArms(STUDY=exampleStudy)
  
	if(!noResults){
		l1 = load(file=system.file("testData/exampleResults.rda", package = "packageDir"), verbose=T)
		exampleResults = get(l1[1])
		exampleStudy@results = exampleResults
	}

	return(exampleStudy)
}


#'@title Get settings for the three basic study arms: functional, aberration and function-aberration overlap.
#'@description Returns the settings list for the analysis of the test data sets: functional drug screen data, somatic mutation data and overlap analysis, of the functional and somatic mutation (aberration) data.
#'@return \code{list} object with slots named "overlap_analysis", "somatic_mutation_aberration_summary" and "functional_drug_screen_summary" containing settings lists for their respective study arms. 
#'@export
getTestStudySettings<-function(){
	# 	testStudySettings = STUDY@studyMetaData@settings
	# 	save(testStudySettings, file="../packageDir/inst/testData/testStudySettings.rda")
	lres = load(file="../packageDir/inst/testData/testStudySettings.rda", verbose=T)
	dols = get(lres[1])
	return(dols)
}

getTestPaths<-function(path_file=system.file("extdata/Reactome.2014.04.06.12.52.27.txt", package = "packageDir"), 
											 forceReload=F){
	
	if(!exists("testEnv")){
		testEnv<<-new.env(parent=globalenv())
	}
	
	if( !forceReload & exists("testEnv") ){
		if(exists("path_detail",envir=testEnv)){
			return(testEnv$path_detail)
		}
	}
	
	testEnv$path_detail<-getPaths(path_file=path_file, verbose=F)
	return(testEnv$path_detail)
}


#'@title Returns an example of the PathSummaryRunner reference class object used by 
#'@description Returns an object of the PathSummaryRunner reference class used by summaryTable. 
#'@param pths A \code{Path_Detail} object. If left null, (the default), the default reactome pathways will be automatically provided. 
#'@param pathIds A \code{character} vector containing the name(s) of paths to be artificially enriched.
#'@return An instance of the \code{PathSummaryRunner} class, loaded with data indicating enrichment of the "Abacavir metabolism" Reactome pathway.
#'@export
getTestPathSummaryRunner<-function(pths=NULL, pathIds="Abacavir metabolism"){
	
	if(is.null(pths)) pths = getTestPaths()
	
	geneSet = getGenesFromPaths(pids=pathIds, STUDY=pths)
	
	tm = getTargetMatrix(tgenes=geneSet, paths=pths$paths)
	
	psr = PathSummaryRunner$new(.verbose=F, path_summary_each_patient=list())
	psr$targetname = "testTarget"
	psr$dataSetName = "unit test, test data set"
	psr$.targetMatrix = tm
	psr$patientGeneMatrix = getTestPGM()
	
	return(psr)
}


# psepFile = "./output/pathSummaryRunnerObject_psep.rda"
# save(psep,file=psepFile)
# psrFile= "./output/pathSummaryRunnerObject.rda"
# save(psr,file=psrFile)

#'@title Get \code{path_summary_runner} object. 
#'@description The  \code{path_summary_runner} object returned by this function is one of the two passed by summaryTable to the path significance tests.
#'@return \code{path_summary_runner} object loaded with the data that is provided for path significance test functions.
#'@seealso \code{\link{summaryTable}}
getLoadedPathSummaryRunner<-function(){
	r1 = load(system.file("extdata/pathSummaryRunnerObject.rda", package="packageDir"), verbose=T)
	psr = get(r1[1])
	return(psr)
}



