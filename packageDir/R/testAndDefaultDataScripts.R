

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
	load(file=defname, verbose=T)
	return(pgm)
}



#'@title Get a Study object containing moc data. 
#'@description Returns a \code{Study} object loaded with moc somatic mutation data which indicates the Reactome pathway, Abacavir metabolism, is enriched in somatic mutation. 
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
