library('RUnit')

getTestPaths<-function(path_file=system.file("extdata/Reactome.2014.04.06.12.52.27.txt", package = "packageDir"), 
											 forceReload=F){
	
	
	if( !forceReload & exists("testEnv") ){
		if(exists("path_detail",envir=testEnv)){
			return(testEnv$path_detail)
		}
	}
	
	testEnv<<-new.env(parent=globalenv())
	testEnv$path_detail<-getPaths(path_file=path_file)
	return(testEnv$path_detail)
	
}

getTestPathSummaryRunner<-function(pths=NULL){
	
	if(is.null(pths)) pths = getTestPaths()
	
	geneSet = getGenesFromPaths(pids="Abacavir metabolism", STUDY=pths)
	
	tm = getTargetMatrix(tgenes=geneSet, paths=pths$paths)
	
	psr = PathSummaryRunner$new(.verbose=F, path_summary_each_patient=list())
	psr$targetname = "testTarget"
	psr$dataSetName = "unit test, test data set"
	psr$.targetMatrix = tm
	psr$patientGeneMatrix = getTestPGM()
	
	return(psr)
}

getTestPGM<-function(){
	#the pgm used is the abacavir metabolism pgm from the moc somatic mutation data
	#pgm=STUDY@results$somatic_mutation_aberration_summary$patientGeneMatrix
	#save(pgm, file="./testData/pgm.rda")
	defname = system.file("testData/pgm.rda", package = "packageDir")
# 	load(file="./testData/pgm.rda", verbose=T)
	return(pgm)
}

message("Running Tests")

test.suite <- defineTestSuite("example",
															dirs = system.file("unitTests", package = "packageDir"),
															testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)