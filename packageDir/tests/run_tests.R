library('RUnit')

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


getTestPathSummaryRunner<-function(pths=NULL){
	
	if(is.null(pths)) pths = getTestPaths()
	
	geneSet = packageDir:::getGenesFromPaths(pids="Abacavir metabolism", STUDY=pths)
	
	tm =  packageDir:::getTargetMatrix(tgenes=geneSet, paths=pths$paths)
	
	psr = PathSummaryRunner$new(.verbose=F, path_summary_each_patient=list())
	psr$targetname = "testTarget"
	psr$dataSetName = "unit test, test data set"
	psr$.targetMatrix = tm
	psr$patientGeneMatrix = getTestPGM()
	
	return(psr)
}

BiocGenerics:::testPackage("packageDir")

message("Running Tests")

test.suite <- defineTestSuite("example",
															dirs = system.file("unitTests", package = "packageDir"),
															testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)