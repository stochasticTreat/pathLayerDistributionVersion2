library('RUnit')

source('~/tprog/main_131219/pathway_functions.R')
source('~/tprog/main_131219/acc_functions.R')
source('~/tprog/main_131219/summaryTable5.R')
source('~/tprog/main_131219/summaryTableFunctions.R')

getTestPaths<-function(path_file="./reference_data/paths/Reactome 2014.04.06 12.52.27.txt", 
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
	load(file="./testData/pgm.rda", verbose=T)
	return(pgm)
}

test.suite <- defineTestSuite("example",
															dirs = file.path("tests"),
															testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)