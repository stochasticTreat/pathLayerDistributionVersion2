#integration tests

makeRandomPgm<-function(npat=5, probabilities=c(1,4), targetPath="ABCA transporters in lipid homeostasis"){
	
	pths = getDefaultPaths()
	genes = packageDir:::getGenesFromPaths(pids=targetPath, STUDY=pths)
	
	patients = paste("p", 1:npat, sep="")
	
	ppgm = matrix(data=sample(x=c(T,F),
														prob=probabilities,
														size=(length(genes)*length(patients)), 
														replace=T), 
								nrow=length(genes), 
								ncol=length(patients), 
								dimnames=list(genes, patients))
	return(ppgm)
}

# test.SummaryTable<-function(){
# 	
# 	
# 	aberration_data_type="generic data type"
# 	dataSetName = "generic data set name"
# 	targetName = "generic target name"
# 	pgmFile = "/Users/samhiggins2001_worldperks/tprog/main_130525/input/test_files/overlaptestGisticFile.txt"
# 	pgm = open.PGM(fname=pgmFile)
# 	pgm = pgm==0
# 	target_matrix = NULL
# 	path_detail  = STUDY@studyMetaData@paths
# 	
# 	dir.create("./testData/")
# 	path_detail_file_name= "./testData/path_detail_object"
# 	pgm_file_name = "./testData/pgm"
# 	save(path_detail, file=path_detail_file_name, compress=T)
# 	save(pgm, file=pgm_file_name, compress=T)
# 	
# 	
# 	gsum = summaryTable4(paths_detail=path_detail, 
# 											 target_matrix=target_matrix,
# 											 individualEnrichment=T,
# 											 verbose=T,
# 											 dataSetName=dataSetName,
# 											 patientGeneMatrix=pgm, 
# 											 targetname=targetName)
# 	
# }

#test.summaryTable()
test.summaryTable<-function(){
	
	testStudyObjectFileNamesomatic = system.file("testData/amlSomaticSum.rda", package = "packageDir")
	
	if(!file.exists("./testData")) dir.create("./testData", recursive=T, showWarnings=F)
	#amlSom = STUDY@results$somatic_mutation_aberration_summary
	#save(amlSom, file=testStudyObjectFileNamesomatic, compress=T)
	load(file=testStudyObjectFileNamesomatic, verbose=T)
	
	stemp = getStudyObject()
	stemp@studyMetaData@paths = getDefaultPaths()
	
	stResSomatic = summaryTable(study=stemp,
															settings=amlSom$settings,
															individualEnrichment=T,
															pgm=amlSom$patientGeneMatrix,
															activeGeneDescription="mutated", 
															dataSetDescription="somatic mutation data")
	
	packageDir:::saveSummary(summ=stResSomatic, study_name="tmpForUnitTest", path="./testData")
	stResSomatic_reopened = loadSummary(study_name="tmpForUnitTest", path="./testData")
	
	#extract the slots that should be the same
	commonSlots = c('pathsummary', 'summarystats', 'patientGeneMatrix', 'patientsums', 
									'genesummary', 'active_genes_ea_path', 'path_summary_each_patient', 
									'active_genes_not_in_paths', 'patientList', 'settings')
	
	amlSom = amlSom[commonSlots]
	stResSomatic_reopened = stResSomatic_reopened[commonSlots]
	
	checkEquals(target=amlSom, 
							current=stResSomatic_reopened, 
							msg="Checking whole summary Table result set")
	
	unlink(x="./testData/tmpForUnitTest", recursive=T, force=T)
	#CheckAll(old=amlSom, novo=stResSomatic_reopened, verbose=T)
}


#test.summaryTable()
test.coverage.summaryTable<-function(){
	
	testStudyObjectFileNamesomatic = system.file("testData/amlDrugScreenSum.rda", package = "packageDir")
	if(!file.exists("./testData")) dir.create("./testData", recursive=T, showWarnings=F)
	
	#amlFun = stResFunctional
	#save(amlFun, file="../packageDir/inst/testData/amlDrugScreenSum.rda", compress=T)
	load(file=testStudyObjectFileNamesomatic, verbose=T)
	
	stemp = getStudyObject()
	stemp@studyMetaData@paths = getDefaultPaths()
	
	stResFunctional = summaryTable(study=stemp, 
															coverageDataSetDescription="Drug screen coverage",
															coverageGeneDescription="drug_targeted",
															coverage=rownames(amlFun$coverage_summary$genesummary),
															settings=amlFun$settings,
															individualEnrichment=T,
															pgm=amlFun$patientGeneMatrix,
															activeGeneDescription="drug_sensitive", 
															dataSetDescription="Drug screen output data")
	
	packageDir:::saveSummary(summ=stResFunctional, study_name="tmpForUnitTest", path="./testData")
	stResFunctional_reopened = loadSummary(study_name="tmpForUnitTest", path="./testData")

	#extract the slots that should be the same
	
	commonSlots = intersect(names(amlFun), names(stResFunctional_reopened))
	
	# 	commonSlots = c('pathsummary', 'summarystats', 'patientGeneMatrix', 'patientsums', 
	# 									'genesummary', 'active_genes_ea_path', 'path_summary_each_patient', 
	# 									'active_genes_not_in_paths', 'patientList', 'settings')
	
	amlFun = amlFun[commonSlots]
	stResFunctional_reopened = stResFunctional_reopened[commonSlots]
	badslots=c()
	for(sn in names(amlFun)){
		print(sn)
		tres = try(checkEquals(target=amlFun[[sn]], 
											current=stResFunctional_reopened[[sn]], 
											msg="Checking whole summary Table result set"),silent=TRUE)
		if(class(tres)=="try-error"){
			print(tres)	
			badslots = c(badslots, sn)
		} 
	}
	if(length(badslots)){
		checkEquals(target=T,current=F, msg=paste("These are the bad slots: ",paste(badslots,sep=" ",collapse="")))
	}
# 	covT = amlFun$path_summary_each_patient$AMLadult07335
# 	covC = stResFunctional_reopened$path_summary_each_patient$AMLadult07335
# 	badslots=c()
# 	for(sn in names(covT)){
# 		tres = try(checkEquals(target=covT[[sn]], 
# 													 current=covC[[sn]], 
# 													 msg="Checking whole summary Table result set"),silent=TRUE)
# 		if(class(tres)=="try-error"){
# 			print(tres)	
# 			badslots = c(badslots, sn)
# 		} 
# 	}

	unlink(x="./testData/tmpForUnitTest/", recursive=T, force=T)
}

test.matchPath.SummayTable<-function(){
	
	stemp = getStudyObject()
	stemp@studyMetaData@paths = getDefaultPaths()
	
	targetPath="ABCA transporters in lipid homeostasis"
	pgm1=makeRandomPgm(npat=6, probabilities=c(1,4),targetPath=targetPath)
	
	stRes = summaryTable(study=stemp, 
											 individualEnrichment=T,
											 pgm=pgm1,
											 activeGeneDescription="testerGenes", 
											 dataSetDescription="test data set")
	
	checkTrue(expr=targetPath%in%rownames(stRes$pathsummary), msg="Checking randomized path analysis for target path inclusion")
	
}