

#test.hypergeometricPathEnrichment
test.hypergeometricPathEnrichment<-function(){
	
	pths = getDefaultPaths()
	
	psr = getTestPathSummaryRunner(pths=pths)
	
	hpres1 = hypergeometricPathEnrichment(pathSigRunner=psr, 
																					 paths_detail=pths) 
	
# 	save(hpres, file="./testData/hypergeometricPathEnrichmentRes.rda")
	load(file=system.file("testData/hypergeometricPathEnrichmentRes.rda", package = "packageDir"), verbose=T)
	
	checkEquals(target=hpres, current=hpres1)
	
}

#test.getActiveGenesEachPath()
test.getActiveGenesEachPath<-function(){
	
	pths = getDefaultPaths()
	psr = getTestPathSummaryRunner(pths=pths)
	
	agep1 = packageDir:::getActiveGenesEachPath(psr=psr)
	agepTestDataFile = system.file("testData/abacavirActiveGenesEachPath.rda", package = "packageDir")
	#save(agep, file=agepTestDataFile) #created on 4/16/14 with data from 
	load(file=agepTestDataFile, verbose=T)
	
	checkEquals(target=agep, current=agep1)
	
}

#test.getBasicPathInformation()
test.getBasicPathInformation<-function(){
	
	pths = getDefaultPaths()#path_file="./reference_data/paths/Reactome 2014.04.06 12.52.27.txt",force=T)
	
	psr = getTestPathSummaryRunner(pths=pths)
	
	psr$study = getStudyObject(study.name="testerStudy", 
														 path_detail=pths,
														 settings=NULL)
	
	bpi_current = packageDir:::getBasicPathInformation(paths_detail=pths, pathNames=rownames(psr$.targetMatrix), psr=psr)
	tdfn = system.file("testData/basicPathInfoAbacavir.rda", package = "packageDir")
	#save(bpi, file=tdfn)
	load(file=tdfn, verbose=T)
	checkEquals(target=bpi, current=bpi_current)
	
# 	CheckAll(old=bpi, novo=bpi_current)
	print("Trying with reduced coverage..")
	remset = c("ADAL",  "ADH1A", "GUK1")
	psr$coverage = setdiff(x=colnames(pths$paths), remset)
	ccres = packageDir:::checkCoverage(psr=psr, paths_detail=pths, geneDescription="covered", "coverage analysis")
	psr = ccres$psr
	pths = ccres$paths_detail
	geneSet = setdiff(packageDir:::getGenesFromPaths(pids="Abacavir metabolism", STUDY=pths), 
										remset)
	psr$.targetMatrix = packageDir:::getTargetMatrix(tgenes=geneSet, paths=pths$paths)
	
	bpiCov_current = packageDir:::getBasicPathInformation(paths_detail=pths, 
																					 psr=psr,
																					 pathNames=rownames(psr$.targetMatrix))
	
	checkTrue(expr=sum(bpiCov_current$full_path_length != bpiCov_current$testable_path_length)>0)
	tfn2 = system.file("testData/basicPathInfoAbacavirCoverage.rda", package = "packageDir")
	# 	save(bpiCov, file=tfn2)
	load(file=tfn2, verbose=T)
	
	checkEquals(target=bpiCov, current=bpiCov_current)
	
}

#test.generalSummary()
test.generalSummary<-function(){
	
	pths = packageDir:::getDefaultPaths()
	psr = packageDir:::getTestPathSummaryRunner(pths=pths)
	
	
	psr$patientsum = t((rep(x=1,times=nrow(psr$patientGeneMatrix))%*%psr$patientGeneMatrix)) #patientsum: matrix with the number of active genes found in each patient
	colnames(	psr$patientsum )<-"count_per_patient"
	
	gs1 = packageDir:::generalSummary(psr=psr, paths_detail=pths)
	#gs = gs1
	#save(gs, file="../packageDir/inst/testData/generalSummary.rda")
	load(system.file("testData/generalSummary.rda", package = "packageDir"), verbose=T)
	
	checkEquals(target=gs, current=gs1)
	
}

test.loadSettings<-function(){
	fname = "./testSaveSettings.txt"
	file.remove(fname)
	study = getTestStudyObject()
	originalSettings = study@studyMetaData@settings$defaultSummaryTable
	dflist = packageDir:::saveSettings(set=originalSettings)
	write.table(x=dflist, file=fname, sep="\t")
	reopened = packageDir:::loadSettings(fname=fname)
	checkEquals(target=originalSettings, current=reopened)
	file.remove(fname)
}
