#3.R


test.checkForceRowNames<-function(){
	fname = system.file("testData/test.checkForceRowNames.rda", package = "packageDir")
	fname2 = system.file("testData/testTarget.checkForceRowNames.rda", package = "packageDir")
	#save(incTab,file=fname)
	#save(tab, file=fname2)
	load(fname, verbose=T)
	load(fname2, verbose=T)
	wrn = checkForceRowNames(tab=incTab)
	
	checkEquals(target=tab, current=wrn, msg="regular row names check")
	wrn2 = checkForceRowNames(tab=incTab)
	checkEquals(target=wrn, current=wrn2, msg="row names when there shouldn't be a change")
	
}


test.loadPathsAsSets<-function(){
	
	#fname1 = "./reference_data/paths/Reactome.2013.12.27.18.00.18.txt"
	#psets1 = loadPathsAsSets(firstGeneColum=3, fname=fname1)
	
	fname2 = "./reference_data/paths/Reactome.2014.04.06.12.52.27.txt"
	psets2 = loadPathsAsSets(firstGeneColum=3, fname=fname2)
	
	checkEquals(target=1459, current=length(psets2))
	checkEquals(target=4, current=length(psets2[[1]]))
}