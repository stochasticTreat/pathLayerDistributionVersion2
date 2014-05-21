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


