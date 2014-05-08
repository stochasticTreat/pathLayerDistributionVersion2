#3.R
source('~/tprog/main_131219/main_functions.R')

test.checkForceRowNames<-function(){
	fname = "./testData/test.checkForceRowNames.rda"
	fname2 = "./testData/testTarget.checkForceRowNames.rda"
	#save(incTab,file=fname)
	#save(tab, file=fname2)
	load(fname, verbose=T)
	load(fname2, verbose=T)
	wrn = checkForceRowNames(tab=incTab)
	
	checkEquals(target=tab, current=wrn, msg="regular row names check")
	wrn2 = checkForceRowNames(tab=incTab)
	checkEquals(target=wrn, current=wrn2, msg="row names when there shouldn't be a change")
	
}