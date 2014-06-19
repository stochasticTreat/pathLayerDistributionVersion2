#unit tests for drugDbInterface.R


test.getDrugTargetData<-function(){
	require(HGNChelper)
	
	checkTrue(expr=is.null(packageDir:::getDrugTargetData(hugo=hgnc.table, fname="notAFileName")) )

	tmpFileName = "./unitTestTmpFile.txt"
	
	#targtab = read.table(file=fname, sep=sep, header=header, quote=quote, stringsAsFactors=F)[1:10,]
	#save(targtab, file="./inst/testData/getDrugTargetDataTestData.rda")
	varname = load(system.file("extdata/abacavirSettings.rda",package="packageDir"), verbose=T)
	
	load(file=system.file("extdata/getDrugTargetDataTestData.rda", package="packageDir"), verbose=T)
	
	write.table(x=targtab,file=tmpFileName, sep=",")
	resNew = packageDir:::getDrugTargetData(hugo=hgnc.table, fname=tmpFileName)
	#save(res1, file="inst/testData/getDrugTargetDataTestDataTrueOutput.rda")
	load(file=system.file("extdata/getDrugTargetDataTestDataTrueOutput.rda", package="packageDir"), verbose=T)
	checkEquals(target=res1, current=resNew)
	
	file.remove(tmpFileName)
	checkEquals(target=T, current=T)
}