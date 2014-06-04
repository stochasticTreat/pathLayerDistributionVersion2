#unit tests for drugDbInterface.R


test.getDrugTargetData<-function(){
	require(HGNChelper)
	
	checkTrue(expr=is.null(getDrugTargetData(hugo=hgnc.table, fname="notAFileName")) )

	tmpFileName = "./unitTestTmpFile.txt"
	
	#targtab = read.table(file=fname, sep=sep, header=header, quote=quote, stringsAsFactors=F)[1:10,]
	#save(targtab, file="./inst/testData/getDrugTargetDataTestData.rda")
	load(file="inst/testData/getDrugTargetDataTestData.rda", verbose=T)
	write.table(x=targtab,file=tmpFileName, sep=",")
	resNew = getDrugTargetData(hugo=hgnc.table, fname=tmpFileName)
	#save(res1, file="inst/testData/getDrugTargetDataTestDataTrueOutput.rda")
	load(file="inst/testData/getDrugTargetDataTestDataTrueOutput.rda", verbose=T)
	checkEquals(target=res1, current=resNew)
	
	file.remove(tmpFileName)
}