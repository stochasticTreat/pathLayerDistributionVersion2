#4.R: unit tests for acc_functions.R
#test.addPidColumn()
test.addPidColumn<-function(){
	
	dfin = data.frame(Tumor_Sample_Barcode=c("TCGA-AB-2988-03A-01D-0739-09", "TCGA-AB-2989-03A-01D-0739-09"), 
										Matched_Sample_Barcode=c("TCGA-AB-2988-11A-01D-0739-09", "TCGA-AB-2989-12A-01D-0739-09"),
										Hugo_Symbol=c("TP53","FLT3"))
	dfout = packageDir:::addPidColumn(tcga_data=dfin)
	checkEquals(target=c("TCGA.AB.2988","TCGA.AB.2989"), current=as.character(dfout$pid))
	
}

#test.emptySignature()
test.emptySignature<-function(){
	
	tdf = data.frame(matrix(data="V1", nrow=1, ncol=1, dimnames=list(NULL,"V1")),stringsAsFactors=F)
	
	checkTrue(expr=packageDir:::emptySignature(testdf=tdf))
	
	tdf[1,1] = "testPid"
	checkTrue(expr=packageDir:::emptySignature(testdf=tdf, pid="testPid"))
	
	tdf[1,1] = 12
	checkTrue(expr=!packageDir:::emptySignature(testdf=tdf))
	
}