#4.R: unit tests for acc_functions.R
test.addPidColumn<-function(){
	
	dfin = data.frame(Tumor_Sample_Barcode=c("TCGA-AB-2988-03A-01D-0739-09", "TCGA-AB-2989-03A-01D-0739-09"), 
										Matched_Sample_Barcode=c("TCGA-AB-2988-11A-01D-0739-09", "TCGA-AB-2989-12A-01D-0739-09"),
										Hugo_Symbol=c("TP53","FLT3"))
	dfout = addPidColumn(tcga_data=dfin)
	checkEquals(target=c("TCGA.AB.2988","TCGA.AB.2989"), current=as.character(dfout$pid))
	
}