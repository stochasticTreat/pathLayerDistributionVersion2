
filterMafToMutsig<-function(){
	
	readline("press any key to continue and select a .maf file to be filtered..")
	fname = file.choose()
	readline("press any key to continue and select a MutSig results file to use for filtering..")
	fname2 = file.choose()
	mafin = read.delim(file=fname, header=T, sep="\t", stringsAsFactors=F, na.strings="-")
	msigGeneTab = read.delim(file=fname2, header=T, sep="\t", stringsAsFactors=FALSE, comment.char="")

	msigGenes = msigGeneTab[msigGeneTab$q<.1,"gene",drop=TRUE]
	
	mafin$Hugo_Symbol = corsym(symbol_set=mafin$Hugo_Symbol, symref=STUDY)
	msigGenes = corsym(symbol_set=msigGenes, symref=STUDY)
	
	msigmaf = mafin[mafin$Hugo_Symbol%in%msigGenes,,drop=FALSE]
	
	fname3 = paste0(fname,".",basename(fname2),".mutsigex.txt" )
	cat("\nSaving filtered .maf\n")
	write.table(x=msigmaf, file=fname3, quote=T, sep="\t", row.names=F, col.names=T)
	cat("\nFiltered .maf saved to file:", fname3, "\nin folder:", getwd())
}


# filterMafToMutsig()
