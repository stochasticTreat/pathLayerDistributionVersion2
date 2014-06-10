#3.R

test.fail<-function(){
	
	stop("Checking that the unit tests work and will report a failure")
	
}

test.UserProvidedBiopax<-function(){
	
	#clean everything up and get it all ready: delete everything in the biopax directory
	pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	for(fn in dir( dirname(path=pwrecord.fileName) )){ file.remove( paste0(dirname(path=pwrecord.fileName), "/", fn) ) }
	
	#user provided biopax tests:
	#1) just supply path names and require the user select the files to go with
	checkEquals(msg="Checking copy of single path via user input",
							target=character(0), current=UserProvidedBiopax(pathNames="Abacavir metabolism"))
	
	#clean it all up again
	for(fn in dir( dirname(path=pwrecord.fileName) )){ file.remove( paste0(dirname(path=pwrecord.fileName), "/", fn) ) }
	
	#check if the pathway file name is provided manually
	checkEquals(target=character(0), 
							current=UserProvidedBiopax(pathNames="Abacavir metabolism", 
																				 pathMetaDataFile="~/tprog/main_131219/reference_data/paths/biopax/record_of_biopax_pathways.txt"))
	pwrecord = getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	
	abacavirFname = paste0(pwrecord$dbID[pwrecord$path_name=="Abacavir metabolism"],".owl")
	allFiles = dir(dirname(pwrecord.fileName))
	allFiles = allFiles[grepl(pattern=".owl", x=allFiles)]
	
	checkTrue(expr=(abacavirFname%in%allFiles), 
						msg="checking that the biopax file is in the list and in the correct folder")
	
	#check that there's one row in the records file for each biopax file
	checkTrue(expr=(sum(paste0(pwrecord$dbID,".owl")%in%allFiles)==nrow(pwrecord)))
	
	#	 	write.table(x=pwrecord, file=pwrecord.fileName, quote=F,sep="\t", row.names=F, col.names=T)
}


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


test.getDefaultPaths<-function(){
	cellular_pathways = getDefaultPaths()
	checkEquals(target="Path_Detail", current=class(cellular_pathways)[1])
	checkEquals(target="Reactome", current=cellular_pathways$source)
	checkTrue(expr=(nrow(cellular_pathways$paths)>1 & ncol(cellular_pathways$paths)>1))
	checkEquals(target="matrix", current=class(cellular_pathways$paths))
}


.test.addBiopaxPath<-function(){
	
	addBiopaxPath(pname="Abacavir metabolism", 
								dbid="p53-Dependent G1 DNA Damage Response.owl", 
								biopaxDat="/Users/samhiggins2001_worldperks/tprog/main_131219/reference_data/paths/biopax/p53-Dependent G1 DNA Damage Response.owl")
	pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	pwrecord = getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	
	
	checkTrue(file.exists("./reference_data/paths/biopax/p53-Dependent G1 DNA Damage Response.owl"))
	checkTrue("p53-Dependent G1 DNA Damage Response"%in%pwrecord$dbID)
	
	pwrecord = pwrecord[pwrecord$dbID!="p53-Dependent G1 DNA Damage Response",]
	write.table(x=pwrecord, file=pwrecord.fileName, quote=F,sep="\t", row.names=F, col.names=T)
	
	file.remove("./reference_data/paths/biopax/p53-Dependent G1 DNA Damage Response.owl")
	
}



test.biopaxFileNameFromPathName<-function(){
	library(RUnit)
	#unit test test.biopaxFileNameFromPathName()
	print("test.biopaxFileNameFromPathName()")
	#positive control
	tpnames = c("Abacavir metabolism", 
							"Transcriptional Regulation of White Adipocyte Differentiation")
	bpnames = biopaxFileNameFromPathName(pathNames=tpnames)
	checkEquals(target=2, current=length(bpnames))
	
	#negative control
	tpnames2 = c("Abacavir metabolism", 
							 "not a pathway")
	bpnames2 = biopaxFileNameFromPathName(pathNames=tpnames2)
	
	checkEquals(target=1, current=length(bpnames2))
	
}


test.checkCurrentBiopax<-function(){
	biopax.dir = "./reference_data/paths/"
	pathNames = c("Abacavir metabolism", "Not a pathway name")
	# 	pathNames="Abacavir metabolism"
	checkEquals(target="Not a pathway name", 
							current=checkCurrentBiopax(pathNames=pathNames))
}


test.getReactomeBiopax<-function(){
	# 	tpnames = c("Abacavir metabolism", "Transcriptional Regulation of White Adipocyte Differentiation")
	tpnames = c("Abacavir metabolism")
	getReactomeBiopax(study=NULL, pathNames=tpnames)
}


test.correctByHgncHelper<-function(){
	
	# 	symbol_set = rownames(STUDY@results$somatic_mutation_aberration_summary$genesummary)
	symbol_set = c("p53", "flt3", "oct5", "BIKE", "LOK")
	res = correctByHgncHelper(symbol_set=symbol_set, verbose=F)
	
	checkEquals(target=c("TP53", "FLT3", "OCT5", "BMP2K", "STK10"), current=res)
	
}