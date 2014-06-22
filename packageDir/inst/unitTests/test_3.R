#3.R

# test.UserProvidedBiopax<-function(){
# 	
# 	checkEquals(target="Checking that the unit tests work and will report a failure", current="no")
# 	
# }
#test.UserProvidedBiopax()
test.UserProvidedBiopax<-function(userInput=FALSE){
	print("testing UserProvidedBiopax...")
	defaultFile = system.file("extdata/record_of_biopax_pathways.txt",package="packageDir")
	#clean everything up and get it all ready: delete everything in the biopax directory
	pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"#system.file("extdata/record_of_biopax_pathways.txt",package="packageDir") #
	
	print("test 1, assuring abacavir metabolism pathway's name will not be returned if the path is present.. ")
	packageDir:::checkFileCopyDefault(fname=pwrecord.fileName)
	pwrecord = packageDir:::getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	
	biopaxFname = paste0(dirname(pwrecord.fileName),"/",pwrecord$dbID[1],".owl")
	packageDir:::checkFileCopyDefault(fname=biopaxFname)
	
	#check if the pathway file name is provided manually
	checkEquals(target=character(0), 
							current=packageDir:::UserProvidedBiopax(pathNames="Abacavir metabolism", 
																											pathMetaDataFile=pwrecord.fileName ))
	
	print("test 2, seeing if abcavir metabolism's biopax file is there... ")
	biopaxFname = paste0(pwrecord$dbID[pwrecord$path_name=="Abacavir metabolism"],".owl")
	allFiles = dir(dirname(pwrecord.fileName))
	allFiles = allFiles[grepl(pattern=".owl", x=allFiles)]
	checkTrue(expr=(biopaxFname%in%allFiles), 
						msg="checking that the biopax file is in the list and in the correct folder")
	print("checking that there's one row in the records file for each biopax file")
	checkTrue(expr=(sum(paste0(pwrecord$dbID,".owl")%in%allFiles)==nrow(pwrecord)))
	
	print("cleaning out the records for the next tests")
	#clean out the directory
	if(file.exists(dirname(pwrecord.fileName))){
		for(fn in dir( dirname(path=pwrecord.fileName) )){ file.remove( paste0(dirname(path=pwrecord.fileName), "/", fn) ) }
	}else{
		dir.create(dirname(pwrecord.fileName))
	}
	
	print("seeing if a remote biopax path meta data can be pointed to, to move the pathways to the local dir.")
	#now supply the file name as the remote repository
	checkEquals(target="not a path", 
							current=packageDir:::UserProvidedBiopax(pathNames=c("Abacavir metabolism","not a path"), pathMetaDataFile=defaultFile))
	
	#user provided biopax tests:
	#1) just supply path names and require the user select the files to go with
	if(userInput){
		print("Checking procurement of biopax with user input.. ")
			checkEquals(msg="Checking copy of single path via user input",
							target=character(0), current=packageDir:::UserProvidedBiopax(pathNames="Abacavir metabolism"))
			#clean it all up again
			for(fn in dir( dirname(path=pwrecord.fileName) )){ file.remove( paste0(dirname(path=pwrecord.fileName), "/", fn) ) }
	}
	
	print("cleaning up..")
	#clean out the directory
	if(file.exists(dirname(pwrecord.fileName))){
		for(fn in dir( dirname(path=pwrecord.fileName) )){ file.remove( paste0(dirname(path=pwrecord.fileName), "/", fn) ) }
	}
	#	 	write.table(x=pwrecord, file=pwrecord.fileName, quote=F,sep="\t", row.names=F, col.names=T)
}

#test.checkForceRowNames()
test.checkForceRowNames<-function(){
	fname = system.file("testData/test.checkForceRowNames.rda", package = "packageDir")
	fname2 = system.file("testData/testTarget.checkForceRowNames.rda", package = "packageDir")
	#save(incTab,file=fname)
	#save(tab, file=fname2)
	load(fname, verbose=T)
	load(fname2, verbose=T)
	wrn = packageDir:::checkForceRowNames(tab=incTab)
	
	checkEquals(target=tab, current=wrn, msg="regular row names check")
	wrn2 = packageDir:::checkForceRowNames(tab=incTab)
	checkEquals(target=wrn, current=wrn2, msg="row names when there shouldn't be a change")
}

#test.loadPathsAsSets()
test.loadPathsAsSets<-function(){
	
	#fname1 = "./reference_data/paths/Reactome.2013.12.27.18.00.18.txt"
	#psets1 = loadPathsAsSets(firstGeneColum=3, fname=fname1)
	fname2 = system.file("extdata/Reactome.2014.04.06.12.52.27.txt",package="packageDir")
	psets2 = packageDir:::loadPathsAsSets(firstGeneColum=3, fname=fname2)
	
	checkEquals(target=1459, current=length(psets2))
	checkEquals(target=4, current=length(psets2[[1]]))
}

#test.getDefaultPaths()
test.getDefaultPaths<-function(){
	cellular_pathways = getDefaultPaths()
	checkEquals(target="Path_Detail", current=class(cellular_pathways)[1])
	checkEquals(target="Reactome", current=cellular_pathways$source)
	checkTrue(expr=(nrow(cellular_pathways$paths)>1 & ncol(cellular_pathways$paths)>1))
	checkEquals(target="matrix", current=class(cellular_pathways$paths))
}

#test.addBiopaxPath()
test.addBiopaxPath<-function(){
	
	if(file.exists("./reference_data/paths/biopax/2161541.owl")) file.remove("./reference_data/paths/biopax/2161541.owl")
	if(file.exists("./reference_data/paths/biopax/record_of_biopax_pathways.txt")) file.remove("./reference_data/paths/biopax/record_of_biopax_pathways.txt")
	
	
	packageDir:::addBiopaxPath(pname="Abacavir metabolism", 
								dbid="2161541", 
								biopaxDat=system.file("extdata/2161541.owl", package="packageDir"))
	
	pwrecord.fileName = system.file("extdata/record_of_biopax_pathways.txt",package="packageDir")
	#pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	pwrecord = packageDir:::getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	
	checkTrue(file.exists(paste0("./reference_data/paths/biopax/",pwrecord$dbID,".owl")))
	checkTrue("2161541"%in%pwrecord$dbID)
	
	write.table(x=pwrecord, file=pwrecord.fileName, quote=F,sep="\t", row.names=F, col.names=T)
	
	file.remove("./reference_data/paths/biopax/2161541.owl")
	file.remove("./reference_data/paths/biopax/record_of_biopax_pathways.txt")
}


#test.biopaxFileNameFromPathName()
test.biopaxFileNameFromPathName<-function(){
	# 	library(RUnit)
	#unit test test.biopaxFileNameFromPathName()
	print("test.biopaxFileNameFromPathName()")
	#positive control
	tpnames = c("Abacavir metabolism")
	bpnames = packageDir:::biopaxFileNameFromPathName(pathNames=tpnames)
	checkEquals(target=1, current=length(bpnames))
	
	#negative control
	tpnames2 = c("Abacavir metabolism", 
							 "not a pathway")
	bpnames2 = packageDir:::biopaxFileNameFromPathName(pathNames=tpnames2)
	
	checkEquals(target=1, current=length(bpnames2))
	
}

# #test.getReactomeBiopax()
# test.getReactomeBiopax<-function(){
# 	# 	tpnames = c("Abacavir metabolism", "Transcriptional Regulation of White Adipocyte Differentiation")
# 	tpnames = c("Abacavir metabolism")
# 	getReactomeBiopax(study=NULL, pathNames=tpnames, verbose=F)
# 	
# }

#test.correctByHgncHelper()
test.correctByHgncHelper<-function(){
	
	# 	symbol_set = rownames(STUDY@results$somatic_mutation_aberration_summary$genesummary)
	symbol_set = c("p53", "flt3", "oct5", "BIKE", "LOK")
	res = correctByHgncHelper(symbol_set=symbol_set, verbose=F)
	
	checkEquals(target=c("TP53", "FLT3", "OCT5", "BMP2K", "STK10"), current=res)
	
}

#test.corsym()
#param userInput Logical flag indicating if unit test should require user input
test.corsym<-function(userInput=FALSE){
	# 		symtabFile = "~/tprog/main_131219/reference_data/current_hugo_table.txt"
	# 		symtab = read.delim(file=symtabFile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
	# 		
	# 		#first remove all the unapproved symbols
	# 		withRows = grepl(pattern="~", x=symtab$Approved.Symbol)
	# 		
	# 		apptab = symtab[!withRows,]
	# 		
	# 		extab = apptab[,c("Approved.Symbol","Previous.Symbols")]
	# 	
	
	
	# 	syms2  = rownames(STUDY@results$somatic_mutation_aberration_summary$genesummary)
	# 	save(syms2, file="./testData/")
	#		save(syms2, file="~/tprog/distribution/pathLayerDistributionVersion2/packageDir/inst/testData/amlSymbolsUncorrected.rda")
	# 	system.file("extdata/gene_symbol_corrections_list.txt", package = "packageDir")
	# 	load("./testData/amlSymbolsUncorrected.rda", verbose=T)
	# 	
	# 	load("./testData/simplesymtable.rda", verbose=T)
	
	data("hgnc.table")
	extab = hgnc.table
	appsyms = extab$Approved.Symbol[1:10]
	
	res = corsym(symbol_set=appsyms, symref=extab, verbose=F)
	checkEquals(target=appsyms, current=res)
	
	htab = getHugoSymbols(curhugofname="./reference_data/current_hugo_table_slim.txt", verbose=F)
	prevSyms = c("CPAMD9","ALD", "ALDL1", "ABC50", "AADACL2")
	correctSymbols = c("A2ML1","ABCD1","ABCD2","ABCF1","AADACL2")
	if(userInput){
		tmp = readline("The next test will require input of the letter y")
		res1 = corsym(symbol_set=prevSyms, symref=htab, verbose=T)
		checkEquals(target=correctSymbols, current=res1)
	}
	
	res2 = corsym(symbol_set=prevSyms, symref=htab, verbose=F)
	checkEquals(target=4, current=sum(correctSymbols!=res2))
	
	#checking corsym using hgnchelper
	load(system.file("testData/amlSymbolsUncorrected.rda", package = "packageDir"), verbose=T)
	path_detail = getDefaultPaths()
	path_detail$symtable = hgnc.table
	res3 = corsym(symbol_set=syms2, symref=path_detail, verbose=F)
	checkEquals(target=15, current=sum(res3!=syms2), msg="check of HGNChelper-based script")
}