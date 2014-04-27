#accessory functions
if(!require("RUnit")){
	install.packages("RUnit")
	library("RUnit")
}
if(!require("tcltk")){
	install.packages("tcltk")
	library("tcltk")
}

acc_loaded=T
source('./summaryTable5.R')
source('./pathway_functions.R')
source('./InitiateDataStructures.R')
source('./OverlapAnalysisFunctions.R')
source('./path_paint.R')
source('./save_and_load_data.R')
source('./SettingsObjectInterface.R')

# source("http://bioconductor.org/biocLite.R")
# biocLite("RCytoscape")
library("RCytoscape")
library("graphite")

if(!exists("VERBOSE")) VERBOSE = F

.pardefault <- par(no.readonly = T)
.parpin <-par()$pin
# .pardefault <- par()

promptNumeric <- function (prompt) {
	while(T){
		line=readline(prompt)
		if(line=="") break
		line = suppressWarnings(expr=as.numeric(line))
		if(!is.na(line)) break
		print("Sorry, your input could not be understood, please try again")
	}
	return(line)
}

longTextBarPlot<-function(data, lab, main=""){
	
	bpdata  = barplot(data, horiz=T, main=main)
	
	xmin = par("usr")[1]
	xrange = par("usr")[2] - par("usr")[1]
	xcoord = xmin + (xrange/40)
	
	nlev = length(data)
	
	text(x=xcoord, y=bpdata, labels=lab, pos=4)
	
}

save.plot<-function(pname){
	return(AutoSavePlot(pname))
}

AutoSavePlot<-function(pname){
	if(is.null(pname)){
		pname = "tmp_image_file"
	}
	if(!grepl(pattern=".png$", x=pname)) pname = paste(pname,".png",sep="")
	pfname = plotFileName(pname)
	dev.copy(png,pfname)
	dev.off()
	return(pfname)
}

plotFileName<-function(pname){
	#takes a name, makes sure there's a directory there, adds time stamp to make it unique, adds a postfix to make sure  
	autoincrament = F
	root="./output/imageTemp/"
	if(grepl(pattern="[/]",x=pname)){#the a file path was passed, find the root
		root = paste(paste(strsplit(x=pname, split="[/]")[[1]][1:(length(strsplit(x=pname, split="[/]")[[1]])-1)],collapse="/"), 
								 "/", 
								 sep="")
		pname = strsplit(x=pname, split="[/]")[[1]][length(strsplit(x=pname, split="[/]")[[1]])]#the last part of the path is the actual file name
	}
	dir.create(path=root, recursive=T, showWarnings=F)
	
	fullPath = paste(root, "graphic",gsub(pattern="[-:]", replacement=".", x=as.character(Sys.time())),".", pname, sep="")
	#check if it exists, if it does, append a number
	while(file.exists(fullPath)&autoincrament){
		pname = paste("1", pname, sep=".")
		fullPath = paste(root, "graphic",gsub(pattern="[-:]", replacement=".", x=as.character(Sys.time())),".", pname, sep="")
	}  
	return(fullPath)
}

#pgvmFromStacked makes a patient gene value matrix from data in stacked format
#rows = genes, columns = patients, values = value for gene, ie, for a gene which has a missense mutation, 
# the value would be "missense"
#if more than one value are given for a particular patient-gene combination 
#takes: stackedData: data frame with patient gene data in stacked format
#				valueCol: 		the name of the column with the gene values
#				patientIdCol: the name of the column with patient ids
#				geneCol: 			the name of the column with the gene identifiers
#				blankSymbol: 	what is placed in the pvm when the stacked file does not give a value
#											for a given patient-gene combination
#returns: patient gene value matrix
pgvmFromStacked<-function(stackedData, valueCol, patientIdCol, geneCol, blankSymbol=NA){
	#1: find list of unique pids
	pids = unique(stackedData[,patientIdCol])
	#2: find lis of unique genes
	ugenes = unique(stackedData[,geneCol])
	#case to handle: if patient does not have a value for a gene, 
	#it will be set to NULL
	pvm = matrix(data = blankSymbol, ncol=length(pids), nrow=length(ugenes), dimnames=list(ugenes, pids))
	#3: put values 
	for(p in pids){
		prows = stackedData[stackedData[,patientIdCol] == p,]
		for(g in ugenes){
			pgrows = prows[prows[,geneCol]==g,]
			if(nrow(pgrows) == 1){
				pvm[g,p] = pgrows[1,valueCol]
			}else if(nrow(pgrows) > 1){
				#case to handle: if there are more than one value for a gene.. 
				cat("\nNote: multiple values found for this patient-gene combo:", p, g, "\n")
				pvm[g,p] = paste(pgrows[,valueCol], collapse=" ")
			}
		}
	}
	return(pvm)
}

PGMFromVector<-function(genevector){
	if(is.data.frame(genevector)) genevector = as.matrix(genevector)
	if(is.matrix(genevector)) genevector = as.vector(genevector)
	pgmout = matrix(data=T, 
									nrow=length(genevector), 
									ncol=1, 
									dimnames=list(genevector, "Genes_in_set"))
	return(pgmout)
}

#selectionList()
#prompts the user with a numbered list of selection options
#takes a vector of options to choose from (does not have to be unique set)
#returns logic vector indicating which items in the input vector match the/those 
#					item(s) chosen by user
selectionList<-function(valcol){
	usel = unique(valcol)
	useltmp = gsub(pattern="_",replacement=" ",x=usel)
	uselmat = matrix(data=1:length(usel), ncol=1, dimnames=list(useltmp, "selection number"))
	print(uselmat)
	uin = readline("Please enter the number(s) corresponding to your selection(s)\n")
	uin = as.integer(strsplit(x=uin, split=" ")[[1]])
	lvout = valcol%in%usel[uin]
	return(lvout)
}

filePrompt<-function(defaultfile){
	#prompts user to select file for data input
	#provides default file option
	#returns file name
	fsel = readline(paste("\nTo select a data file, Enter s\n",
												"To load the default data from \n",defaultfile,",\njust press enter \n",sep=""))
	if(fsel=="s"){
		pfile = file.choose()
	}else{
		pfile = defaultfile
	}
	cat("\nLoading data from:\n",pfile,"\n")
	return(pfile)	
}

open.PGM<-function(fname = NULL){
	#opens a patient gene matrix file
	#takes: 1) file name
	#				or 2) no arg (in this case, user will be prompted to select a file)
	#returns: patient gene matrix
	
	if(is.null(fname)){
		fname = file.choose()
	}
	pgmdat = read.table(file=fname,check.names=T,
											header=T,
											row.names=1,
											sep="\t",
											stringsAsFactors=F)
	pgmdat = as.matrix(pgmdat)
	return(pgmdat)
}

#extract_pid_w_matchnormal
#like extract_pid but leaves the sample type code at the end (match/normal)
#takes: n: a string that might contain a TCGA id, ex: TCGA-IQ-7632-01 in unc.edu__IlluminaHiSeq_RNASeqV2__TCGA-IQ-7632-01A-11R-2081-07__expression_rsem_isoforms_normalized.txt
#returns: the TCGA id (with -'s replaced by .'s) if an id is found; NULL if there is not a TCGA id in the input string
extract_pid_w_matchnormal<-function(n=NULL)
{
	if(is.null(n)){n="unc.edu__IlluminaHiSeq_RNASeqV2__TCGA-IQ-7632-01A-11R-2081-07__expression_rsem_isoforms_normalized.txt"}
	#first replace all special characters with spaces
	#en = gsub(pattern="[[:punct:]]", replacement=" ", x=n)
	sen = strsplit(n, split="[[:punct:]]")[[1]]
	#find where "TCGA" is
	i = grep(pattern="TCGA", x=sen)
	if(length(i)<1){
		return(NULL)
	}
	#take out the TCGA part, plus the next two parts
	sub = sen[c(i, i+1, i+2, i+3)]
	sub[4] = gsub(pattern="[a-zA-Z]*", "", x=sub[4])
	out = paste(sub, sep=".", collapse=".")
	return(out)
}

#cleanGeneSymbols
#cleanes input gene symbols, removing leading and trailing spaces, quotes, 
#repaces spaces with dashes and periods with dashes
#takes & returns vector of gene symbols
cleanGeneSymbols<-function(genes){
	#remove leading and trailing spaces
	ltspace2 = gsub(pattern="^ ",replacement="",x=genes)
	#remove trailing spaces
	ltspace2 = gsub(pattern=" $",replacement="",x=ltspace2)
	#replace " " with "-"
	spfixed2 = gsub(pattern=" ",replacement="-",x=ltspace2)
	spfixed2 = gsub(pattern="\\.",replacement="-",x=spfixed2)
	spfixed2 = gsub(pattern="\"",replacement="",x=spfixed2)
	out = toupper(spfixed2)
	return(out)
}

#corListCheck
#if the HUGO symbols reference file has changed, 
#corListCheck will assure the symbol changes are consistent with the new reference file
#takes:		cl:		corrections list
#					htab:	hugo reference table 
#returns: corrections list data frame
#
#works with files: "./reference_data/gene_symbol_corrections_list.txt"
#
corListCheck<-function(cl=NULL, htab=NULL){
	cat("\nScreening symbol correction table...\n")
	curhugofname = "./reference_data/current_hugo_table.txt"
	correctionsfile="./reference_data/gene_symbol_corrections_list.txt"
	if(is.null(cl)){
		cl = read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
	}
	if(is.null(htab)){
		cat("\nOpening HUGO symbol reference table...\n")
		htab=read.table(file=curhugofname,#"./reference_data/current_hugo_table.txt","./reference_data/slim_current_hugo_table.txt"
										header=T, 
										sep="\t", 
										quote="", 
										comment.char="", 
										stringsAsFactors=F,na.strings="-")
	}
	#handle three cases:
	#1: symbol added to HUGO ref
	#2: symbol removed from HUGO ref
	#3: symbol changed
	
	#for 2 and 3: 
	#check if any of the "new" are withdrawn:
	wnew = paste(cl$new_symbol, "~withdrawn", sep="")
	
	iswithdrawn = wnew%in%htab$Approved.Symbol
	if(sum(iswithdrawn)){#if some of the symbols have been withdrawn, see if replacement symbols can be found
		cltmp = cbind.data.frame(cl[iswithdrawn,], wnew[iswithdrawn], stringsAsFactors=F)
		colnames(cltmp)[3] = "withdrawn"
		#for each row in cltmp, merg it's corresponding row in htab
		htabtmp = htab[htab$Approved.Symbol%in%cltmp$withdrawn,1:8]#extract the rows from htab
		chmerge = merge(x=cltmp, y=htabtmp, by.x="withdrawn", by.y="Approved.Symbol")
		####### issue 2
		###### check if some symbols were withdrawn and there is no replacement
		if(sum(chmerge$Status=="Entry Withdrawn")){
			cat("\nOf the symbols in the symbol correction table, these entries were\nwithdrawn as official HUGO sybols, and no replacements were provided:\n")
			print(chmerge$new_symbol[chmerge$Status%in%"Entry Withdrawn"])
		}
		chmerge = chmerge[chmerge$Status == "Symbol Withdrawn",]
		if(nrow(chmerge)){
			tmp= sapply(chmerge$Approved.Name, function(x) strsplit(x=x, split="symbol withdrawn, see ")[[1]][2])
			chmerge$new_symbol = tmp
			clt = rbind(cl, chmerge[,c("old_symbol","new_symbol")])
			clt = clt[!duplicated(x=clt$old_symbol, fromLast=T),]
			cl = clt
			cl = corListCheck(cl=cl,htab=htab)#recursive call, 
			#handles case where symbol was changed more than once
		}
	}
	### now scan the old_symbols for appoved symbols
	#remove rows where the old_symbol is an approved symbol
	if(sum(cl$old_symbol%in%htab$Approved.Symbol)){
		cat("\nSome rows in the symbol correction table were ",
				"found to correct symbols that were already approved.",
				"These rows will be removed to prevent data inconsistencies\n")
		cl = cl[!cl$old_symbol%in%htab$Approved.Symbol,]
	}
	return(cl)
}

toPGMWithCoverage<-function(sds, coverageSet){
	
	pgm = toPGM(sds)
	
	additionalRowNames = setdiff(coverageSet, rownames(pgm))
	newChunk = matrix(data=F, 
										nrow=length(additionalRowNames), 
										ncol=ncol(pgm), 
										dimnames=list(additionalRowNames, colnames(pgm)))
	
}

##toPGM
# converts a stacked format to a patient gene matrix
#				
#takes: 		a data frame with columns "Symbol" and "Ids"
#returns: 	patient gene matrix: rownames = gene names, 
#																	colnames = patientIDs
#																cell values: logical, T if gene is active in that patient
toPGM<-function(sds){
	#find unique patient ids
	pids = unique(sds$Ids)
	pids = make.names(pids)
	sds$Ids = make.names(sds$Ids)
	#find unique genes affected
	syms = unique(sds$Symbol)
	#make the matrix
	out = matrix(data=F, nrow=length(syms), ncol=length(pids), dimnames=list(syms,pids))
	#for each patient, find the set of genes
	for(p in pids){
		#get the current patient's active genes
		cursyms = sds$Symbol[sds$Ids == p]
		#assign the cells for that patient, for those genes
		out[cursyms, p] = T
	}
	return(out)
}



# toBipartateGraph
# takes a data frame with two columns
# makes a bipartate graph from the two columns

toBipartateGraph<-function(dfin){
	colnames(dfin)<-c("Ids", "Symbol")
	#find unique patient ids
	pids = unique(dfin$Ids)
	pids = make.names(pids)
	dfin$Ids = make.names(dfin$Ids)
	#find unique genes affected
	syms = unique(dfin$Symbol)
	#make the matrix
	out = matrix(data=F, nrow=length(syms), ncol=length(pids), dimnames=list(syms,pids))
	#for each patient, find the set of genes
	for(p in pids){
		#get the current patient's active genes
		cursyms = dfin$Symbol[dfin$Ids == p]
		#assign the cells for that patient, for those genes
		out[cursyms, p] = T
	}
	return(out)
}

paths<-function(smd){
	return(smd@paths)
}

#rounds columns in table to "figs" number of significant digits
#and replaces underscores in col names with spaces
cleanTables<-function(tab,figs=3, verbose=F){
	cat("\nCleaning table..")
	if(!(is.data.frame(tab)|is.matrix(tab))){
		cat("WAIT, that's not a table..\n")
		return(tab)
	}
	if(nrow(tab)){
		
		for(i in 1:ncol(tab)){
			ccol = tab[,i]
			cname = colnames(tab)[i]
			if(verbose){
				print(cname)
				print(tab[1,i])
				print(typeof(ccol))
				print(mode(ccol))
			}
			if(typeof(tab[,i])%in%c("integer","double")){
				tab[,i] = round(tab[,i],digits=figs)
			}
		}
	}else{
		print("table has no rows")
	}
	if(ncol(tab)){
		#replace underscores:
		colnames(tab)<-gsub(pattern="[_]", replacement=" ", x=colnames(tab))
	}else{
		print("table has no columns")
	}
	
	cat("cleaned\n")
	return(tab)
}


#getHugoSymbols()
#paths_detail: paths object: if passed, this function will excise the currently used hugo symbol set and return them.<them?>
#curhugofname: the name of the 
#returns HUGO lookup table
#
#function allows re-download of hugo cross ref file. 
getHugoSymbols<-function(paths_detail=NULL, 
												 curhugofname="./reference_data/current_hugo_table_slim.txt",
												 verbose=F){
	if(class(paths_detail)=="Study"){
		cat("\nGetting HUGO gene symbols from study..\n")
		paths_detail=paths_detail@studyMetaData@paths
	}
	if(is.null(paths_detail)){
		cat("\nLoading official HUGO gene symbols from file..\n")
		
		cref = read.table(file=curhugofname,
											sep="\t",
											comment.char="",
											header=T,
											quote="", 
											stringsAsFactors=F, 
											na.strings="-")
		
		cref$Approved.Symbol = cleanGeneSymbols(genes=cref$Approved.Symbol)
		cat("\nUsing a symbol correction file downloaded from http://www.genenames.org/",
				"on",as.character(file.info(curhugofname)$mtime),".\n")
		ginfo = ""
		if(verbose){
			ginfo=readline(paste("To get info on how to update this file enter \"i\"\n",
													 "Otherwise, just press enter to continue"))
		}else{
			cat("\nTo get info on how to update HUGO symbol cross reference file \n",
					"or automatically re-download cross ref file, re-run the getHugoSymbols() function\n",
					"using the argument verbose=T (default: verbose=F)\n")
		}
		if(ginfo=="i"){
			cref = getHugoDownloadInfo(curhugofname = curhugofname)
		}
		return(cref)
	}else{
		cat("\nGetting HUGO symbols from Path_Detail reference calss object")
		return(paths_detail$HUGOtable)
	}#if/else
}#getHugoSymbols


getHugoDownloadInfo <- function (curhugofname) {
	
	if("r"==readline("To attempt to download a current cross reference table of HUGO symbols enter r \n(note: this can take more than 10 minutes to download)\nIf you tried this once and it didn't work, press enter to get other options. ")){
		full_hurl = "http://www.genenames.org/cgi-bin/hgnc_downloads?title=HGNC+output+data&hgnc_dbtag=on&preset=all&status=Approved&status=Entry+Withdrawn&status_opt=2&level=pri&=on&where=&order_by=gd_app_sym_sort&limit=&format=text&submit=submit&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag"
		cat("\nConnecting to HGNC website...\n")
		reopenedfurl=try(expr=read.table(file=full_hurl,sep="\t",comment.char="",header=T,quote="", stringsAsFactors=F,na.strings="-"),
										 silent=T)
		if(sum(grep(pattern="error", x=class(reopenedfurl), ignore.case=T))){
			cat("\nError text:\n")
			print(reopenedfurl)
			cat("\nError, could not download Hugo reference table from HGNC..\n")
		}else{
			cat("\nTable downloaded, writing to file...\n")
			write.table(x=reopenedfurl,file=curhugofname,quote=F,sep="\t")
			cref = read.table(file=curhugofname,sep="\t",comment.char="",header=T,quote="", stringsAsFactors=F, na.strings="-")
			cref$Approved.Symbol = cleanGeneSymbols(cref$Approved.Symbol)
			return(cref)
		}
	}
	
	cat("Go to http://www.genenames.org/ and find the biomart interface.\n",
			"Make sure to check the boxes to download Approved symbol, previous symbol\n",
			"synonyms, status, approved name, date approved and name synonyms.\n")
	cat("Make sure column names match those in the file",curhugofname,"\n",
			"and replace that file with the one downloaded (it is suggested that\n",
			"you rename the old file so as not to loose it if anything goes wrong\n")
	cat("These are the column names from the current hugo cross reference file:")
	print(colnames(cref))
	readline("Press enter to continue.\nPress escape to exit the program so that you can update the cross ref file.")
	cref = getHugoDownloadInfo(curhugofname)
	return(cref)
}



#corsym
#corrects symbols in somatic data
#takes		symbol_set: 1 of 2 options
#																1)table with columns: 
#																"Hugo_Symbol": the set of symbols which should be checked for symbols that need correction
#																col2: 			 The chromosome in the genome that the gene symbols are associated with
#																2)vector of symbols to be corrected
#					curhugofname: the relative file path to a hugo lookup table
#					verbose: if this is set to true, no text output will be given and previously official symbols will not be checked. 
#
#returns: vector: set of gene symbols, corrected to hugo symbols
#
#usage: som_select[,"Hugo_Symbol"] = corsym(som_select[,c("Hugo_Symbol", col2)], "./reference_data/hugo_dl.txt", verbose=verbose)
corsym<-function(symbol_set, hugoref=NULL, verbose=T, col2="Chrom", correctionsfile="./reference_data/gene_symbol_corrections_list.txt"){
	
	if(verbose) cat("\n\nChecking that gene symbols match official HUGO gene symbols. . . . . \n")
	curhugofname="./reference_data/current_hugo_table.txt"
	if(is.null(hugoref)){
		cref=getHugoSymbols(curhugofname=curhugofname)
	}else if(is.character(hugoref)){#if hugoref is a character, it is a file name; open the file
		cref=getHugoSymbols(curhugofname=hugoref)#the direct passing of a hugo file name.  . .not sure if we're doing that any more...
	}else if(class(hugoref)=="Study"){
		if(STUDY@studyMetaData@geneIdentifierType!="HUGO"){
			message("Checking and correction of non-HUGO symbols not implemented.")
			return(symbol_set)
		} 
		cref=hugoref@studyMetaData@paths$HUGOtable
		hugoref=cref
	}else if(class(hugoref)=="data.frame"){#the hugo table was passed directly
		cref = hugoref
	}
	
	#make sure the input symbol set is a table
	if(is.vector(symbol_set)){
		symbol_set = cbind(symbol_set,rep("",times=length(symbol_set)))
	}
	#correct the column names
	if(!sum(c("Hugo_Symbol",col2)%in%colnames(symbol_set))){
		colnames(symbol_set)<-c("Hugo_Symbol",col2)
	}
	
	if(verbose){
		cat(sep="","\n\n",nrow(cref)," symbol records loaded from ",curhugofname,
				"\nNote: these do not all correspond to currently approved HUGO symbols.",
				"\nSome are place holders from previously used or withdrawn symbols\n")
	}
	
	###########################################################
	# Correct gene names
	###########################################################
	#check how many symbols are in both TCGA somatic data and HUGO
	#cross refernce table : if chasm output cant so easily be matched to chasm input
	
	#change all symbols to upper case: 
	cref[,"Approved.Symbol"] = toupper(cref[,"Approved.Symbol"])
	symbol_set[,"Hugo_Symbol"] = toupper(symbol_set[,"Hugo_Symbol"])
	if(verbose){
		cat("\nTo conduct symbol comparrisons and corrections, these reformattings were made:")
		cat("\nConversion to upper case.\nRemoval of leading and trailing spaces.",
				"\nConversion of spaces and periods to dashes.\n",
				"***NOTE: These reformattings are not saved or recorded to any file!!\n")
	}
	
	symbol_set[,"Hugo_Symbol"]= cleanGeneSymbols(symbol_set[,"Hugo_Symbol"])
	
	tcga_hugo = intersect(cref[,"Approved.Symbol"], symbol_set[,"Hugo_Symbol"])
	if(verbose){
		cat(paste("\nOut of the ", as.character(length(unique(symbol_set[,"Hugo_Symbol"])))," symbols in current data set,\n", 
							length(tcga_hugo), " were found to be currently approved HUGO symbols.\n",sep=""))
	}
	not_approvedi = which(!symbol_set[,"Hugo_Symbol"] %in% cref[,"Approved.Symbol"])
	not_approved = symbol_set[not_approvedi,c("Hugo_Symbol", col2), drop=F]
	not_approved = unique(not_approved)
	
	if(!max(0,nrow(not_approved))){#if all symbols match approved hugo symbols, this will be true
		return(symbol_set[,"Hugo_Symbol"])
	}
	if((max(0,nrow(not_approved))<50) & verbose){
		cat("\nThese are the symbols that were not found to be approved HUGO symbols:\n")
		print(not_approved)
	}else if(verbose){
		cat("\nThese are the first 50 out of",max(nrow(not_approved),0),"symbols that were not found to be approved HUGO symbols:\n")
		print(not_approved[1:50,])
	}
	
	#check how many appear to be Micro RNA genes
	possible_miRNAs = grep(pattern="^MIR|-MIR",x=not_approved,ignore.case=T)
	if(length(possible_miRNAs)&verbose){
		cat("\nOf the aformentioned symbols,",length(possible_miRNAs),"symbols begin with \"MIR-\", or contain the string \"-MIR\"",
				"and thus appear to be symbols for microRNAs, as opposed to protein-coding genes.\n",
				"Note: it is possible that gene symbols not containing \"MIR\" correspond to microRNAs as well.\n")
	}
	
	if(verbose){cat("\nA full list of symbols found not to be approved, before attempted",
									"correction, can be found at ./output/not_approved_symbols_from_last_run.txt\n")}
	
	write.table(x=not_approved,file="./output/not_approved_symbols_from_last_run.txt",sep="\t",row.names=F)
	
	correction_set = NULL#the set of mappings from old names to new names, which is retreived from the corrections file
	correctedSymbols = NULL#this will contain mappings to be made in this run of the program
	new_corrections = NULL
	
	###################################  Check Corrections File
	if(file.exists(correctionsfile)){
		raw_correction_set = read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
		resave = F
		#check for repeated rows
		dupcheck = duplicated(raw_correction_set)
		if(sum(dupcheck)){
			raw_correction_set = raw_correction_set[!dupcheck,]
			resave=T
		}
		
		#check that individual old symbol do not have multiple entries
		dupcheck2 = duplicated(raw_correction_set$old_symbol)
		if(sum(dupcheck2)){			
			allolddup = dupcheck2 | duplicated(raw_correction_set$old_symbol,fromLast=T)
			cat("\nWarning, multimapping issue found!\n",
					"These symbol corrections indicate multiple, different corrections for the same symbols:\n")
			print(raw_correction_set[allolddup,])
			readline("Please edit these symbols in the symbol correction file and re-run the program to prevent errors\nPress enter to continue")
		}
		
		if(sum(raw_correction_set$old_symbol == raw_correction_set$new_symbol)){#if there are symbols that are changed to the same thing
			#remove symbols that are the same btwx old and new
			chsymindex = which(raw_correction_set$old_symbol != raw_correction_set$new_symbol)
			correction_set = raw_correction_set[chsymindex,]
			resave=T
		}
		raw_correction_set_tmp=corListCheck(cl=raw_correction_set, htab=hugoref)
		if(sum(!all.equal(target=raw_correction_set_tmp, current=raw_correction_set)==T)){
			resave = T
		}
		raw_correction_set = raw_correction_set_tmp
		correction_set = raw_correction_set
		
		numcor = intersect(x=correction_set[,1], not_approved[,"Hugo_Symbol"])#numcor is the number of symbols in symbols set 
		#that can be corrected from the symbol correction file 
		if(length(numcor)){
			if(verbose){
				cat("\nA previously made corrections file was found at",correctionsfile,"\nThis file contains corrections for ", 
						as.character(length(numcor)), "of the ",nrow(not_approved),"unapproved gene symbols.\n")
			}
			use_previous=""
			if(verbose){
				use_previous = readline(paste("Press enter to use these corrections.\nEnter anything else to skip using these corrections\n", 
																			"(You will be provided with a chance to select your own corrections): "))
			}
			if(use_previous==""){
				symbol_set[,"Hugo_Symbol"] = swapsymbols2(corrected=correction_set, genelist=symbol_set[,"Hugo_Symbol"])
				not_approved = which(!symbol_set[,"Hugo_Symbol"] %in% cref[,"Approved.Symbol"])#temporary state of not_approved
				not_approved = symbol_set[not_approved,c("Hugo_Symbol", col2), drop=F]#not approved now has two columns
				not_approved = unique(not_approved)
				if(verbose) cat("\nThere are now", as.character(max(0,nrow(not_approved))), "symbols remaining which do not match approved HUGO symbols.")
			}
			print("not approved 1.2:")
			print(not_approved)
			if(resave){
				write.table(x=correction_set, file=correctionsfile, 
										quote=F, sep="\t", row.names=F, col.names=c("old_symbol", "new_symbol"))
				raw_correction_set = read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
			}
		}
		
		if(max(0,nrow(not_approved))==0){
			return(symbol_set[,"Hugo_Symbol"])
		}
	}
	
	if(verbose){
		if(max(0,nrow(not_approved))){
			###################################   Check previously used symbols
			checkprev = readline("\nWould you like to check previously official HUGO symbols for the remaining unmatching symbols? \n(enter y or n)")
			if(checkprev=="y"){
				######## Check previous HUGO symbols
				switches = checkPreviousSymbols(symbols=not_approved, indexes=1:nrow(not_approved), hugolookup=cref, col2=col2)
				if(nrow(switches)){
					print(switches)
					symbol_set[,"Hugo_Symbol"]=swapsymbols2(corrected=switches, genelist=symbol_set[,"Hugo_Symbol"]) 
					not_approved = which(!symbol_set[,"Hugo_Symbol"] %in% cref[,"Approved.Symbol"])#temporary state of not_approved
					not_approved = symbol_set[not_approved,c("Hugo_Symbol", col2)]#not approved now has two columns
					not_approved = unique(not_approved)
					new_corrections = rbind(switches, new_corrections)
					cat("\n",nrow(switches),"mystery symbols were found in the previously official hugo symbols.\n")
					cat("There is/are now", as.character(nrow(not_approved)), "symbol(s) remaining which do/does not match approved HUGO symbols.\n")
					print(not_approved)
				}else{
					cat("\nNo matches were found in the previously used HUGO symbols\n")
				} 
			}	
		}
		print(not_approved)
		if(max(nrow(not_approved),0)){
			###################################   Check synonyms
			checksyn = readline("\nWould you like to check synonyms for the remaining unmatching symbols? (enter y or n) ")
			if(checksyn=="y"){
				
				syncor = checkSynonyms(symbols=not_approved, indexes=1:nrow(not_approved), hugolookup=cref, col2=col2)
				if(nrow(syncor)){
					colnames(syncor)<-c("old_symbol","new_symbol")
					symbol_set[,"Hugo_Symbol"]=swapsymbols2(corrected=syncor, genelist=symbol_set[,"Hugo_Symbol"])
					not_approved = which(!symbol_set[,"Hugo_Symbol"] %in% cref[,"Approved.Symbol"])#temporary state of not_approved
					not_approved = symbol_set[not_approved,c("Hugo_Symbol", col2), drop=F]#not approved now has two columns
					not_approved = unique(not_approved)
					new_corrections = rbind(new_corrections, syncor)
				}
				cat("\nThere is/are now", as.character(nrow(not_approved)), "symbol(s) remaining which do not match approved HUGO symbols.\n")
				print(not_approved)
			}
			
			if(verbose){cat("\nA full list of symbols which were found not to be approved and which could\n",
											"not be corrected can be found at ./output/not_approved_not_correctable_symbols_from_last_run.txt\n")}
			write.table(x=not_approved,file="./output/not_approved_not_correctable_symbols_from_last_run.txt",sep="\t",row.names=F)
			if(length(new_corrections)){
				cat("\nThese are the new corrections that were just added to the symbol corrections file:\n")
				print(new_corrections)
				if(readline("Would you like to save the set of corrections just made to the corrections file (y/n)")=="y"){
					addCorrections(new_corrections=new_corrections, correction_set=correction_set, correctionsfile=correctionsfile)
				}
			}
		}
	}#if verbose
	return(symbol_set[,"Hugo_Symbol"])
}#corsym function


#'@title addCorrections()
#'@description Adds corrections to the corrections file
#'@param new_corrections: a two column matrix; column 1 = old, incorrect symbols, column 2 = new, corrected symbols
#'@param correctionsfile: the file name of the corrections file
#'@param correction set: the original contents of teh corrections file before new corrections were made
addCorrections<-function(new_corrections, correctionsfile="./reference_data/gene_symbol_corrections_list.txt", correction_set=NULL){
	if(is.null(correction_set)){
		if(file.exists(correctionsfile)){
			correction_set = as.matrix(read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-"))
		}else{
			correction_set = matrix(data="", ncol=2, nrow=0, dimnames=list(NULL, c("old_symbol", "new_symbol")))
		}
	} 
	colnames(new_corrections)<-colnames(correction_set)
	final_corrections = rbind(correction_set, new_corrections)
	#screen for NA values
	final_corrections = final_corrections[!is.na(final_corrections[,1]),,drop=F]
	final_corrections = final_corrections[!is.na(final_corrections[,2]),,drop=F]
	#now attempt to clean the final_collections of rows where no changes are made (ie, old symbol is the same as new symbol)
	final_corrections=final_corrections[final_corrections[,1]!=final_corrections[,2],,drop=F]
	final_corrections = unique(final_corrections)
	write.table(x=final_corrections, file=correctionsfile, 
							quote=F, sep="\t", row.names=F, col.names=c("old_symbol", "new_symbol"))
	cat("\nSymbol corrections recorded in file:", correctionsfile, "\n")	
}

#examineHugoSet
#get summary information on a set of HUGO symbols
#takes: symbol_set: the list of symbols
#study_name and data_type describe the source of the symbols, for display on the read out
#
examineHugoSet<-function(symbol_set,study_name,data_type,curhugofname="./reference_data/current_hugo_table.txt"){
	
	cref = read.delim(file=curhugofname, header=T, stringsAsFactors=F, na.strings="-")
	approvedHugoFile = paste("./output/",study_name,"_approved_hugo_w_annotation_from_",data_type,"_data.txt",sep="")
	
	cat("\nA table of approved HUGO symbols, including those just corrected, can be found in the file",approvedHugoFile,"\n")
	
	in_study = cref[cref[,"Approved.Symbol"]%in%symbol_set,]
	not_hugo = setdiff(x=symbol_set, y=cref[,"Approved.Symbol"])
	nhtab = matrix(nrow=length(not_hugo),ncol=ncol(in_study))
	colnames(nhtab)<-colnames(in_study)
	nhtab[,"Approved.Symbol"] = not_hugo
	nhtab[,"Locus.Type"] = rep("Not HUGO symbol",times=length(not_hugo))
	nhtab[,"Locus.Group"] = rep("Not HUGO symbol",times=length(not_hugo))	
	
	in_study = rbind(in_study, nhtab)
	
	write.table(x=in_study,file=approvedHugoFile,quote=F,sep="\t",row.names=F,col.names=T)
	
	return(in_study)
}

#takes: indexes: indexes of symbols to be checked
#       symbols: list of symbols which the indexes refer to 
#                             (those to be switched and those to remain)
#       hugolookup: hugo look up table with columns : <Approved.Symbols>, <Previous.Symbols>
#returns: table: <index to be switched> <what it should be switched to>
checkPreviousSymbols<-function(symbols, indexes, hugolookup, col2){
	
	hugolookup[,"Previous.Symbols"] = toupper(hugolookup[,"Previous.Symbols"])
	
	oldsyms = NULL #output set of those which are previous symbols
	ssymbols = NULL #
	for(i in 1:length(indexes))
	{
		cur = symbols[indexes[i],"Hugo_Symbol"]
		curchrom = symbols[indexes[i], col2]
		cursym = gsub(pattern="[[:punct:]]", replacement=" ", x=cur)
		
		#grep the row against the previous symbols column
		pat = paste("(^|[[:blank:]])", cursym, "(,|$)", collapse="", sep="")
		ind = grep(pattern=pat, x=hugolookup[["Previous.Symbols"]],ignore.case=T)
		if(length(ind)==1){#if grep caught something
			ssymbols = c(ssymbols, as.character(hugolookup[ind,"Approved.Symbol"]))
			oldsyms = c(oldsyms, cur)
		}else if(length(ind)>1){#if grep found more than one matching previous symbol
			dcols = colnames(hugolookup)[grep("date", colnames(hugolookup), ignore.case=T)]
			dcols = c(dcols,colnames(hugolookup)[grep(col2, colnames(hugolookup), ignore.case=T)])
			cat("\nCorrecting gene symbol", as.character(i), "of", as.character(length(indexes)), "symbols that need correction.\n")
			cat("\nMultiple previous symbols were found to match", as.character(cursym), " (from Chrom ", as.character(curchrom), "). \n")
			cat("\nThese are the symbols that were found to match:\n")
			print(hugolookup[ind,c("Approved.Symbol", "Previous.Symbols", dcols)])
			line = readline("Please enter the correct HUGO name if it appears above under Approved.Symbol.\nPress enter to continue with out making any symbol corrections.")
			if(line != ""){
				ssymbols = c(ssymbols, line)
				oldsyms = c(oldsyms, cur)
			}
		}
	}
	#cat("\n",length(oldsyms),"symbols were found to match previously official symbols, and are being corrected.\n")
	out = cbind.data.frame(oldsyms, ssymbols, stringsAsFactors=F)
	return(out)
}


#takes: indexes: indexes of symbols to be checked
#       symbols: list of symbols which the indexes refer to 
#                             (those to be switched and those to remain)
#       hugolookup: hugo look up table with columns : <Approved.Symbols>, <Previous.Symbols>
#returns: table: <index to be switched> <what it should be switched to>
checkSynonyms<-function(symbols, indexes, hugolookup, col2){
	sindexes = NULL #output set of those which are previous symbols
	ssymbols = NULL
	for(i in 1:length(indexes))
	{
		osym = symbols[indexes[i],"Hugo_Symbol"]
		curchrom = symbols[indexes[i], col2]
		cursym = gsub(pattern="[[:punct:]]", replacement=" ", x=osym)
		sterm=cursym
		acc = NULL
		while(T){
			print("start while(T)")
			#grep the row against the previous symbols column
			pat = paste("(^|[[:blank:]])", sterm, "(,|$)", collapse="", sep="")
			res1 = grep(pattern=pat, x=hugolookup[["Synonyms"]])
			res2 = grep(pattern=sterm, x=hugolookup[["Synonyms"]], ignore.case=T)
			res3 = grep(pattern=pat, x=hugolookup[["Approved.Symbol"]], ignore.case=T)
			cat("\n###################################################### Searching synonyms for", cursym,"\n")
			cat("#####in the data being processed, this symbol is associated with chromosome", curchrom, "\n\n########## Original query:", osym,
					"\n########## Current search term: ", sterm, "\n")
			if(length(res2)>0){
				cat("\n ############ These near matches were found: \n")
				print(hugolookup[res2, c("Approved.Symbol","Approved.Name","Status", "Synonyms", 
																 colnames(hugolookup)[grep("date", colnames(hugolookup), ignore.case=T)], 
																 colnames(hugolookup)[grep(col2, colnames(hugolookup), ignore.case=T)])])
			}else{cat("\nNo synonyms were found for the search term.\n")}
			if(length(res3)==1){
				cat("\n ************ This exactly matching approved HUGO symbol was found: \n")
				print(hugolookup[res3, c("Approved.Symbol","Approved.Name", 
																 colnames(hugolookup)[grep(col2, colnames(hugolookup), ignore.case=T)])])
				if(readline(prompt="If you would like to accept this exactly matching approved HUGO symbol, press ENTER.")==""){
					acc = as.character(hugolookup[res3, "Approved.Symbol"]) 
					break
				}
			}
			if(length(res1)==1){
				cat("\n ************ This exactly matching synonym was found: \n")
				print(hugolookup[res1, c("Approved.Symbol","Approved.Name","Status", "Synonyms", 
																 colnames(hugolookup)[grep("date.approved", colnames(hugolookup), ignore.case=T)], 
																 colnames(hugolookup)[grep(col2, colnames(hugolookup), ignore.case=T)])])
				if(""==readline("If you would like to accept the exactly matching synonym above, please press ENTER\nIf you enter anything else, more options will be provided.")){
					acc = as.character(hugolookup[res1, "Approved.Symbol"]) 
					break
				}
			}
			
			line = readline(paste("\nIf you would like to accept the current search term, \"",
														sterm,
														"\" \nas the symbol correction, please press ENTER.\nTo continue with out making a correction, enter c\nIf you would prefer to enter a correction, please key it in now: "))
			if(line ==""){
				acc = sterm
				break
			}
			if(line == "c"){
				acc = ""
				break
			}
			sterm = line
		}#while
		print("Out of while loop")
		if(acc != ""){
			ssymbols = c(ssymbols, acc)
			sindexes = c(sindexes, osym)
		}
	}#for
	out = cbind.data.frame(sindexes, ssymbols, stringsAsFactors=F)
	if(nrow(out)) colnames(out)<-c("old_symbol","new_symbol")
	print("about to return")
	return(out)
}#checkSynonyms

#qinfo: gives you quick info about a data structure
#takes: ds: vector, dataframe or matrix
qinfo<-function(ds,mx=25){
	print(mode(ds))
	print(typeof(ds))
	print(class(ds))
	if(is.vector(ds)){
		print(head(ds))
		print(summary(ds))
	}else if(is.matrix(ds)|is.data.frame(ds)){
		print(dim(ds))
		if(ncol(ds)<5){
			print(head(ds))
			print(summary(ds))
		}else{
			print(ds[1:min(nrow(ds),5),1:4])
			print(summary(ds[,1:4]))
			print(colnames(ds)[1:min(ncol(ds),mx)])
			print(rownames(ds)[1:min(nrow(ds),mx)])
		}
	}else if(!is.list(ds)){
		print(ds)
	}else if(is.list(ds)){
		print(names(ds)[1:min(mx,length(ds))])
		for(i in 1:length(ds)){
			print(names(ds)[i])
			qinfo(ds[[i]])
		}
	}
}


#extract_pid
#takes: n: a string that might contain a TCGA id, ex: TCGA-IQ-7632-01 in unc.edu__IlluminaHiSeq_RNASeqV2__TCGA-IQ-7632-01A-11R-2081-07__expression_rsem_isoforms_normalized.txt
#returns: the TCGA id (with -'s replaced by .'s) if an id is found; NULL if there is not a TCGA id in the input string
extract_pid<-function(n=NULL)
{
	if(is.null(n)){
		n="unc.edu__IlluminaHiSeq_RNASeqV2__TCGA-IQ-7632-01A-11R-2081-07__expression_rsem_isoforms_normalized.txt"
		cat("\nNote: extracting example pid from\n", n, "\n")
	}
	#first replace all special characters with spaces
	#en = gsub(pattern="[[:punct:]]", replacement=" ", x=n)
	sen = strsplit(n, split="[[:punct:]]")[[1]]
	#find where "TCGA" is
	i = grep(pattern="TCGA", x=sen)
	if(length(i)<1){
		print("** Note: TCGA id not found, returning barcodes as pids.")
		return(n)
	}
	#take out the TCGA part, plus the next two parts
	sub = sen[c(i, i+1, i+2)]
	out = paste(sub, sep=".", collapse=".")
	return(out)
}

test.summarize_by<-function(){
	
	tset1 = c(rep("thirty", times=30), rep("twenty", times=20), rep("not ten", times=11))
	tset1res = summarize_by(col=tset1, display=T, barPlotTitle="test summarize by", left_margin_factor=1.1)
	checkEquals(checkNames=F, 
							current=tset1res, 
							target=cbind.data.frame(stringsAsFactors=F, 
																			types=c("not ten", "thirty","twenty"), counts=c(11,30,20)))
	
}
###summarize_by
#summarize table by categories in a particular column
#used for making a quick summary of different types of somatic mutations in TCGA somatic mutation tables
#takes:   the table
#         the column to be summarized
#         whether or not you want to display the summary to screen and bar graph
#returns: a summary table: category column, count column
#'@title summarize_by
#'@description summarizes the data levels in a column of categorical data and optionally provides a barplot output
#'@param col Vector, the categorical values to be summarized
#'@param display A flag indicating if the data should be displayed in a bar plot. 
#'@param barPlotTitle The title that should be given to the bar plot
#'@param left_margin_factor A numerical factor expanding the left margin of the bar plot. Used to adjust in situations where there is insufficient room to display the text of each bar's label
#'@return A two column data.frame with columns "types", the different levels and "counts", the different numer of times each level occurs
summarize_by<-function(col, display=F, barPlotTitle="Counts across types", left_margin_factor=1)
{
	
	ctab = table(col)
	
	out = cbind.data.frame(rownames(ctab), as.vector(ctab), stringsAsFactors=F)
	colnames(out)<-c("types", "counts")
	out = out[order(out$types),]
	out$counts = as.numeric(out$counts)
	types = out$types
	counts = out$counts
	if(display){
		oldmar <- par()$mar
		while(T){
			res = try({
				par(mar=c(5.1, max(4.1,max(left_margin_factor*nchar(types))/2.5) ,4.1 ,2.1))
				#	try(displayGraph(w), silent=T)
				barplot(counts, horiz=T, las=2, main = barPlotTitle, xlab="Number found in data set", names.arg=types)
				par(oldmar)
			}, silent=T)
			# 			if(!grepl(pattern="Error", x=res)) break
			if(is.null(res)) break
			par(oldmar)
			readline(prompt="There seems to have been an error with plotting the bar graph.\nPlease increase the size of the plot window, the press enter")
		}
		par(oldmar)
		print(out)
	}
	
	return(out)
}

#join_tables
#takes: orig: 				a table to which columns from new_t are to be joined (if this is null t2 will be returned with adjusted names)
#       new_t: 				joins columns from new_t to orig, filling 0s in where new_t doesn't have a value
#				name_prefix: 	prefix to be appended to columns from new_t
#				fill: 				value to be placed where new_t contains no value (typically 0 or NA)
#returns: t1 columns joined to t2 columns, new rows filled with the value in "fill"
joinTables<-function(t1,t2, name_prefix=NULL, fill=0){
	print("inside joinTables()")
	if(is.null(name_prefix)){#invent a name prefix
		name_prefix = as.character(sum(ncol(t1), ncol(t2)))
	}
	#change t2 column names
	t2n = colnames(t2)
	newt2n = paste(name_prefix, t2n, sep=".")
	print(newt2n)
	colnames(t2)<-newt2n
	if(is.null(t1)){
		return(t2)
	}
	allrownames=union(rownames(t1), rownames(t2))
	#make new chunk
	newchunkrows = setdiff(allrownames, rownames(t2))
	tmpnew = data.frame(matrix(ncol=ncol(t2), nrow=length(newchunkrows), data=fill), row.names=newchunkrows)
	colnames(tmpnew)<-colnames(t2)
	newchunk = rbind(t2, tmpnew)
	
	rowstoold = setdiff(allrownames, rownames(t1))
	
	tmpold = data.frame(matrix(ncol=ncol(t1), nrow=length(rowstoold), data=fill), row.names=rowstoold)
	
	colnames(tmpold)<-colnames(t1)
	updatedold = rbind(t1, tmpold)
	
	#join the old and the new
	out = cbind(updatedold, newchunk[rownames(updatedold),])
	cat("\nTables Joined.. \n")
	return(out)	
}

#like swapsymbols, but makes a dictionary out of corrected, not genelist, allowing genelist to contain duplicates
#makes a dictionary out of corrected, to correct the items in genelist
#takes: corrected: two column table, columns: <old symbols> <new symbols>
#       genelist: list containing symbols to be corrected
#returns: genelist with symbols corrected
swapsymbols2<-function(corrected, genelist){
	
	#corrected needs to be subseted so it doesn't add anything
	corrected = corrected[corrected[,1] %in% genelist,]
	
	#figure out which genes were actually corrected; if they were not correted, corrected[,2] or [,1] will contain an NA
	#this gets only the lines that dont have NA in column 1 or 2
	index = which(!(is.na(corrected[,2]) | is.na(corrected[,1])))
	#make switch dictionary; names = old symbol, values = new symbol
	cor = corrected[index,2]
	names(cor) <- corrected[index,1]
	#now switch the uncorrected unique targets to the corrected unique targets
	out = genelist
	rows = which(genelist%in%corrected[,1])#make sure to only correct the ones we have corrections for
	for(i in rows){
		out[i] = as.character(cor[genelist[i]])
	}
	print(length(rows))
	return(out)
}



#getEmpP
#get empirical p-value
getEmpP<-function(path,pgm,reps=1000,paths){
	#some number of genes
	
	#draw random genes for pathway
	#for each gene: see if it's in the pgm
	#								if it is, set it to the score from the pgm
	#								if it's not, set it to zero
	#score = sum for all genes in pathway
	#add score to set of scores for ECDF
	
	flatpgm = pgm%*%rep(T, ncol(pgm))
	ngenes = length(path) #the number of genes in the path
	
	vals = rep(0,times=reps)
	for(i in 1:reps){
		res = rep(x=0,times=ngenes) # results from this iteration
		#select random genes for pathway
		curPath = sample(colnames(paths),size=ngenes,replace=F)
		curInSomatic = rownames(pgm)%in%curPath
		#now get the scores for the genes that are actually in the pathway
		curScores = flatpgm[curInSomatic] #extract the scores for all the genes in the current randomized pathway
		curScore = sum(curScores)
		vals[i] = curScore
	}
	
	genesInPath = rownames(flatpgm)%in%path
	trueScore = sum(flatpgm[genesInPath])
	
	curFunction <- ecdf(vals)
	pout = 1-curFunction(trueScore)
	return(pout)
}#getEmpP()

#goal: provided this set of data for each piece of the project: 
# Panel used: 									Described in Tyner et al Oct. 2012*
# Total number of unique genes targeted by drug or mutation : 				290
# Unique genes from targeted genes found in the provided set of pathways : 		180 
# genes targeted, but not found in the current set of pathways : 	110 
# Pathways used : 								KEGG PATHWAYS / 2011-Mar-14** 
# Number of unique genes found in all pathways : 							5894 
# Total number of unique pathways examined : 								229 
# Number of unique pathways targeted by aberrations or drugs : 							98 
# List of genes targeted, but not found in current set of pathways
# AAK1 ACVRL1 ADCK3 ADCK4 ALK ANKK1 NUAK1 AURKB AURKC AXL BMP2K BLK BMX BRSK1 BRSK2 CAMK1 CAMK1D CAMK1G CDK11B CDK11A CDK19 CDK3 CDK8 CDK9 CIT CLK1 CLK2 CLK3 CLK4 DCLK1 DCLK2 DCLK3 DDR1 DDR2 DMPK CDC42BPG STK17A STK17B DYRK1B MAPK6 MAPK4 MAPK15 FRK GAK LATS1 LATS2 STK10 LTK MAP4K5 MARK1 MARK2 MARK3 MARK4 MELK MERTK MAP3K9 MAP3K10 CDC42BPA CDC42BPB MST1 STK24 MST4 MUSK MYO3A MYO3B STK38L NEK1 NEK2 NEK5 NEK6 NEK7 NEK9 CDK16 CDK17 CDK18 CDK14 PIM3 PKN1 PKN2 PLK3 PLK4 PRKD1 PRKD2 PRKD3 PTK6 RIOK3 ROS1 MYLK4 NUAK2 SIK1 SIK2 SRMS SRPK1 SRPK2 STK16 STK33 TESK1 TIE1 TLK1 TLK2 TNIK TNK1 TNK2 TNNI3K TSSK1B TYRO3 STK32B STK32C STK25 ZAK

startHTMLPlugIns<-function(){
	if(!require("xtable")){
		print("Trying to install xtable so that HTML output can be generated")
		install.packages("xtable")
		if(require("xtable")){
			print("xtable installed and loaded")
		} else {
			stop("could not install xtable")
		}
	}
	library("xtable", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
	if(!require("hwriterPlus")){
		print("Trying to install hwriterPlus so that HTML output can be generated")
		install.packages("hwriterPlus")
		if(require("hwriterPlus")){
			print("hwriterPlus installed and loaded")
		} else {
			stop("could not install hwriterPlus")
		}
	}
	library("hwriterPlus", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
}

WriteSettings<-function(sumset, p){
	dwu = sumset[["settings"]]
	
	if(is.null(dwu)) return()
	if(is.list(dwu)&!is.data.frame(dwu)) dwu = saveSettings(set=dwu)
	print("Inside WriteSettings()... ")
	#WriteDataWorkUp()
	#builds the data work up table and outputs it to the HTML page
	#takes: sumset: the summary set
	#				p: the handle for the html page being created
	#returns: nothing
	
	hwrite(paste("<A name=\"section",".002","\">","Settings used for this data set", "</A>", sep=""),p,heading=2)
	hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
	cat("\nWriting settings to HTML..\n")
	hwrite(dwu, p, table.style= 'border-collapse:collapse', table.cellpadding='5px', col.names=F)
}

WriteDataWorkUp<-function(sumset, p, fpath){
	dwu = sumset[["Data_work_up_notes"]]
	if(is.null(dwu)|!length(dwu)) return()
	
	if(!is.data.frame(dwu)&is.list(dwu)){
		dwu = saveDataWorkUpNotes(set=dwu, fpath=fpath)
	}
	print("Inside WriteDataWorkUp()... ")
	#WriteDataWorkUp()
	#builds the data work up table and outputs it to the HTML page
	#takes: sumset: the summary set
	#				p: the handle for the html page being created
	#returns: nothing
	
	hwrite(paste("<A name=\"section",".001","\">","Data work-up notes", "</A>", sep=""),p,heading=2)
	hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
	cat("\nOutputting data work up notes.\n")
	print(class(dwu))
	print(dim(dwu))
	if(is.data.frame(dwu)|is.matrix(dwu)){
		dwu_tab = dwu
	}else{
		print(dwu)
		dwu_names = names(dwu)
		dwu_vals =c()
		
		for(i in 1:length(dwu)){
			if(length(dwu[[i]])){
				if(length(dwu[[i]][[1]])==1&!is.matrix(dwu[[i]])){
					# 					dwu_names = c(dwu_names,names(dwu)[i])
					dwu_vals = c(dwu_vals, dwu[[i]][[1]])
				}
			}
		}
		print("bout to bind the data fram")
		dwu_tab = cbind.data.frame(dwu_names, dwu_vals,stringsAsFactors=F)
		print("bound")
		
	}
	
	cat("\nwriting work up notes to HTML..\n")
	hwrite(dwu_tab, p, table.style= 'border-collapse:collapse', table.cellpadding='5px', col.names=F)
	
}

#htmlSummary outputs summary information from a summaryTable output list to an HTML file
#the input sumset must have, at the very least:
#															$summarystats: list of summary stats (path coverage, etc. .)
#															$pathsummary: the generalized output table
#															$genesummary: counts of times each gene is targeted across cohort 
htmlSummary<-function(sumset,
											fname="./test.html",
											plimit=.05, 
											sortColumn=5, 
											pagetitle=NULL, 
											path_detail=NULL){
	cat("\nCreating HTML summary. . . \n")
	cat("\nFile name of HTML summary:\n")
	print(fname)
	
	startHTMLPlugIns()
	targetType=sumset$summarystats[2,2]
	headings = c()
	
	#check if there's a coverage summary
	baseDir = paste(strsplit(fname, "/")[[1]][1:(length(strsplit(fname, "/")[[1]])-1)], collapse="/")
	if("coverage_summary"%in%names(sumset)){
		print("Found coverage summary")
		htmlSummary(path_detail=path_detail, sumset=sumset$coverage_summary, fname=paste(baseDir, "coverage_summary.html", sep="/"), path_detail)
		if(!"summarystats"%in%names(sumset)){#if a coverage summary is all there is, then exit
			return()
		}
		cat("Back to writing",fname,"\n" )
	}
	
	if(!file.exists(fname)){
		cat("\nCreating file folder..\n")
		dir.create(path=baseDir, recursive=T,showWarnings=F)
	}
	
	cat("names in sumset:\n")
	print(names(sumset))
	print(ncol(sumset$patientGeneMatrix))
	
	if(is.null(names(sumset$summarystats))) return() #if it's blank, get outta there
	
	print("just about to look at the image slots...")
	multiPatients = F
	if(!is.null(sumset$patientsums)){
		multiPatients = nrow(as.matrix(sumset$patientsums))>1
	}
	
	cat("\nOpening HTML page..\n")
	print(fname)
	p = openPage(fname)
	if(!is.null(pagetitle)) hwrite(gsub(pattern="_",replacement=" ",x=pagetitle), p, heading=1)
	
	headings = c(paste("Summary of",sumset$summarystats[1,2],"used in enrichment testing"),
							 paste("Distribution of",targetType,"genes across patients"), 
							 paste("Top 20 most commonly",targetType,"genes, across cohort"), 
							 "Analysis of pathways", 
							 "Individual patient enrichments for the patients with the top 5 highest mutation rates ", 
							 headings)
	print("headings:")
	print(headings)
	print(multiPatients)
	if(("imageSlots"%in%names(sumset)&!is.null(sumset$imageSlots$imageSlots))|multiPatients){
		#establish image directory: any images will go here
		imageDir = paste(baseDir, "/imageSlots/", sep="")
		cat("\nCreating image directory..\n")
		dir.create(imageDir,recursive=T,showWarnings=F)
		headings = c(headings, "Path images")
	} 
	
	#for all the headings, turn the heading into an anchor
	links = headings
	hwrite("<A name=\"section0.1\">Table of contents</A>", p, heading=2)
	
	if("Data_work_up_notes"%in%names(sumset)) hwrite(paste("<A href=\"#section",".001","\">","Data work-up notes", "</A>", sep=""),p,heading=3)
	if("settings"%in%names(sumset)) hwrite(paste("<A href=\"#section",".002","\">","Settings", "</A>", sep=""),p,heading=3)
	
	for(i in 1:length(headings)){
		links[i] = paste("<A href=\"#section",i,"\">",headings[i], "</A>", sep="")
		headings[i] = paste("<A name=\"section",i,"\">",headings[i], "</A>", sep="")
		hwrite(links[i],p,heading=3)
	}
	
	WriteDataWorkUp(sumset, p, baseDir)
	
	WriteSettings(sumset, p=p)
	
	cat("\nTarget type:",targetType,"\n")
	
	hwrite(headings[1],p,heading=2)
	hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
	print("just wrote heading 1")
	#output the summary stats
	sum_stats=sumset$summarystats
	print(names(sumset))
	print(sum_stats)
	hwrite(sum_stats, p, table.style= 'border-collapse:collapse', table.cellpadding='5px')
	print("Just wrote sum_stats")
	#output most targeted genes
	if(multiPatients){#if more than one patient's worth of gene data is being summarized here, 
		cat("\nOutputting distribution of",targetType,"genes across all patients\n")
		#output information about the counts for each gene across the cohort
		hwrite(headings[2],p,heading=2)
		hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
		#put a histogram into the summary
		imageName = paste("dist_",targetType,"_genes_per_patient.jpg",sep="")
		imgloc = paste(baseDir,"/imageSlots/",imageName,sep="")
		jpeg(imgloc)
		hist(sumset$patientsums[,1],breaks=20,main=paste("Distribution of",targetType,"gene counts for all patients:"),xlab=paste("Number of",targetType,"genes in patient"))
		dev.off()
		hwriteImage(paste("./imageSlots/", imageName, sep=""), p, br=TRUE)
		genelist = sumset$genesummary
		geneindex = order(genelist[,1],decreasing=T)
		# 		ordgenelist = as.matrix(genelist[geneindex[ 1:min( 20,length(geneindex) ) ],],ncol=1)
		ordgenelist = genelist[geneindex[ 1:min( 20,length(geneindex) ) ],,drop=F]
		if(!is.null(sumset$preFilteringGeneSummary)){
			print("!is.null(sumset$preFilteringGeneSummary)")
			genelistpre = sumset$preFilteringGeneSummary
			geneindexpre = order(genelistpre[,2],decreasing=T)
			ordgenelistpre = as.matrix(genelistpre[geneindexpre[1:20],],ncol=1)
			ordgenelist = cbind(rownames(ordgenelist),ordgenelist[,1])
			ordgenelistpre = GenesInPaths(ordgenelistpre, path_detail)
			ordgenelist = GenesInPaths(ordgenelist, path_detail)
			ordgenelist = cbind(ordgenelistpre, ordgenelist)
			print(head(ordgenelist))
			colnames(ordgenelist)<-c("Top mutated genes, before filtering", 
															 "Number of patients with mutated gene, pre filtering",
															 "In the current set of pathways?",
															 "Top mutated genes, after filtering", 
															 "Number of patients with mutated gene, after filtering",
															 "In the current set of pathways?")
		}else{
			ordgenelist = GenesInPaths(ordgenelist, path_detail)
		}
		print("Bout to print heading 3")
		hwrite(headings[3],p,heading=2)
		hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
		hwrite(ordgenelist, p, table.style= 'border-collapse:collapse', table.cellpadding='1px', row.names=is.null(sumset$preFilteringGeneSummary),col.names=T)
	}
	
	hwrite("<br><br>",p)
	#pathway analysis/coverage/enrichment
	hwrite(headings[4],p,heading=2)
	hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
	
	#output the significant paths
	ptable = sumset$pathsummary
	
	tmpTest = F
	if(is.null(ptable)){
		tmpTest = T
	}else if(!ncol(ptable)|!nrow(ptable)){
		tmpTest=T
	}
	
	if(tmpTest){
		print(sumset$summarystats[1,2])
		# 		readline("empty table")
		hwrite(paste("There are no",targetType,"genes in this patient.",p))
	}else{
		ptable = cleanTables(ptable)
		colnames(ptable) = gsub(pattern="_",replacement=" ",x=colnames(ptable))
		if(ncol(ptable)<6){
			print("Coverage summary.")
			sortColumn = 4
		}
		hwrite(paste("Pathways sorted by",colnames(ptable)[sortColumn],"."),p,heading=2)
		torder = order(ptable[,sortColumn],decreasing=T)
		hwrite(ptable[torder,], p, table.style= 'border-collapse:collapse', table.cellpadding='1px', row.names=F)
		
	}
	
	#determine, for the patients with the most mutation, in the paths that are most significantly 
	#mutated, which genes are mutated the most
	#find genes and pathways mutated in patients with highest overall mutational load
	if(length(sumset$path_summary_each_patient)){
		hwrite(headings[5],p,heading=2)
		hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
		print("About to write info about top 5 most affected patients.. . ")
		if(length(sumset$patientsums)>1){
			#first: select to 5 patients
			ordmurates = sumset$patientsums[order(sumset$patientsums,decreasing=T),]
			print(head(ordmurates))
			top5pat = names(ordmurates[1:5])#patients with the top 5 highest mutation rates
			#second: select pathways most enriched paths in those patients
			for(ci in 1:length(top5pat)){
				cp = top5pat[ci]
				cat("\nPatient: ",cp,"\n")
				#Error in if (nrow(sigpaths)) { : argument is of length zero
				cpat = sumset$path_summary_each_patient[[cp]]
				if(!is.null(cpat$pathsummary)){
					hwrite(paste("<br><br>In patient",cp,"(",ci,"of 5 patients with the largest number of",targetType,"genes):"),p,heading=2)
					#third: find genes in those pathways which are significantly aberrational
					sigpaths = cpat$pathsummary[cpat$pathsummary$hyperg_p_w_FDR<0.05,,drop=F]
					if(nrow(sigpaths)){
						hwrite(paste(targetType,"genes in pathways significantly enriched (hypergeometric, with BD FDR correction) for",targetType,"genes"),p,heading=2)
					}else{
						hwrite(paste("No paths were found to be significantly enriched for",targetType,"genes as to the hypergeometric test with BH FDR"),p,heading=2)
						sigpaths = cpat$pathsummary[cpat$pathsummary$hyperg_p_w_FDR<0.05,]
						if(nrow(sigpaths)){
							hwrite(paste(targetType,"genes in pathways significantly enriched (hypergeometric test) for",targetType,"genes"),p,heading=2)
						}else{
							hwrite(paste("No paths were found to be significantly enriched for",targetType,"genes as to BH FDR"),p,heading=2)
							hwrite(paste(targetType,"genes in top 5 pathways most enriched (proportion enrichment) for ",targetType,"genes"),p,heading=2)
							sigpaths = cpat$pathsummary[order(cpat$pathsummary[,5],decreasing=T),][1:5,]
						}
					}
					sigpaths2 = cbind(sigpaths,cpat$active_genes_ea_path[rownames(sigpaths),])
					
					colnames(sigpaths2) = gsub(pattern="_",replacement=" ",x=colnames(sigpaths2))
					colnames(sigpaths2)[ncol(sigpaths2)] = paste(targetType,"genes in this pathway")
					hwrite(sigpaths2, p, table.style= 'border-collapse:collapse', table.cellpadding='1px', row.names=F)
					
					hwrite(paste(nrow(sigpaths2),"pathways total."),p,heading=2)
					
					hwrite("<br>",p)
					
					#find commonality of genes in top paths
					tmp = strsplit(cpat$active_genes_ea_path[rownames(sigpaths),],split=" ")
					names(tmp) = rownames(sigpaths)
					
					tmpmat = list_to_table(tmp)
					
					gcounts = matrix(rep(T,nrow(tmpmat))%*%tmpmat,ncol=1)
					rownames(gcounts) = colnames(tmpmat)
					hwrite(paste("Number of the above",nrow(tmpmat),"pathways each gene is found",targetType,"in."),p,heading=2)
					colnames(gcounts) = paste("Number of the above pathways the gene is found in")
					hwrite(gcounts, p, table.style= 'border-collapse:collapse', table.cellpadding='1px', row.names=T)
					hwrite(paste(ncol(tmpmat),"genes total."),p,heading=2)
				}
			}
		}
		
	}
	
	print("About to write image tags... ")
	print(names(sumset))
	#this is defined above: imageDir = paste(baseDir, "/imageSlots/")
	#search the summary for imageSlots
	if("imageSlots"%in%names(sumset)){
		if(length(sumset$imageSlots)){
			print("looking at image slots")
			print(fname)
			hwrite(x=headings[7], page=p, heading=2)
			hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
			print(sumset$imageSlots)
			images = as.matrix(sumset$imageSlots$imageSlots, ncol=1)
			print(images)
			if(nrow(images)){
				for(i in 1:nrow(images)){
					iname = rownames(images)[i]
					print(iname)
					ifile = images[i,1]
					hwrite(iname,p,heading=2)
					tmpImageLoc = paste("./imageSlots/",ifile, sep="")
					hwriteImage(tmpImageLoc, p, br=TRUE)
				}
			}
		}else{
			hwrite(x="Sorry, no images were found", page=p, heading=2)
		}
		
	}
	
	closePage(p)
	cat("\nSaved summary to HTML file: ",fname,"\n")
	return()
}# htmlSummary()


# htmlSummary(sumset=somatic_summary,
# 						fname="./quick_check_aml_somatic_nci.html")
#htmlSummary(sumset=somatic_summary, fname="./output/testSavesomaticsummary.html") 

#toHTML outputs an HTML file
#the input table_list cat have: 
#															tables
#																			
#															image file names
#
toHTML<-function(table_list,
								 fname="./test.html",
								 plimit=.05,
								 maxrows=100,
								 reorder=NULL, 
								 pagetitle=NULL, 
								 limit_col=NULL, 
								 path_detail=NULL){
	cat("____________________________Saving HTML page to :",fname,"\n")
	startHTMLPlugIns()
	fnameroot = paste(strsplit(fname, split="/")[[1]][1:(length(strsplit(fname, split="/")[[1]])-1)],sep="/",collapse="/")
	print(fnameroot)
	print(fname)
	
	p = openPage(fname)
	
	if(!is.null(pagetitle)) hwrite(pagetitle, p, heading=1)
	
	sl1 = "<A href=\"#section"
	sl2 = "\">"
	sl3 = "</A><BR>"
	hwrite("<A name=\"section0.1\">Table of contents</A>", p, heading=2)
	#first, print the list of elements
	hwrite("These are the sections of this HTML page:", p, heading=3)
	for(i in 1:length(table_list)){
		curlink = paste(sl1,i,sl2, names(table_list)[i], sl3, sep="")
		hwrite(curlink, p, heading=3)
	}
	trg1 = "<A name=\"section"
	trg2 = "\">"
	trg3 = "</A>"
	if(length(table_list)){
		
		if("settings"%in%names(table_list)){
			
			#if it's from the overlap analysis
			if(class(table_list$settings)=="list"){
				if(sum( names(table_list$settings)%in%c("basic_overlap", "defaultSummaryTable", "combinedAberrations") ) ){
					nsnames = paste0("setting_for_",names(table_list$settings))
					names(table_list$settings)<-nsnames
					table_list[nsnames] = table_list$settings
					table_list$settings=NULL
				}
			}
		}
		
		for(i in 1:length(table_list)){
			if(!is.list(table_list)) print(table_list)
			curname = names(table_list)[i]
			if(VERBOSE) print(curname)
			if(is.null(curname)) print(table_list)
			if(VERBOSE) print("names(table_list)")
			if(VERBOSE) print(names(table_list))
			if(VERBOSE) print(i)
			if(VERBOSE) print("stops here boo hoo")
			postfix = strsplit(curname,split="\\." )[[1]][ length(strsplit(curname,split="\\.")[[1]]) ]
			
			if( grepl(pattern="setting", x=curname) & class(table_list[[i]])!="data.frame" ){
				table_list[[i]] = listToDf(lst=table_list[[i]])
			}
			
			if( postfix %in% c("png","jpeg","bmp")){#if it's an image
				curnameanchor = paste(trg1, i, trg2, curname, trg3, sep="")
				if(VERBOSE) print(curnameanchor)
				hwrite(curnameanchor, p, heading=2)
				hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
				
				cat("\nImage found")
				hwriteImage(curname, p, br=TRUE)
			}else if(curname=="imageSlots"){
				#add anchor
				curnameAnchor = paste(trg1, i, trg2, "Network diagrams of pathways which are significantly enriched in aberration and contain drug targets", trg3, sep="")
				if(VERBOSE) print(curnameAnchor)
				hwrite(curnameAnchor, p, heading=2)
				hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
				
				cat("Image list found...\n")
				cat("(file:",fname,")\n")
				cat("\nCurrent item being output:",curname,"\n")
				
				#if the curname is imageSlots, 
				if(VERBOSE) print(names(table_list))
				#get the matrix out and print the images
				imageSet = table_list[["imageSlots"]][["imageSlots"]]
				#check if there are images in there
				if(is.null(imageSet)){
					if(VERBOSE) print("no images found in image folder")
					hwrite("Sorry, images not available", p, heading=2)
				}else{
					if(VERBOSE) print("images found..")
					if(VERBOSE) print(imageSet)
					imageSetNames = names(imageSet)
					if(is.matrix(imageSet))  imageSetNames = rownames(imageSet)
					for(im in 1:length(imageSet)){
						imageName = imageSetNames[im]
						if(VERBOSE) print(imageName)
						imageFile = paste("./imageSlots/",imageSet[im], sep="")
						if(VERBOSE) print(imageFile)
						hwrite(imageName, p, heading=3)
						hwriteImage(imageFile, p, br=TRUE)
					}
				}
			}else if(grepl("_each_patient",x=curname, ignore.case=T)){#it's the ... for each patient
				print("OLA eap")
				eaPatSet = table_list[[i]]
				for(patsetname in names(eaPatSet)){
					print("**********************************************************")
					ePatRoot = paste(fnameroot,curname,patsetname,paste("overlap_analysis/"),sep="/")
					SaveToHTML(results=eaPatSet[[patsetname]], 
										 study_name=patsetname,
										 fileRoot=ePatRoot,
										 selNames=names(eaPatSet[[patsetname]]))
				}
			}else if(grepl("_summary",x=curname,ignore.case=T)){
				print("Found a summary:::::::::::::::::::::::::::+:::::::::::::::::::::::::")
				
				sumfname =paste(fnameroot,curname,paste(curname,"_.html",sep=""),sep="/")
				ctab = table_list[[curname]]
				
				htmlSummary(path_detail=path_detail, 
										sumset=ctab,
										fname=sumfname, 
										pagetitle=curname)
			}else if(is.list(table_list[[i]])&!is.data.frame(table_list[[i]])){
				print("is.list(table_list[[i]])&!is.data.frame(table_list[[i]])")
				#if it's a list, then call toHTML again
				listfname = paste(fnameroot,curname,sep="/")
				toHTML(path_detail=path_detail, 
							 table_list=table_list[[i]], 
							 limit_col=limit_col,
							 pagetitle=paste(pagetitle, curname,sep=":"),
							 reorder=reorder,
							 fname=listfname, 
							 plimit=plimit, 
							 maxrows=maxrows)
				
			}else{ # it's not an image, summary table, folder of images or 'analysis each patient' set
				curname = paste(trg1, i, trg2, curname, trg3, sep="")
				hwrite(curname, p, heading=2)
				
				hwrite("<A href=\"#section0.1\">Back to top</A><BR>", p, heading=4)
				
				initial_table = table_list[[i]]#if it's turned into a vector clean table will not work; needs to be data frame
				if(is.vector(initial_table)&!is.list(initial_table)){
					initial_table = as.matrix(initial_table,ncol=1)
				}
				print(dim(initial_table))
				if(nrow(initial_table)){
					curtab = cleanTables(tab=initial_table)#makes sure there are only 3 sig figs instead of 9
					print(dim(initial_table))
					if(!is.null(limit_col)&sum(limit_col%in%colnames(curtab))){#if there are columns to limit the output by
						cat("There are columns by which to limit the number of table rows printed.\n")
						for(c in limit_col){#for each type of limit that is to be set, scan the table for the limit's column,
							if(c%in%colnames(curtab)){#if the column is found, limit and sort the table based on that column
								curindex = as.double(curtab[,c])<plimit
								tmptab=curtab[curindex,]
								curorder = 1:length(tmptab)
								if(!is.null(reorder)){
									print("reorder")
									curorder = order(tmptab[,c],decreasing=reorder)
								}
								colnames(tmptab)<-gsub(pattern="_",replacement=" ",x=colnames(tmptab))
								hwrite(paste(curname,"limited to ",c,"<",as.character(plimit),"(", as.character(nrow(tmptab)), "records exist, showing top", as.character(min(maxrows,nrow(tmptab))) ," )"), p, heading=2)
								
								hwrite(x=tmptab[curorder[1:min(maxrows,nrow(tmptab))],], page=p, 
											 row.names=!sum(rownames(tmptab)==tmptab[,1]),
											 table.style= 'border-collapse:collapse', 
											 table.cellpadding='1px')
								print("header and table written")
							}
						}
					}else{ #if there are no columns to limit the output by
						
						curtab = curtab[1:min(maxrows,nrow(curtab)),,drop=F]#limit the number of rows
						if(!is.null(reorder)){#order the rows if wanted
							curtab = curtab[order(curtab[,5],decreasing=T),,drop=F]
						}
						title = paste("(", as.character(nrow(initial_table)), "records, showing top", as.character(nrow(curtab)) ,")")
						hwrite(title, p, heading=4)
						colnames(curtab)<-gsub(pattern="_",replacement=" ",x=colnames(curtab))
						hasRowNames = !sum(rownames(curtab)==curtab[,1])
						hwrite(curtab,row.names=hasRowNames, 
									 p, 
									 table.style= 'border-collapse:collapse', 
									 table.cellpadding='1px')
						#readline("press enter to continue")
						cat("~")
					}
				}else{
					hwrite(paste("Zero",curname,"were found."), p, heading=4)
				}#if/else there are rows in curtab
			}#if/else is image
		}#for each element in the table_list
	}
	
	cat("\nHTML summary written to",fname,"\n")
	closePage(p)
}


SaveToHTML<-function(results,
										 study_name, 
										 selNames=NULL,
										 fileRoot=NULL,
										 path_detail=NULL){
	print("Saving summary to HTML..")
	
	if(!grepl(pattern="^test", ignore.case=T, x=study_name)){
		study_name = paste0("study_", study_name)
	}
	
	if(is.null(fileRoot)) fileRoot = paste("./output/",study_name,"/results",sep="")
	print("file root:")
	print(fileRoot)
	
	isels = names(results)
	badi = grep(pattern="file", x=names(results), ignore.case=T)
	if(length(badi)){
		isels = names(results)[-badi]
	}
	#pick out all the indexes except the path_file_name
	
	if(is.null(selNames)){
		#let the user decide
		selvector.li = selectionList(valcol=isels)	
		#pullOutTheCorrectNames
		selNames = isels[selvector.li]
	}
	
	for(n in selNames){
		print("selnames:")
		print(selNames)
		print(n)
		#save em
		ssum = results[[n]]
		if(grepl(pattern="summary", x=n, ignore.case=T)){
			print("chose summary")
			root2 = paste(fileRoot,"/",n,"/",sep="")
			sumfname = paste(root2,n,"_ViewableSummary.html",sep="")
			# 			readline("enter to continue")
			htmlSummary(path_detail=path_detail, sumset=ssum, fname=sumfname, pagetitle=n)
			htmlPerPatient(results, summ=n, fileroot=root2, overlap=F, path_detail=path_detail)
		}else if(grepl(pattern="overlap_analysis$", x=n, ignore.case=T)){
			print("chose overlap")
			
			root2 = paste(fileRoot,"/",n,"/",sep="")
			dir.create(root2,recursive=T,showWarnings=F)
			outHTMLname = paste(root2, "overlap_analysis.html",sep="")
			print(outHTMLname)
			print(n)
			# 			readline("enter to continue")
			toHTML(path_detail=path_detail, 
						 table_list=results$overlap_analysis, 
						 fname=outHTMLname, 
						 maxrows=10000, 
						 pagetitle=n)
		}else if(grepl(pattern="overlap_analysis_each_patient$", x=n, ignore.case=T)){
			
			root2 = paste(fileRoot,"/", n, "/", sep="")
			for(pn in names(summ)){
				print(pn)
				print(root2)
				root3 =paste(root2, pn, "/overlap_analysis/", sep="")
				print(root3)
				dir.create(root3, recursive=T, showWarnings=F)
				outHTMLname = paste(root3, "overlap_analysis.html",sep="")
				print(outHTMLname)
				# 				readline("enter to continue")
				#for each name in overlap analysis each patient
				#pull the patient
				curpat = summ[[p]]
				toHTML(path_detail=path_detail, table_list=curpat, fname=outHTMLname, maxrows=10000, pagetitle=n)
			}
		}
	}
}

htmlPerPatient<-function(results, summ, fileroot, overlap=T, path_detail=NULL){
	cat("\nWriting an HTML summary for each patient...........\n")
	fileroot = paste(fileroot,"path_summary_each_patient/",sep="")
	psep = results[[summ]]$path_summary_each_patient
	for(patient in names(psep)){
		cat("\n***************************",patient,"*******************************\n")
		#for each patient, establish a summary
		perPatRoot = paste(fileroot,patient,"/",summ,"_ViewableSummary.html",sep="")
		print(summ)
		ss = psep[[patient]]
		if(overlap) ss = psep[[patient]][[summ]]
		htmlSummary(path_detail=path_detail, sumset=ss, fname=perPatRoot)
	}
}

#twoHistOnePlot
#puts two histograms on the same plot, shows each with a different color, and overlap with a third color
twoHistOnePlot<-function(dataset1,dataset2,frequency_not_density=T,main_title="main_title not set",x_label="x_label not set", y_label="Frequency",legend_titles=c("legend_titles two item vector","not included"),breaks=c(30,30)){
	
	pmin = min(c(dataset1, dataset2))
	pmax = max(c(dataset1, dataset2))
	
	p1 <- hist(dataset1, freq=frequency_not_density,breaks=breaks[1],)                     
	p2 <- hist(dataset2, freq=frequency_not_density,breaks=breaks[2])
	
	hist = plot( p2, col=rgb(0,0,1,1/4), xlim=c(pmin,pmax), 
							 freq=frequency_not_density, 
							 main=main_title, 
							 xlab=x_label, 
							 ylab = y_label)  # first histogram
	hist2 = plot( p1, col=rgb(1,0,0,1/4), xlim=c(pmin,pmax), 
								add=T,
								freq=frequency_not_density)  # second
	
	if(!is.null(legend_titles)){
		legend("topright", inset=.05,
					 legend_titles, 
					 fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), 
					 horiz=F)
	}
	return(list(p1=p1, p2=p2))
}

scatterhist = function(x, y, xlab="", ylab="", main=""){
	zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
	layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
	xhist = hist(x, plot=FALSE)
	yhist = hist(y, plot=FALSE)
	top = max(c(xhist$counts, yhist$counts))
	par(mar=c(3,3,1,1))
	plot(x,y, main="")
	par(mar=c(0,3,1,1))
	barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
	par(mar=c(3,0,1,1))
	barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
	par(oma=c(3,3,0,0))
	
	mtext(text=main,side=3,outer=T,padj=2)
	mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
				at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
	mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
				at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
	par(mfrow = c(1,1))
	par(mar=c(5, 4, 4, 2) + 0.1 )
}

scatterhist2 = function(x1, y1, x2, y2, xlab="", ylab="", main="", legendTxt=NULL){
	
	zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
	layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
	
	xhist = twoHistOnePlot(dataset1=x1, dataset2=x2, 
												 legend_titles=NULL,
												 main_title="", 
												 y_label="",
												 x_label="")#hist(x, plot=FALSE)
	
	yhist = twoHistOnePlot(dataset1=y1, dataset2=y2,
												 legend_titles=NULL,
												 main_title="", 
												 y_label="",
												 x_label="")#hist(y, plot=FALSE)
	
	top = max(c(xhist$p1$density, yhist$p1$density, xhist$p2$density, yhist$p2$density))
	par(mar=c(3,3,1,1))
	
	plot(x1,y1, main=main, col=rgb(0,0,1,1/4), pch=46)
	points(x2,y2,col=rgb(1,0,0,1/4), pch=46)
	
	
	if(!is.null(legendTxt)){
		legend("topright", c(legendTxt[1], legendTxt[2]), fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)))
	}
	
	par(mar=c(0,3,1,1))
	
	barplot(xhist$p1$density, axes=FALSE, ylim=c(0, top), space=0, col=rgb(0,0,1,1/4))
	barplot(xhist$p2$density, axes=FALSE, ylim=c(0, top), space=0, add=T, col=rgb(1,0,0,1/4))
	par(mar=c(3,0,1,1))
	
	barplot(yhist$p1$density, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col=rgb(0,0,1,1/4))
	barplot(yhist$p2$density, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, add=T, col=rgb(1,0,0,1/4))
	par(oma=c(3,3,0,0))
	
	mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
				at=.8 * (mean(x1) - min(x1))/(max(x1)-min(x1)))
	mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
				at=(.8 * (mean(y1) - min(y1))/(max(y1) - min(y1))))
	par(mfrow = c(1,1))
	par(mar=c(5, 4, 4, 2) + 0.1 )
	
}

addGeneInfo<-function(genelist, path_detail){
	#adds additional gene info to a list of genes
	#indicates genes which did not have hugo records
	#initilize variables
	htab = path_detail$HUGOtable
	
	# 	#assure there are no symbols given as "\\N" 
	# 	genelist = genelist[genelist!="\\N"]
	
	#attempt to make sure the symbols match 
	genelist = corsym(genelist, hugoref=htab, verbose=T)
	
	#extract the needed columns
	xhtab = htab[,c("Approved.Symbol", "Approved.Name", "Locus.Type", "Chromosome", "Previous.Names")]
	rownames(xhtab) = xhtab$Approved.Symbol
	
	#figure out which gene symbols we actually have info on:
	have = genelist%in%xhtab$Approved.Symbol
	
	#reduce the gene list to those we have info on and report those we dont have info on
	goodset = genelist[have]
	badset = genelist[!have]
	
	out=list()
	out$'Genes and their records' = xhtab[goodset,]
	out$'Gene symbols for which records could not be found' = badset
	
	return(out)
}

stackedGeneBar<-function(tcga_som){
	#stackedGeneBar
	#takes: tcga_som: initial input .maf table, after cleaning of repeat rows and gene symbols
	#outputs: stacked bar plot
	ufilt = tcga_som
	gs = summarize_by(col=tcga_som[,"Hugo_Symbol"], display=F)
	top20 = head(gs[order(gs[,2],decreasing=T),],20)
	gsnames = top20$types
	
	toprows = ufilt[ufilt$Hugo_Symbol%in%gsnames,]
	
	slimrows = toprows[,c("Hugo_Symbol","Variant_Classification")]
	
	#set up the output matrix
	mtypes = unique(slimrows$Variant_Classification)
	mmat = matrix(data=0,ncol=length(mtypes), nrow=length(gsnames), dimnames=list(gsnames, mtypes))
	
	for(n in gsnames){
		#get the rows for those names
		msum = summarize_by(col=slimrows[slimrows$Hugo_Symbol==n,"Variant_Classification"], display=F)
		mmat[n,msum$types] = msum$counts
	}
	mdf  = cbind.data.frame(gene = rownames(mmat),mmat)
	
	par(las=2) # make label text perpendicular to axis
	par(mar=c(4,8,3,3)) # increase y-axis margin.
	oldmar =c(5,4,4,2)
	oldlas = 0
	
	scol = 2
	colseq = (scol:(scol + ncol(mmat)))*10 + 3
	barplot(t(mmat), 
					legend.text = colnames(mmat),
					xlab="Number of mutations found in gene",
					main="Mutation types for the top 20 most mutated genes", 
					col=colors(distinct=T)[colseq],
					horiz=TRUE, 
					cex.names=0.7)
	par(las=oldlas)
	par(mar=oldmar)
	#barplot(mmat, horiz=T)
}

OneColStackedBar<-function(incol=tcga_som$Variant_Classification, mainText = "Counts of mutation type across cohort"){
	#makes stacked bar plot with one column from one row
	#takes: one column of data
	#outputs: stacked bar plot
	
	gs = summarize_by(col=incol, display=F)
	mmat = gs[,"counts",drop=F]
	rownames(mmat) = gs$types
	colnames(mmat)<-""
	mmat = t(mmat)
	scol = 2
	colseq = (scol:(scol + ncol(mmat)))*10 + 3
	barplot(t(mmat), width=2, ylab="Number of mutations", args.legend = list(inset=c(.05,.02), x="bottomright", bg="white"),
					legend = colnames(mmat),
					main=mainText, 
					col=colors(distinct=T)[colseq],
					horiz=F, 
					cex.names=0.8)
	
}

