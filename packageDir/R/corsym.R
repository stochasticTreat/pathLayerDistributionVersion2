


getOfficialGeneIds<-function(idSource){
	
	tab =NULL
	if(class(idSource)=="Path_Detail"){
		tab = idSource$symtable
	}else if( class(idSource)=="Study" ){
		tab = idSource@studyMetaData@paths$symtable
	}
	
	return(tab[,"Approved.Symbol"])
}


#'@title Correct unofficial gene symbols. 
#'@description Function to correct and unify gene identifiers across a study arms. 
#'@param symbol_set the set of gene identifiers to be corrected
#'@param symref The gene identifier lookup table containing a set of official gene symbols. If this contains two columns, it should be in the format of the hgnc.table data set provided by the CRAN package, HGNChelper. Else, this should be a whole data table of gene symbol look ups as provided by genenames.org
#'@param verbose If set to TRUE, this flag will cause the gene symbol corrections to be conducted automatically. 
#'@param col2 If two columns of data are provided as the symbol_set input to corsym, col2 should provide the name of a second column containing additional data bout the gene identifier to be checked, such as a chromosomal location.
#'@param correctionsfile character string giving the file path to the a gene symbol corrections file. This file should contain two columns: old_symbol and new_symbol containing the errant and the correct symbols, respectively. This is the primary set of symbols to be used in coordinating gene symbol corrections between the multiple arms of the study. 
#'@return The character vector of gene identifiers provided as an input, with any possible corrections made. 
#'@details This function allows a user to interactively correct gene identifiers so they can be coordinated between study arms. All corrections can optionally be recorded to the corrections file so that at latter time they can be run automatically, without user interaction, and so that they can be automatically re-used in correcting gene symbols from other study arms. 
#'@export
#'@examples
#'toCheck = c("p53", "FLT3", "ASM1", "ASML3B","ARF","APOBEC3C")
#'corrected = corsym(symbol_set=toCheck) #using HGNChelper and previously made corrections
#'
#' symref  = getHugoSymbols()
#'
#'\dontrun{
#'symref  = getHugoSymbols(curhugofname="./reference_data/current_hugo_table_slim.txt")
#'corrected = corsym(symbol_set=toCheck, symref=symref, verbose=T)
#'}
corsym<-function(symbol_set, 
								 symref=NULL, 
								 verbose=T, 
								 col2="Chrom", 
								 correctionsfile="./reference_data/gene_symbol_corrections_list.txt"){

	isHugo = T
	hgnchelper = F
	symtab = NULL
	
	if(is.null(symref)){
		data("hgnc.table", package="HGNChelper", envir=environment())
		symref = hgnc.table
	}
	
	#option 0: no symbol correction available for current data type
	if(class(symref)=="Study"){
		symtab = symref@studyMetaData@paths$symtable
		isHugo = "HUGO"==symref@studyMetaData@paths$symbol_type
	}else if(class(symref)=="Path_Detail"){
		symtab=symref$symtable
		isHugo = "HUGO"==symref$symbol_type
	}else if(class(symref)=="data.frame"){
		symtab = symref
	}
	
	if(!isHugo){
		print("No corrections available for this symbol type")
		return(symbol_set)
	}

	#check if hgnchelper should be used
	#use it if the symtab has two columns, symbol and approved.symbol, and if the symbol type is HUGO
	hgnchelper = isHugo & ncol(symtab)==2 & !sum( !c("Symbol","Approved.Symbol")%in%colnames(symtab) ) 
	
	if(hgnchelper){ 	#option 1: use HGNChelper package
		print("Correcting by HGNChelper")
		res = correctByHgncHelper(symbol_set=symbol_set, symtab=symtab, correctionsfile=correctionsfile)
		
	}else{	#option 3: use corsym_full()
		print("Correcting with the advanced symbol correction script")
		res = corsym_full(symbol_set=symbol_set, 
											symref=symref, 
											verbose=verbose, 
											col2=col2, 
											correctionsfile=correctionsfile)
	}
	
	if(!is.null(dim(res))){#if the results have more than one dimension, reduce to only a vector
		res = as.vector(res[,1])
	}
	return(res)
}


checkPreviousCorrections<-function(bsub, ctab){

	pcori = bsub%in%ctab$old_symbol
	cat(sum(pcori),"symbols could be corrected using previous symbol corrections.\n")
	if(sum(pcori)){
		#make a dictionary
		tmp = ctab[ctab$old_symbol%in%bsub,"new_symbol"]
		names(tmp)<-ctab[ctab$old_symbol%in%bsub,"old_symbol"]
		#swap out the symbols
		bsub[pcori] = tmp[bsub[pcori]]
	}
	return(bsub)
}


correctFromHelperPrevious<-function(bsub, verbose){
	
	cres = checkGeneSymbols(x=bsub)
	
	fixi = !is.na(cres$Suggested.Symbol)
	multi = grepl(pattern=" /// ",x=cres$Suggested.Symbol)
	
	cat("Found",sum(multi),"symbol(s) that map to more than one approved symbol.\n")
	if(verbose&sum(multi)){
		if(readline("Would you like to attempt to fix multi-mapping symbols? (y/n) ")=="y"){
			
			for(i in which(multi)){
				cat("For the provided symbol,",cres$x[i],
						"\nthese possible approved symbols were found:\n",
						gsub(pattern="///",  replacement=" ", x=cres$Suggested.Symbol[i]),"\n")
				uin = readline("Please key in the gene symbol you would like to use: ")
				cres$Suggested.Symbol[i] = uin
				multi[i] = T
			}
			
		}
	}
	
	cat("Found appropriate corrections for",sum(fixi),"symbols.\n")
	#now correct the symbols in symbol_set
	#first make the dictionary
	tmp = cres$Suggested.Symbol[fixi&!multi]
	names(tmp)<-cres$x[fixi&!multi]
	#now swap the symbols
	bsub[fixi&!multi] = tmp[bsub[fixi&!multi]]
	cat("Corrected",sum(fixi&!multi),"symbols.\n")
	return(bsub)
}


#'@title Correct gene symbols using the HGNChelper package
#'@description Correct gene symbols using the HGNChelper package. Provides interactive or automatic correction of gene identifiers, and corrdination of any correction between multiple calls to this functions.  
#'@param symbol_set Vector, matrix or data.frame (if more than 2 dimensions, the symbols must be in the first column)
#'@param correctionsfile A file name for a file of gene symbol corrections. This file should be tab delimited and contain two named columns: old_symbol and new_symbol. The old_symbol column contains incorrect identifiers to be corrected, and the new_symbol column contains the corrections. There should be no row names. 
#'@param symtab The table of official, orthodox gene identifiers. If not provided, the \code{hgnc.table data frame} object from HGNChelper will be used. 
#'@param verbose A logical flag indicating if interactive mode should be run. 
#'@return A character vector containing the identifiers for the set of genes provided in the symbol_set function argument, but with any applicable corrections made. 
#'@import HGNChelper 
#'@export correctByHgncHelper
correctByHgncHelper<-function(symbol_set, 
															correctionsfile="./reference_data/gene_symbol_corrections_list.txt", 
															symtab=NULL, 
															verbose=F){

	#1: symbol correct : do nothing
	#2: symbol incorrect, no correction available : report to user
	#3: symbol incorrect, simple correction available : correct
	#4: symbol incorrect, multiple correction options : allow user to correct or do nothing, (use the first?)
	
	if(!is.vector(symbol_set)) symbol_set  = symbol_set[,1,drop=T]
	#clean the gene symbols
	if(verbose){
		cat("\nTo conduct symbol comparrisons and corrections, these reformattings were made:")
		cat("\nConversion to upper case.\nRemoval of leading and trailing spaces.",
				"\nConversion of spaces and periods to dashes.\n",
				"***NOTE: These reformattings are not saved or recorded to any file!!\n")
	}
	symbol_set= cleanGeneSymbols(symbol_set)
	
	
	if(is.null(symtab)){
		if( !exists("hgnc.table") ) data("hgnc.table", package=HGNChelper)
	}else{
		hgnc.table = symtab
	}
	
	ctab = getSymbolCorrectionTable(correctionsfile=correctionsfile)

	#check which if the input are correct. 
	goodsymsi  = symbol_set%in%hgnc.table$Approved.Symbol
	if(sum(!goodsymsi)) cat("\n",sum(!goodsymsi), "gene symbols were not found in the approved HUGO symbols.\n")
	bsub = symbol_set[!goodsymsi]
	
	cat("Checking previously made corrections...\n")
	bsub = checkPreviousCorrections(bsub=bsub, ctab=ctab)
	
	if(verbose){
		cat("Checking for Excel-morgified symbols...\n")
		corbsub = findExcelGeneSymbols(x=bsub)
		
		cat("Checking previously official symbols...\n")
		corbsub = correctFromHelperPrevious(bsub=corbsub, verbose=verbose)
		
		checkAddToCorrectionsFile(bsub=bsub, 
															corbsub=corbsub, 
															ctab=ctab, 
															correctionsfile=correctionsfile, 
															verbose=verbose)
		bsub = corbsub
	}

	symbol_set[!goodsymsi] = bsub
	return(symbol_set)
}

checkAddToCorrectionsFile<-function(bsub, corbsub, ctab, correctionsfile,verbose){
	difi = bsub!=corbsub
	if(sum(difi)&verbose){
		if("y"==readline("New symbol corrections were made.\nWould you like to add them to the corrections file? (y/n) ")){
			corblock = cbind(bsub[difi], corbsub[difi])
			colnames(corblock)<-c("old_symbol","new_symbol")
			ctab = rbind(ctab, corblock)
			addCorrections(new_corrections=corblock, correctionsfile=correctionsfile, correction_set=ctab)
		}
	}
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
corsym_full<-function(symbol_set, symref=NULL, verbose=T, col2="Chrom", correctionsfile="./reference_data/gene_symbol_corrections_list.txt"){
	
	if(verbose) cat("\n\nChecking that gene symbols match official HUGO gene symbols. . . . . \n")
	curhugofname="./reference_data/current_hugo_table.txt"
	if(is.null(symref)){
		cref=getHugoSymbols(curhugofname=curhugofname)
	}else if(class(symref)=="data.frame"){#the hugo table was passed directly
		cref = symref
	}else if(is.character(symref)){#if symref is a character, it is a file name; open the file
		cref=getHugoSymbols(curhugofname=symref)#the direct passing of a hugo file name.  . .not sure if we're doing that any more...
	}else if(class(symref)=="Study"){
		study=symref
		if(study@studyMetaData@geneIdentifierType!="HUGO"){
			message("Checking and correction of non-HUGO symbols not implemented.")
			return(symbol_set)
		} 
		cref=symref@studyMetaData@paths$symtable
		symref=cref
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
	resave = F
	###################################  Check Corrections File
	raw_correction_set = getSymbolCorrectionTable(correctionsfile=correctionsfile)
	if(file.exists(correctionsfile)){
		
		raw_correction_set_tmp=corListCheck(cl=raw_correction_set, symtab=symref)
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
			cat("\nnot approved 1.2:\n")
			print(not_approved)
			if(resave){
				dir.create(path=dirname(path=correctionsfile), showWarnings=F, recursive=T)
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
				switches = checkPreviousSymbols(symbols=not_approved, indexes=1:nrow(not_approved), symlookup=cref, col2=col2)
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
				
				syncor = checkSynonyms(symbols=not_approved, indexes=1:nrow(not_approved), symlookup=cref, col2=col2)
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
				cat("\nThese are the new corrections that were just added:\n")
				print(new_corrections)
				if(readline("Would you like to save this set of corrections just made to the corrections file (y/n)")=="y"){
					addCorrections(new_corrections=new_corrections, correction_set=correction_set, correctionsfile=correctionsfile)
				}
			}
		}
	}#if verbose
	return(symbol_set[,"Hugo_Symbol"])
}#corsym function

getSymbolCorrectionTable<-function(correctionsfile){
	
	if(!file.exists(correctionsfile)){
		message("\nNote: cannot find gene identifier correction file.\nCopying the default correction file from package data.\n")
		dir.create(path=dirname(path=correctionsfile), recursive=T, showWarnings=F)
		file.copy(from=system.file("extdata/gene_symbol_corrections_list.txt", package = "packageDir"),
							to=correctionsfile)
	}
	
	raw_correction_set = read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")

	return(raw_correction_set)
}

#'@title Function to add corrections to the corrections file. Internal, used by corsym to coordinate symbol corrections. 
#'@description Adds corrections to the corrections file
#'@param new_corrections a two column matrix; column 1 = old, incorrect symbols, column 2 = new, corrected symbols
#'@param correctionsfile the file name of the corrections file
#'@param correction_set the original contents of the corrections file before new corrections were made
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
	dir.create(path=dirname(path=correctionsfile), showWarnings=F, recursive=T)
	write.table(x=final_corrections, file=correctionsfile, 
							quote=F, sep="\t", row.names=F, col.names=c("old_symbol", "new_symbol"))
	cat("\nSymbol corrections recorded in file:", correctionsfile, "\n")	
}



#examineHugoSet
#get summary information on a set of HUGO symbols
#takes: symbol_set: the list of symbols
#study_name and data_type describe the source of the symbols, for display on the read out
#
examineHugoSet<-function(symbol_set, study_name, data_type, curhugofname="./reference_data/current_hugo_table.txt"){
	
	cref = read.delim(file=curhugofname, header=T, stringsAsFactors=F, na.strings="-")
	approvedHugoFile = paste("./output/",study_name,"_approved_hugo_w_annotation_from_",data_type,"_data.txt",sep="")
	
	cat("\nA table of approved HUGO symbols, including those just corrected, can be found in the file",approvedHugoFile,"\n")
	
	in_study = cref[cref[,"Approved.Symbol"]%in%symbol_set,]
	not_hugo = setdiff(x=symbol_set, y=cref[,"Approved.Symbol"])
	nsymtab = matrix(nrow=length(not_hugo),ncol=ncol(in_study))
	colnames(nsymtab)<-colnames(in_study)
	nsymtab[,"Approved.Symbol"] = not_hugo
	nsymtab[,"Locus.Type"] = rep("Not HUGO symbol",times=length(not_hugo))
	nsymtab[,"Locus.Group"] = rep("Not HUGO symbol",times=length(not_hugo))	
	
	in_study = rbind(in_study, nsymtab)
	
	write.table(x=in_study,file=approvedHugoFile,quote=F,sep="\t",row.names=F,col.names=T)
	
	return(in_study)
}


#takes: indexes: indexes of symbols to be checked
#       symbols: list of symbols which the indexes refer to 
#                             (those to be switched and those to remain)
#       symlookup: hugo look up table with columns : <Approved.Symbols>, <Previous.Symbols>
#returns: table: <index to be switched> <what it should be switched to>
checkPreviousSymbols<-function(symbols, indexes, symlookup, col2){
	
	symlookup[,"Previous.Symbols"] = toupper(symlookup[,"Previous.Symbols"])
	
	oldsyms = NULL #output set of those which are previous symbols
	ssymbols = NULL #
	for(i in 1:length(indexes))
	{
		cur = symbols[indexes[i],"Hugo_Symbol"]
		curchrom = symbols[indexes[i], col2]
		cursym = gsub(pattern="[[:punct:]]", replacement=" ", x=cur)
		
		#grep the row against the previous symbols column
		pat = paste("(^|[[:blank:]])", cursym, "(,|$)", collapse="", sep="")
		ind = grep(pattern=pat, x=symlookup[["Previous.Symbols"]],ignore.case=T)
		if(length(ind)==1){#if grep caught something
			ssymbols = c(ssymbols, as.character(symlookup[ind,"Approved.Symbol"]))
			oldsyms = c(oldsyms, cur)
		}else if(length(ind)>1){#if grep found more than one matching previous symbol
			dcols = colnames(symlookup)[grep("date", colnames(symlookup), ignore.case=T)]
			dcols = c(dcols,colnames(symlookup)[grep(col2, colnames(symlookup), ignore.case=T)])
			cat("\nCorrecting gene symbol", as.character(i), "of", as.character(length(indexes)), "symbols that need correction.\n")
			cat("\nMultiple previous symbols were found to match", as.character(cursym), " (from Chrom ", as.character(curchrom), "). \n")
			cat("\nThese are the symbols that were found to match:\n")
			print(symlookup[ind,c("Approved.Symbol", "Previous.Symbols", dcols)])
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
#       symlookup: hugo look up table with columns : <Approved.Symbols>, <Previous.Symbols>
#returns: table: <index to be switched> <what it should be switched to>
checkSynonyms<-function(symbols, indexes, symlookup, col2){
	sindexes = NULL #output set of those which are previous symbols
	ssymbols = NULL
	if(!"Synonyms"%in%colnames(symlookup)){ #if the curr
		return(data.frame())
	}
	
	for(i in 1:length(indexes))
	{
		osym = symbols[indexes[i],"Hugo_Symbol"]
		curchrom = symbols[indexes[i], col2]
		cursym = gsub(pattern="[[:punct:]]", replacement=" ", x=osym)
		sterm=cursym
		acc = NULL
		while(T){
			#print("start while(T)")
			#grep the row against the previous symbols column
			pat = paste("(^|[[:blank:]])", sterm, "(,|$)", collapse="", sep="")
			res1 = grep(pattern=pat, x=symlookup[["Synonyms"]])
			res2 = grep(pattern=sterm, x=symlookup[["Synonyms"]], ignore.case=T)
			res3 = grep(pattern=pat, x=symlookup[["Approved.Symbol"]], ignore.case=T)
			cat("\nSearching synonyms for", cursym,"###################################################### \n")
			cat("#####in the data being processed, this symbol is associated with chromosome", curchrom, "\n\n########## Original query:", osym,
					"\n########## Current search term: ", sterm, "\n")
			if(length(res2)>0){
				cat("\n ############ These near matches were found: \n")
				print(symlookup[res2, c("Approved.Symbol","Approved.Name","Status", "Synonyms", 
																colnames(symlookup)[grep("date", colnames(symlookup), ignore.case=T)], 
																colnames(symlookup)[grep(col2, colnames(symlookup), ignore.case=T)])])
			}else{cat("\nNo synonyms were found for the search term.\n")}
			if(length(res3)==1){
				cat("\n ************ This exactly matching approved HUGO symbol was found: \n")
				print(symlookup[res3, c("Approved.Symbol","Approved.Name", 
																colnames(symlookup)[grep(col2, colnames(symlookup), ignore.case=T)])])
				if(readline(prompt="If you would like to accept this exactly matching approved HUGO symbol, press ENTER.")==""){
					acc = as.character(symlookup[res3, "Approved.Symbol"]) 
					break
				}
			}
			if(length(res1)==1){
				cat("\n ************ This exactly matching synonym was found: \n")
				print(symlookup[res1, c("Approved.Symbol","Approved.Name","Status", "Synonyms", 
																colnames(symlookup)[grep("date.approved", colnames(symlookup), ignore.case=T)], 
																colnames(symlookup)[grep(col2, colnames(symlookup), ignore.case=T)])])
				if(""==readline("If you would like to accept the exactly matching synonym above, please press ENTER\nIf you enter anything else, more options will be provided.")){
					acc = as.character(symlookup[res1, "Approved.Symbol"]) 
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
		#print("Out of while loop")
		if(acc != ""){
			ssymbols = c(ssymbols, acc)
			sindexes = c(sindexes, osym)
		}
	}#for
	out = cbind.data.frame(sindexes, ssymbols, stringsAsFactors=F)
	if(nrow(out)) colnames(out)<-c("old_symbol","new_symbol")
	#print("about to return")
	return(out)
}#checkSynonyms


#corListCheck
#if the HUGO symbols reference file has changed, 
#corListCheck will assure the symbol changes are consistent with the new reference file
#takes:		cl:		corrections list
#					symtab:	hugo reference table 
#returns: corrections list data frame
#
#works with files: "./reference_data/gene_symbol_corrections_list.txt"
#
corListCheck<-function(cl=NULL, symtab=NULL, curhugofname = "./reference_data/current_hugo_table.txt"){
	cat("\nScreening symbol correction table...\n")
	
	correctionsfile="./reference_data/gene_symbol_corrections_list.txt"
	if(is.null(cl)){
		cl = read.delim(file=correctionsfile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
	}
	
	#handle three cases:
	#1: symbol added to HUGO ref
	#2: symbol removed from HUGO ref
	#3: symbol changed
	if("Status"%in%colnames(symtab)){
		#for 2 and 3: 
		#check if any of the "new" are withdrawn:
		wnew = paste(cl$new_symbol, "~withdrawn", sep="")
		
		iswithdrawn = wnew%in%symtab$Approved.Symbol
		if(sum(iswithdrawn)){#if some of the symbols have been withdrawn, see if replacement symbols can be found
			cltmp = cbind.data.frame(cl[iswithdrawn,], wnew[iswithdrawn], stringsAsFactors=F)
			colnames(cltmp)[3] = "withdrawn"
			#for each row in cltmp, merg it's corresponding row in symtab
			symtabtmp = symtab[symtab$Approved.Symbol%in%cltmp$withdrawn,1:8]#extract the rows from symtab
			chmerge = merge(x=cltmp, y=symtabtmp, by.x="withdrawn", by.y="Approved.Symbol")
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
				cl = corListCheck(cl=cl,symtab=symtab)#recursive call, 
				#handles case where symbol was changed more than once
			}
		}
	}
	
	### now scan the old_symbols for appoved symbols
	#remove rows where the old_symbol is an approved symbol
	if(sum(cl$old_symbol%in%symtab$Approved.Symbol)){
		cat("\nSome rows in the symbol correction table were ",
				"found to correct symbols that were already approved.",
				"These rows will be removed to prevent data inconsistencies\n")
		cl = cl[!cl$old_symbol%in%symtab$Approved.Symbol,]
	}
	return(cl)
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


#getHugoSymbols()
#paths_detail: paths object: if passed, this function will excise the currently used hugo symbol set and return them.<them?>
#curhugofname: the name of the 
#returns HUGO lookup table
#
#function allows re-download of hugo cross ref file. 

#'@title getHugoSymbols
#'@description obtains HGNC/HUGO symbol look up table containing official set of HUGO symbols.
#'@param paths_detail A \code{Path_Detail} or a \code{Study} object. If a \code{Path_Detail} object is provided, the symbols will be extracted from its \code{symtable} slot. If a \code{Study} object is provided for this argument, the symbol lookup table will be exracted from the \code{Path_Detail} it contains.
#'@param curhugofname If interactive gene symbol correction is to be used, this argument should be the file path to the HUGO table as downloaded from genenames.org.
#'@param verbose Controlls if symbol corrections are to be interactive (if yes, curhugofname file must be supplied as it contains critical information, such as gene symbol status, past identifiers and synonyms)
#'@return Table of symbols: either a two column data.frame, the hgnc.table provided by HGNChelper, or a data frame as provided by genenames.org. Either of which have a column titled Approved.Symbol which contains official, approved symbols.  
#'@import HGNChelper
#'@export
#'@examples
#'#get the default HUGO lookup table, hgnc.table from HGNChelper.
#'hsyms = getHugoSymbols() 
#'#attempt to download and or use a full hugo lookup table from genenames.org
#'hsyms = getHugoSymbols(curhugofname="./reference_data/current_hugo_table_slim.txt") 
getHugoSymbols<-function(paths_detail=NULL, 
												 curhugofname=NULL,#"./reference_data/current_hugo_table_slim.txt",
												 verbose=F){
	# 	curhugofname = "./reference_data/current_hugo_table_slim.txt"

	if(!is.null(paths_detail)){
		if(class(paths_detail)=="Study"){
			cat("\nGetting Path_Detail object from Study object ..\n")
			paths_detail=paths_detail@studyMetaData@paths
		}
		cat("\nGetting HUGO symbols from Path_Detail reference calss object...\n")
		return(paths_detail$symtable)
	}#if/else

	if(is.null(curhugofname)){#if a file name for a hugo lookup table was not passed, return the hgnc.table from HGNChelper
		library("HGNChelper")
		data("hgnc.table", package="HGNChelper")
		if(verbose){
			if(readline("Would you like to use interactive gene symbol correciton? (enter y or n) ")=="n"){
				return(hgnc.table)
			}else{
				curhugofname = "./reference_data/current_hugo_table_slim.txt"
			}
		}else{
			return(hgnc.table)
		}
	}

	if(!file.exists(curhugofname)){
		downloadHugoLookupTable(curhugofname=curhugofname, autoDownLoad=!verbose)
	}
	
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
		cref = downloadHugoLookupTable(curhugofname = curhugofname)
	}
	return(cref)

}#getHugoSymbols


#'@title Obtain full HUGO gene symbol look up table. 
#'@description The function attempts to download the HUGO lookup table (includes previous symbols, synonyms and chromosome corrdinates) from the HGNC website genenames.org.
#'@param curhugofname \code{character string} The location the symbol lookup table should be saved. Default is the local directory at ./reference_data/current_hugo_table_slim.txt. For this program to use the full HUGO table, the file needs to be left here. 
#'@param autoDownLoad \code{logical} A flag indicating if the program should automatically attempt to download the table (T) or if it should just display information on how to manually download the table (F).
## '@examples
## 'curhugofname = "./reference_data/current_hugo_table_slim.txt"
## 'downloadHugoLookupTable(curhugofname=curhugofname, autoDownload=TRUE) #attempt to download the HUGO table
## 'hsyms = getHugoSymbols(curhugofname=curhugofname, verbose=F)
## 'file.remove(curhugofname)
downloadHugoLookupTable <- function (curhugofname="./reference_data/current_hugo_table_slim.txt", autoDownLoad=F) {
	
	if(!autoDownLoad){
		autoDownLoad = "r"==readline("To attempt to download a current cross reference table of HUGO symbols enter r \n(note: this can take more than 10 minutes to download)\nIf you tried this once and it didn't work, press enter to get other options. ")
	}
	
	if(autoDownLoad){
		full_hurl = "http://www.genenames.org/cgi-bin/hgnc_downloads?title=HGNC+output+data&hgnc_dbtag=on&preset=all&status=Approved&status=Entry+Withdrawn&status_opt=2&level=pri&=on&where=&order_by=gd_app_sym_sort&limit=&format=text&submit=submit&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag"
		cat("\nConnecting to HGNC website...\n")
		reopenedfurl=try(expr=read.table(file=full_hurl,sep="\t",comment.char="",
																		 header=T,quote="", stringsAsFactors=F,na.strings="-"),
										 silent=T)
		if(sum(grep(pattern="error", x=class(reopenedfurl), ignore.case=T))){
			cat("\nError text:\n")
			print(reopenedfurl)
			cat("\nError, could not download Hugo reference table from HGNC..\n")
		}else{
			cat("\nTable downloaded, writing to file...\n")
			dir.create(path=dirname(path=curhugofname), showWarnings=F, recursive=T)
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
	cat("Place the file in this location:", curhugofname,"\n")
	readline("Press enter to continue.\nPress escape to exit the program so that you can update the cross ref file.")
	cref = downloadHugoLookupTable(curhugofname)
	return(cref)
}


getAndStripHGNC<-function(){
	htabFile = "~/tprog/main_131219/reference_data/current_hugo_table.txt"
	htab = read.delim(file=htabFile, header=T, sep="\t", stringsAsFactors=F,quote="", na.strings="-")
	
	#first remove all the unapproved symbols
	withRows = grepl(pattern="~", x=htab$Approved.Symbol)
	
	apptab = htab[!withRows,]
	
	extab = apptab[,c("Approved.Symbol","Previous.Symbols")]
	
	write.table(x=extab, file="./htabExtract.txt", sep="\t",row.names=F, col.names=T)
	
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
