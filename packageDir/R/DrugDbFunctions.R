
getDrugEav<-function(fname){
	curDrug = ""
# 	readline("Check for situations where the attribute name is blank... press enter to continue.")
	#first scan the file to determine how many rows
	cat("\nReading file", fname,"\n")
	rawtab = read.table(file=fname, 
											header=F, 
											sep="\n", 
											quote="", 
											comment.char="", 
											check.names=F, 
											stringsAsFactors=F, 
											na.strings="")
	#first count how many unique drugs
	drugcount = 0
	hashcount = 0
	print("Surveying input file... ")
	for(i in 1:nrow(rawtab)){
		if(grepl(pattern="^#", x=rawtab[i,1])){
			hashcount=hashcount+1
			if(grepl(pattern="^#BEGIN_DRUGCARD", x=rawtab[i,1])){
				cat(".",i,".")
				drugcount=drugcount+1
			}
		}
	}
	print("Scan complete")
	#make the receptical data frame
	outmat = matrix(data="", nrow=(hashcount-drugcount), ncol=3, dimnames=list(NULL, c("entity", "attribute", "value")))
	
	#fill the outmat
	outi=1
	inputi=1
	endline=""
	curdrug=""
	duratt = ""
	while(T){
		
		next_line=rawtab[inputi,1]
		print(inputi)
		if(grepl(pattern="^#", x=next_line)){#control: set drug or attribute
			if(grepl(pattern="^#BEGIN_DRUGCARD",x=next_line, ignore.case=T)){
				#parse the drug name out and put the records up
				curdrug = gsub(pattern="#BEGIN_DRUGCARD ", replacement="", x=next_line)
				print(curdrug)
				print(inputi)
				print(outi)
				endline = paste("#END_DRUGCARD", curdrug)
			}else if(next_line!=endline){
				curatt = gsub("# ", replacement="", x=next_line)
			}
			
			inputi=inputi+1
			if(inputi == nrow(rawtab)) {
				print("nrow(rawtab")
				break
			}
		}else{#set values
			curval=next_line
			sepChar = ifelse(test=grepl(pattern="sequence", x=curatt, ignore.case=T), 
											 no="|",
											 yes="")
			while(T){
				#incrament the inputi counter
				inputi=inputi + 1
				next_line=rawtab[inputi,1]
				#check if there's a hash tag
				if(grepl(pattern="^#", x=next_line)) break
				if(inputi>=nrow(rawtab)) break
				#if there's not a has tag, paste it onto curval
				curval = paste(curval, next_line, sep=sepChar)
			}
			print("Attribute aquired.")
			if(outi>nrow(outmat)) break
			outmat[outi,]=c(curdrug, curatt, curval)
			outi=outi+1
			print("next")
		}
	}
	print("Finished transformation.")
	return(outmat)
}

removeTabs<-function(vec){
	
	vec = gsub(pattern="\t", replacement=" ", x=vec)
	
}


getDrugTargetRecords<-function(dnames="*", eavtab){
	
	ptm<-proc.time()
	targrows = eavtab[grepl(pattern="Drug_Target_[0-9]+_", x=eavtab$attribute),]
	totalTime = proc.time()  - ptm
	cat("Operation took", totalTime["elapsed"], "seconds.\n")
	
	ptm<-proc.time()
	drugTargetRows = targrows[grepl(pattern=dnames, x=targrows$entity, ignore.case=T),]
	totalTime = proc.time()  - ptm
	cat("Operation took", totalTime["elapsed"], "seconds.\n")
	
	return(drugTargetRows)
}


loadDrugEAVtable<-function(fileName="./reference_data/drugDB/drugbank/drugbankEAV.txt"){
	cat("Loading the drug data file", fileName,"... \n")
	ptm<-proc.time()
	tabIn = read.table(file=fileName, header=T, sep="\t", quote="", check.names=F, comment.char="", stringsAsFactors=F)
	totalTime = proc.time()  - ptm
	cat("File load took", totalTime["elapsed"], "seconds.\n")
	return(tabIn)
}

# 
# test.getAtts<-function(){
# 	cat("\nThis should take about 20 seconds \n")
# 	
# 	eavtab=loadDrugEAVtable()
# 
# 	# 	> dim(eavtab)
# 	# 	[1] 933299      3
# 	
# 	ptm2<-proc.time()
# 	targrec = getDrugTargetRecords(eavtab=eavtab)
# 	totalTime = proc.time() - ptm2
# 	cat("Getting drug target records took", totalTime["elapsed"], "seconds.\n")
# 	# 	> dim(targrec)
# 	# 	[1] 452414      3
# 	
# 	retAtts = getAtts(mat=targrec, atts=c("^ID:", 
# 																				"^Name", 
# 																				"Pfam_Domain_Function", 
# 																				"HGNC_ID",
# 																				"Gene_Name", 
# 																				"General_Function"))
# 	# 	> dim(retAtts)
# 	# 	[1] 87564     3
# 	
# 	#for each 
# 	#get table of targets
# 	#needs: drug, target gene symbol, reference
# 	
# 	#first: pull out the individual target symbols for all the individual drugs
# 	targSyms = getAtts(mat=targrec, atts="Gene_Name:$")
# 	tsyms=targSyms
# }


# getDrugTargetTable<-function(tsyms){
# 	eavtab=loadDrugEAVtable()
# 	
# 	# 	> dim(eavtab)
# 	# 	[1] 933299      3
# 	
# 	targrec = getDrugTargetRecords(eavtab=eavtab)
# 
# 	# 	> dim(targrec)
# 	# 	[1] 452414      3
# # 	targrec = getDrugTargetRecords(eavtab=eavtab)
# 	
# 	#first get the unique drugs
# 	udrugs = unique(targrec$entity)
# 	
# # 	1999 Jul;22(3):231-8. 1349838 Iwahana
# # 	regpat = "[0-9]{4} [a-zA-Z]+( [0-9]+)*[;]"#"#[(][0-9]+[)]:[0-9]+-[0-9]+[.] [0-9]+ [A-Za-z]+"
# 	regpat = ";"
# 	
# 	for(d in udrugs){
# 		drecs = targrec[targrec$entity==d,]#pull out the records for the drug
# 		#for each drug get all the targets
# 		targSyms = getAtts(mat=drecs, atts="Gene_Name:$")
# 		targets = targSyms$value
# 		#for each target, get any references
# 		targRefs = getAtts(mat=drecs, atts="General_References:$")
# 		refs = strsplit(x=targRefs$value, split=regpat, perl=T)
# 	}
# 	
# 	
# }

# #'@title getAtts
# #'@description Gets all rows in the matrix, mat, with a the attribute strings(s) in atts
# #'@param mat the matrix to be searched through; contains at least 2 columns with column to being the attributes column. 
# #'@param atts a vector of strings; the attributes to be matched. perl compatible regular expressions are implemented for string matching
# #'@param verbose Flag indicating if extra information should be displayed.
# #'@return matrix, the subset of rows matching the strings provided in atts. 
# getAtts<-function(mat, atts, verbose=F){
# 	omat = mat
# 	mat[,2] = gsub(pattern="Drug_Target_[0-9]*_", replacement="", x=mat[,2])
# 	selli = rep(F, times=nrow(mat))
# 	for(a in atts){
# 		if(verbose) cat(a, "retreived:\n")
# 		curi = grepl(pattern=a, x=mat[,2], ignore.case=T, perl=T)
# 		selli = selli | curi
# 		if(verbose) print(head(omat[curi,]))
# 	}
# 	retmat = omat[selli,]
# }