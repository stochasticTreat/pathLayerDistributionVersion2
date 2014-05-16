#PolyPhenProcessing_functions.R

# runOHSUPolyPhenData<-function(fname, path_detail){
# 	#runOHSUPolyPhenData is the main interface for accepting polyphen data from OHSU
# 	#open file
# 	pdat = 
# 		#process file, open up scores
# 		pdatOpen = splitScoresOut(pdat)
# 	#process file, make correct columns
# 	
# 	#process polyphen data
# 	psum = processPolyPhen(polyDat=pdatOpen,paths_detail=path_detail)
# 	return(psum)
# }

getNumericScoreColumn<-function(xcol, pdat){
	#add column indicating numeric score for each polyphen score
	ppmn = c("benign","probably damaging","possibly damaging","unknown", "\\N")
	ppsc = c(1,2,3,0,-1)
	names(ppsc)<-ppmn
	nscol = rep(0, nrow(pdat))
	for(n in ppmn){
		nscol[grep(pattern=n, x=xcol$prediction)] = ppsc[n]
	}
	return(nscol)
}

# PPMultiMapRedux<-function(pdat){
# 	#takes the formatted polyphen output and reduces the rows to only include
# 	# the max score for each index
# 	#return: data frame, same structure as input; same data, but rows removed; no multi mapping to indexes
# 	#extract the needed columns
# 	cat("\nHandling situations where there is more than one variants in the same gene,\nwith those variants recieving more than one unique PolyPhen score.\nIn these cases, the gene will be assigned the score indicating maximum variant affect of all variant affects attibuted to the gene.\n")
# 	
# 	xcol = pdat[,c("index","pid","Symbol","prediction")]
# 	
# 	numscore = getNumericScoreColumn(xcol = xcol, pdat = pdat)
# 	cat("\nPolyphen classifications traslated...")
# 	xcol = cbind(xcol, numscore)
# 	#find the unique set of indexes
# 	ui = unique(xcol$index)
# 	gindexes = c()
# 	cat("\nFinding maximum scores..\n")
# 	for(i in ui){# for each index, find the max polyphen score
# 		seti = which(xcol$index == i)
# 		sel = which(xcol$numscore[seti] == max(xcol$numscore[seti]))
# 		gindexes = c(gindexes, seti[sel[1]])
# 	}
# 	redux = xcol[gindexes,]	
# 	cat("\n.\n")
# 	return(redux)
# }
 
PolyPhenFromMaf<-function(mafFname=NULL, outFname=NULL){
	#orchestrates creation of polyphen input file from .maf file
	#takes: mafData: name of .maf file
	tracker = list()
	
	mafData = read.table(comment.char="",
											 file=mafFname,
											 header=T,
											 sep="\t", 
											 stringsAsFactors=F)
	
	tracker[["Rows of mutation data found in .maf file input"]] = nrow(mafData)
	tracker[["Number of patient IDs found in .maf file"]] = length(unique(mafData$pid))
	datf1 = filtToPolyPhenTypes(mafData=mafData)
	tracker[["Rows of mutation data found after filtering to only PolyPhen types"]] = nrow(datf1)
	tracker[["Number of patient IDs found in mutation data after filtering by mutation types"]] = length(unique(datf1$pid))
	datf1 = addMappingColumn(datf1)
	ppcols = makePolyPhenCols(datf1)
	if(is.null(outFname)){
		outFname = strsplit(fname, split="/")[[1]][length(strsplit(fname, split="/")[[1]])]
		outFname = paste("./output/", outFname, ".polyphenInput.txt", sep="")
	}
	write.table(x=ppcols, file=outFname, sep=" ", quote=F, row.names=F,col.names=F)
	tracker[["PolyPhen input file name"]] = outFname
	cat("\nPolyPhen input file written to",outFname, "\n")
	return(tracker)
}

makePolyPhenCols<-function(datf1){
	#makes the actual columns of data accepted by polyphen
	#takes: datf1: data frame with columns: Chrom, Start_Position, Reference_Allele, Tumor_Sel_Allele2, indexCol
	#returns: data frame ready to be saved as input file for polyphen
	locCol = paste("chr",datf1$Chrom,":",datf1$Start_Position, sep="")
	allele = paste(datf1$Reference_Allele,"/",datf1$Tumor_Seq_Allele2, sep="")
	idcol = paste("#", datf1$indexCol)
	finalPPinputDF = cbind(locCol,allele, idcol)
	return(finalPPinputDF)	
}

addMappingColumn<-function(datf1){
	indexCol = 1:nrow(datf1)
	indexCol = paste(indexCol,datf1$Hugo_Symbol, datf1$Start_Position, datf1$pid, sep="|")
	datf1 = cbind(datf1, indexCol)
	return(datf1)
}

filtToPolyPhenTypes<-function(mafData, mutTypes = c("Missense_Mutation")){
	
	#filt to only the selected mut types
	datf = mafData[mafData$Variant_Classification%in%mutTypes,]
	#remove mitochondrial genes
	datf0 = datf[datf$Chrom!="MT",]
	#remove ins/del
	datf1 = datf0[datf0$Variant_Type == "SNP",]
	#for the time being, just remove the NT_ rows, such as those with chromosome annotation NT_113923
	ntrows = grep(pattern="NT_", x=datf1$Chrom)
	datf1 = datf1[!1:nrow(datf1)%in%ntrows,]
	return(datf1)
}

selectionList<-function(valcol, verbose=T){
	
	usel = unique(valcol)
	uselmat = matrix(data=1:length(usel), ncol=1, dimnames=list(usel, "selection number"))
	uin = readline(paste("Please enter the number(s) corresponding to the PolyPhen values you would like to consider aberrational.\n",
											 "Press enter for the default set: (probably_damaging)"))
	if(uin=="") uin = "1"
	uin = as.integer(strsplit(x=uin, split=" ")[[1]])
	lvout = valcol%in%usel[uin]
	return(lvout)
}

splitScoresOut<-function(seqdat){
	#takes polyphen data (data frame with PolyPhen column) as provided in OHSU data
	#breaks out the score and polyphen category
	#returns table with score and category appended
	polyPhen_full = seqdat$PolyPhen
	if(length(grep(pattern="\\(",x=polyPhen_full))){
		scores  = gsub(x=polyPhen_full,replacement="",pattern=")")
		ppcategory = rep("", times=length(scores))
		ppscore = rep(0, times=length(scores))
		for(i in 1:length(scores)){
			tmpscores = strsplit(x=scores[i], split="\\(")[[1]]
			if(length(tmpscores)==2){
				ppcategory[i] = tmpscores[1]
				ppscore[i] = as.numeric(tmpscores[2])
			}
		}
		seqdat$PolyPhen=ppcategory
		out = cbind.data.frame(seqdat, 
													 polyPhen_full,
													 ppscore, 
													 stringsAsFactors=F)
		return(out)
	}
	return(seqdat)
}

parseConsequence<-function(cons){
	out = c()
	for(i in 1:length(cons)){
		cur = strsplit(x=cons[i], ";")[[1]]
		out = union(out, cur)
	}
	return(out)
}

preProcessPdat<-function(fname){
	cat('\nPreprocessing PolyPhen output from file', fname, "....")
	#removes/corrects unwanted formatting in initial PolyPhen output and prepares file for regular input
	prepdat = read.table(file=fname, header=F,sep="\n",comment.char="", stringsAsFactors=F)
	#first add an extra column header to the first row
	if(!length(grep(pattern="ref.col",prepdat[1,]))){
		prepdat[1,] = paste(prepdat[1,], "pp.ref.dat","map.col", sep="\t")
	}
	postind = grepl(pattern="^## ",x=prepdat[,1])
	prepdat=prepdat[!postind,,drop=F]
	newFname = paste(fname, ".reformatted.txt",sep="")
	write.table(x=prepdat, file=newFname, sep='\n',row.names=F, col.names=F, quote=F)
	cat('Preprocessing complete\n')
	return(newFname)
}

# loadPolyPhenResults<-function(fname){
# 	
# 	newFname = preProcessPdat(fname)
# 	pdat = read.table(file=newFname, sep="\t",header=T,comment.char="", stringsAsFactors=F)
# 	#first parse out the map.col
# 	mc = pdat$map.col
# 	mcs = data.frame(matrix(byrow=T,data=unlist(strsplit(mc, split="\\|")), ncol=4), stringsAsFactors=F)
# 	colnames(mcs)<-c("index", "Symbol", "Start.Pos", "pid")
# 	pdat = cbind(pdat, mcs)
# 	return(pdat)
# }


runMAFPolyPhen<-function(fname){
	
	#get polyphen output file
	#re-formatt it
	#break out the mapping
	pdat = loadPolyPhenResults(fname=fname)
	
	fpdat = PPMultiMapRedux(pdat = pdat)
	cat("\n.\n")
	#find max for each gene/patient: for each index, find the max polyphen score
	#possibly do this again for each patient/gene combo?
	#change column names to make them appropriate for processPolyPhen
	colnames(fpdat)<-c("index", "alias", "Symbol", "PolyPhen", "numscore")
	psum = processPolyPhen(polyDat=fpdat,paths_detail=path_detail)
	return(psum)
}

processPolyPhen<-function(polyDat, paths_detail, disease_type="disease type not given", threshold=NULL){
	#processPolyPhen()
	#orchestrates processing of PolyPhen data
	#takes: polyDat: the table, in stacked format, of polyphen data
	#									must have columns "alias", "Symbol", "PolyPhen"
	#				disease_type: the study name/ ideally, the type of disease being examined
	#				paths_detail: the path_detail object
	#				thresh: 			The polyphen values accepted to indicate genes are "active"
	cat("\nProcessing polyphen data.. \n")
	
	tracker = list()#establish the tracker to document data work-up
	
	tracker[["Number of mutation records in input data set"]] = nrow(polyDat)
	
	#### remove rows with "\N" as their symbol
	badrows = grep("\\N", x=polyDat$Symbol)
	tracker[["Number of rows/records with no gene symbol given (\"\\N\" in the place of gene symbol)"]] = length(badrows)
	cat("\nThere were",
			length(badrows),
			"records which had no HUGO symbol associated,\nand had \"\\N\" in the place of a HUGO symbol.")
	polyDat = polyDat[!(1:nrow(polyDat)%in%badrows),]
	
	##### check / correct gene names
	polyDat$Symbol = corsym(symbol_set=polyDat$Symbol, 
													verbose=T,
													symref=path_detail$symtable)
	
	#### clean polyphen values, removing any that are marked "unknown" or ""
	polyDat = splitScoresOut(seqdat=polyDat)
	
	#reduce coverage, removing all the genes that are "" and "unknown", as to polyphen output
	selcover  = polyDat[!polyDat$PolyPhen%in%c("", "unknown"),]
	
	#select threshold/on/off set
	logicvector = selectionList(valcol=selcover$PolyPhen, verbose=verbose)
	
	#now make a patient gene matrix
	forPGM = selcover[logicvector,]
	colnames(forPGM)[1] = "Ids"
	ppgm = toPGM(sds=forPGM)
	
	targmat = getTargetMatrix(tgenes=selcover$Symbol, paths=paths_detail$paths)
	
	polyPhenSummary = summaryTable4(paths_detail=paths_detail,
																	individualEnrichment=T,
																	verbose=T,
																	target_matrix=targmat,
																	dataSetName=disease_type,
																	patientGeneMatrix=ppgm, 
																	targetname="mutated")
	return(polyPhenSummary)
}

filePrompt<-function(defaultfile = "./input/OHSUseqPolyPhencoding_and_nodbSNP_rows.txt.reformatted.txt"){
	#prompts user to select file for data input
	#provides default file option
	#returns file name
	fsel = readline(paste("\nTo select a file of PolyPhen data, Enter s\n",
												"To load the default AML PolyPhen data from \n",defaultfile,",\njust press enter \n",sep=""))
	if(fsel=="s"){
		pfile = file.choose()
	}else{
		pfile = defaultfile
	}
	cat("\nLoading PolyPhen data from:\n",pfile,"\n")
	return(pfile)	
}