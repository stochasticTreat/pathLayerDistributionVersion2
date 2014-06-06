#PolyPhenProcessing_functions.R

getNumericScoreColumn<-function(xcol, pdat){
	#add column indicating numeric score for each polyphen score
	ppmn = c("benign","probably damaging","possibly damaging","unknown", "\\N")
	ppsc = c(1,3,2,0,-1)
	names(ppsc)<-ppmn
	nscol = rep(0, nrow(pdat))
	for(n in ppmn){
		nscol[grep(pattern=n, x=xcol$prediction)] = ppsc[n]
	}
	return(nscol)
}

#'@title orchestrates creation of polyphen input file from .maf file
#'@description Allows user to create valid polyphen input file from TCGAs somatic mutation data storage format, .maf.
#'@param mafFname character string. The file name and full or relative path to the maf input file
#'@param outFname character string. The name of the file for the polyphen input file to be saved to.
#'@return A list of diagnostic information, with each slot containing a specific data point, such as "Rows of mutation data found in .maf file input"
#'@export
PolyPhenFromMaf<-function(mafFname=NULL, outFname=NULL){
	#orchestrates creation of polyphen input file from .maf file
	#takes: mafData: name of .maf file
	tracker = list()
	mafData = read.delim(file=mafFname, 
						 header=T, 
						 sep="\t", 
						 stringsAsFactors=F, 
						 na.strings="-")
	
	if(!"pid"%in%colnames(mafData)){
		mafData = addPidColumn(tcga_data=mafData)
	}
	
	tracker[["Rows of mutation data found in .maf file input"]] = nrow(mafData)
	tracker[["Number of patient IDs found in .maf file"]] = length(unique(mafData$pid))
	datf1 = filtToPolyPhenTypes(mafData=mafData)
	tracker[["Rows of mutation data found after filtering to only PolyPhen types"]] = nrow(datf1)
	tracker[["Number of patient IDs found in mutation data after filtering by mutation types"]] = length(unique(datf1$pid))
	datf1 = addMappingColumn(datf1)
	ppcols = makePolyPhenCols(datf1)
	if(is.null(outFname)){
		outFname = strsplit(mafFname, split="/")[[1]][length(strsplit(mafFname, split="/")[[1]])]
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
	indexCol1 = 1:nrow(datf1)
	indexCol = paste(indexCol1, datf1$Hugo_Symbol, datf1$Start_Position, datf1$pid, sep="|")
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
	#first add an extra column headers to the first row: "pp.ref.dat" and "map.col"
	if(!length(grep(pattern="ref.col",prepdat[1,]))){
		prepdat[1,] = paste(prepdat[1,], "pp.ref.dat","map.col", sep="\t")
	}
	postind = grepl(pattern="^## ",x=prepdat[,1])
	if(sum(postind)) prepdat=prepdat[!postind,,drop=F]
	
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