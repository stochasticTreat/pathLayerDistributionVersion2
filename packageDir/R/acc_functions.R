
if(!exists("VERBOSE")) VERBOSE = F

.pardefault <- par(no.readonly = T)
.parpin <-par()$pin
# .pardefault <- par()


stackedGeneBar<-function(tcga_som, 
												 title="Mutation types for the top 20 most mutated genes", 
												 colorMatchesFile="./reference_data/MAFcolorMatches.txt"){
	ufilt = tcga_som
	gs = summarize_by(col=tcga_som[,"Hugo_Symbol"], display=F)
	top20 = head(gs[order(gs[,2],decreasing=T),],20)
	gsnames = top20$types
	
	toprows = ufilt[ufilt$Hugo_Symbol%in%gsnames,]
	
	slimrows = toprows[,c("Hugo_Symbol","Variant_Classification")]
	slimrows = merge(x=slimrows, y=top20, by.x="Hugo_Symbol", by.y="types")
	slimrows = slimrows[order(slimrows$counts),]
	slimrows$Hugo_Symbol<-factor(slimrows$Hugo_Symbol, levels=top20$types)
	
	colcols = getColorSequence(cnames=unique(slimrows$Variant_Classification), 
														 fname=colorMatchesFile)
	colnames = colors()[colcols]
	names(colnames)<-names(colcols)
	
	p1 = ggplot(data=slimrows, aes(x=Hugo_Symbol, 
																 fill=factor(Variant_Classification)))+
		geom_bar(color="black")+
		coord_flip()+
		scale_fill_manual(values=colnames)+
		theme_bw()+
		theme(legend.title=element_blank())+
		ggtitle(title)
	
	print(p1)
}

simpleGGHist<-function(dataSet, xlab, ylab, mainTitle, showPlot=T){
	require(ggplot2)
	p = qplot(x=as.vector(dataSet), geom="histogram")+
		theme_bw()+
		xlab(xlab)+
		ylab(ylab)+
		ggtitle(mainTitle)
	
	if(showPlot){
		print(p)
	} else {
		return(p)
	}
}

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

#addPidColumn
#extracts TCGA barcodes, and reformats them so dashes (-) are replaced with periods (.) . 
#param tcga_data data.frame with a column containing the string Sample_Barcode (if more than one are found, the first one will be used), from which the base tcga barcode will be extracted.
#return the tcga_data data.frame with a column appended containing the extracted pid, in the format TCGA.AB.2988, for each row. 
addPidColumn<-function(tcga_data){
	cat("\nFixing patient ids...\n")
	colIndex = grep(pattern="Sample_Barcode", x=colnames(tcga_data), ignore.case=T)[1]
	pid = sapply(as.character(tcga_data[,colIndex]), extract_pid)
	tcga_data = cbind.data.frame(pid, tcga_data, stringsAsFactors=F)#append extracted pids as a sepparated column
	return(tcga_data)
}


longTextBarPlot<-function(data, lab, main=""){
	
	bpdata  = barplot(data, horiz=T, main=main)
	
	xmin = par("usr")[1]
	xrange = par("usr")[2] - par("usr")[1]
	xcoord = xmin + (xrange/40)
	
	nlev = length(data)
	
	text(x=xcoord, y=bpdata, labels=lab, pos=4)
	
}



###############################################################################
#
# helperFunctions.R: 	This file contains the all helper functions not directly related to any other source file.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Replace factors/levels in a data.frame and use plain strings instead 
#' 
#' This function takes a data.frame as argument and returns it with strings instead of factors.
#' 
#' @param df any data.frame with factor levels in at least one column
#' @return The data.frame is returned using strings instead of factors.
#' @author Frank Kramer
unfactorize <- function (df) {
	if(!("data.frame" %in% class(df))) { stop("Error: unfactorize: data.frame as argument expected") }
	
	for (i in 1:ncol(df)) {
		if (class(df[,i]) == "factor") df[,i] <- as.character(df[,i])
	}
	df
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


plotFileName<-function(pname, addTime=F){
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
	ctime = ""
	if(addTime) ctime = gsub(pattern=" ", replacement=".", x=paste0(as.character(Sys.time()),"."))
	
	fullPath = paste(root, "graphic", gsub(pattern="[-:]", replacement=".", x=ctime), pname, sep="")
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
selectionList<-function(valcol, verbose=T){
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

#'@title openPGM
#'@description Opens a patient gene matrix file.
#'@param fname The file path to the patient gene matrix file.
#'@return The patient gene matrix (logic matrix with gene identifiers given as row names and patient identifiers gien as column names)
#'@export
openPGM<-function(fname = NULL){
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
#'@import ggplot2
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
		p = ggplot(out, aes(x=types, y=counts))+
			geom_bar(stat='identity')+
			coord_flip()+
			ggtitle(barPlotTitle)
		print(p)
		# 		oldmar <- par()$mar
		# 		while(T){
		# 			res = try({
		# 				par(mar=c(5.1, max(4.1,max(left_margin_factor*nchar(types))/2.5) ,4.1 ,2.1))
		# 				#	try(displayGraph(w), silent=T)
		# 				barplot(counts, horiz=T, las=2, main = barPlotTitle, xlab="Number found in data set", names.arg=types)
		# 				par(oldmar)
		# 			}, silent=T)
		# 			# 			if(!grepl(pattern="Error", x=res)) break
		# 			if(is.null(res)) break
		# 			par(oldmar)
		# 			readline(prompt="There seems to have been an error with plotting the bar graph.\nPlease increase the size of the plot window, the press enter")
		# 		}
		# 		par(oldmar)
		# 		print(out)
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
# 
startHTMLPlugIns<-function(){

	if(!require("hwriterPlus")){
		print("Trying to install hwriterPlus so that HTML output can be generated")
		install.packages("hwriterPlus")
		if(require("hwriterPlus")){
			print("hwriterPlus installed and loaded")
		} else {
			stop("could not install hwriterPlus")
		}
	}
	if(!require("hwriter")){
		print("Trying to install hwriterPlus so that HTML output can be generated")
		install.packages("hwriter")
		if(require("hwriter")){
			print("hwriter installed and loaded")
		} else {
			stop("could not install hwriter")
		}
	}
	if(!require("xtable")){
		print("Trying to install hwriterPlus so that HTML output can be generated")
		install.packages("hwriterPlus")
		if(require("xtable")){
			print("xtable installed and loaded")
		} else {
			stop("could not install xtable")
		}
	}

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
	
	#startHTMLPlugIns()
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
								 reorderTables=NULL, 
								 pagetitle=NULL, 
								 limit_col=NULL, 
								 path_detail=NULL){
	cat("____________________________Saving HTML page to :",fname,"\n")
	#startHTMLPlugIns()
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
					SaveToHTML_inner(results=eaPatSet[[patsetname]], 
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
							 reorderTables=reorderTables,
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
								if(!is.null(reorderTables)){
									print("reorderTables")
									curorder = order(tmptab[,c],decreasing=reorderTables)
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
						if(!is.null(reorderTables)){#order the rows if wanted
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


#'@title Save results to an HTML page. 
#'@description Writes data from the  \code{results} slot of a \code{Study} object to an HTML document. 
#'@param study A \code{Study} object with results sets from one or more study arm. 
#'@import xtable
#'@import hwriterPlus
#'@import hwriter
#'@export
SaveToHTML<-function(study){
# 	if(!require("hwriterPlus")){
# 		print("Trying to install hwriterPlus so that HTML output can be generated")
# 		install.packages("hwriterPlus")
# 		if(require("hwriterPlus")){
# 			print("hwriterPlus installed and loaded")
# 		} else {
# 			stop("could not install hwriterPlus")
# 		}
# 	}
# 	if(!require("xtable")){
# 		print("Trying to install hwriterPlus so that HTML output can be generated")
# 		install.packages("hwriterPlus")
# 		if(require("xtable")){
# 			print("xtable installed and loaded")
# 		} else {
# 			stop("could not install xtable")
# 		}
# 	}
	SaveToHTML_inner(study_name=study@studyMetaData@studyName,
						 results=study@results, 
						 path_detail=study@studyMetaData@paths)	
}

SaveToHTML_inner<-function(results,
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
			htmlPerPatient(results=results, summ=n, fileroot=root2, overlap=F, path_detail=path_detail)
		}else if(grepl(pattern="overlap_analysis$", x=n, ignore.case=T)){
			print("overlap analysis selected")
			
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


test.twoHistOnePlot<-function(){
	
	prefilt = data.frame(muts = sample(x=100:200, replace=T, size=100), stringsAsFactors=F)
	postfilt = data.frame(muts=sample(x=50:150, replace=T, size=100), stringsAsFactors=F)
	twoHistOnePlot(dataset1=prefilt, 
								 dataset2=postfilt, 
								 x_label="Number of mutations", 
								 y_label="Number of patients",
								 legend_titles=c("Before filtering", "After filtering"),
								 main_title="Distributions of mutations in patients before and after filtering")
}

#twoHistOnePlot
#
#'@title Puts two histograms on the same plot.
#'@description Puts two histograms on the same plot, shows each with a different color, and overlap with a third color.
#'@param dataset1 A numeric data set; \code{data.frame or matrix} with one column.
#'@param dataset2 A numeric data set; \code{data.frame or matrix} with one column.
#'@param frequency_not_density Depricated. 
#'@param main_title The main title for the plot.
#'@param x_label The label to be put along the x axis. 
#'@param y_label The label to be put along the y axis. 
#'@param legend_titles Character vector, length 2, giving the names of the data in dataset1 and dataset2. 
#'@param breaks Depricated. 
#'@import ggplot2
twoHistOnePlot<-function(dataset1, dataset2, 
												 frequency_not_density=T, 
												 main_title="main_title not set",
												 x_label="x_label not set", 
												 y_label="Frequency",
												 legend_titles=c("dataset 1","dataset 2"),
												 breaks=c(30,30)){
 	cat("\nReporting distributions of ", legend_titles[1],"and",legend_titles[2],"...")
	#first combine the two data sets
	colnames(dataset1)<-"dat"
	colnames(dataset2)<-"dat"
	
	stackedSet = rbind.data.frame(dataset1, dataset2)
	idCol = c(rep(legend_titles[1],times=nrow(dataset1)), rep(legend_titles[2], times=nrow(dataset2)))
	stackedSet = cbind.data.frame(stackedSet, idCol)
	colnames(stackedSet)<-c("value","idCol")
	
	# 	pmin = min(c(dataset1, dataset2))
	# 	pmax = max(c(dataset1, dataset2))
	# 	
	# 	p1 <- hist(dataset1, freq=frequency_not_density,breaks=breaks[1],)                     
	# 	p2 <- hist(dataset2, freq=frequency_not_density,breaks=breaks[2])
	# 	
	# 	hist = plot( p2, col=rgb(0,0,1,1/4), xlim=c(pmin,pmax), 
	# 							 freq=frequency_not_density, 
	# 							 main=main_title, 
	# 							 xlab=x_label, 
	# 							 ylab = y_label)  # first histogram
	# 	hist2 = plot( p1, col=rgb(1,0,0,1/4), xlim=c(pmin,pmax), 
	# 								add=T,
	# 								freq=frequency_not_density)  # second
	# 	
	# 	if(!is.null(legend_titles)){
	# 		legend("topright", inset=.05,
	# 					 legend_titles, 
	# 					 fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), 
	# 					 horiz=F)
	# 	}
	bwidth = ifelse(test=length(unique(stackedSet$value))>10, yes=3, no=1)
	p1 = ggplot(stackedSet, aes(x=value, fill=idCol)) + 
							geom_histogram(alpha=0.5, position="identity", binwidth=bwidth)+
							ggtitle(main_title)+
							theme_bw()+
							theme(legend.title=element_blank())+
							xlab(x_label)
# 	if(min(stackedSet[,1])<=0){
# 		rng = (range(stackedSet$value)[2]-range(stackedSet$value)[1])*.1
# 		p1 = p1+scale_x_continuous(limits=c( (min(stackedSet[,"value"])-rng), (max(stackedSet[,"value"])+rng ) ))
# 	}
	print(p1)
	cat("...done\n")
	return(p1)
}

scatterhist <- function(x, y, xlab="", ylab="", main=""){
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

scatterhist2 <- function(x1, y1, x2, y2, xlab="", ylab="", main="", legendTxt=NULL){
	
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
	symtab = path_detail$symtable
	
	# 	#assure there are no symbols given as "\\N" 
	# 	genelist = genelist[genelist!="\\N"]
	
	#attempt to make sure the symbols match 
	genelist = corsym(genelist, symref=symtab, verbose=T)
	
	#extract the needed columns
	xsymtab = symtab[,c("Approved.Symbol", "Approved.Name", "Locus.Type", "Chromosome", "Previous.Names")]
	rownames(xsymtab) = xsymtab$Approved.Symbol
	
	#figure out which gene symbols we actually have info on:
	have = genelist%in%xsymtab$Approved.Symbol
	
	#reduce the gene list to those we have info on and report those we dont have info on
	goodset = genelist[have]
	badset = genelist[!have]
	
	out=list()
	out$'Genes and their records' = xsymtab[goodset,]
	out$'Gene symbols for which records could not be found' = badset
	
	return(out)
}

