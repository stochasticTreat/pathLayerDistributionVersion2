#processOHSUseqDat_functions.R

#'@title Add sequence capture data input arm.
#'@description Add arm onto main study for input of sequence capture data. 
#'@param STUDY A \code{Study} object.
#'@return STUDY A \code{Study} object with the sequence capture data input arm loaded. 
#'@export
#'@examples
#'STUDY = getTestStudyObject()
#'STUDY = addSequenceCaptureArm(STUDY=STUDY)
#'\dontrun{
#'STUDY = allInteractiveMainFunction(additionalArms=addSequenceCaptureArm)
#'}
addSequenceCaptureArm<-function(STUDY){
	arms = STUDY@arms
	arms = loadDataArm(description="Load OHSU sequence capture data",
										 title="sequence_capture_aberration_summary", 
										 mainFunction=processSequenceCaptureData, 
										 arms=arms)
	STUDY@arms = arms
	return(STUDY)
}


loadCoverageData<-function(path_detail, fname, verbose=F){
	#function to load seq.capture  gene coverage data
	tab = read.table(file=fname,header=F,stringsAsFactors=F)
	#check that usym are approved hugo symbols
	usym = corsym(symbol_set=tab[,1], symref=path_detail$symtable, verbose=verbose)
	usym = unique(usym)
	return(list(file=fname, cov=usym))
}#loadCoverageData

appendPPScoresOHSU<-function(oseqDat){
	li = oseqDat$PolyPhen!=""
	npol = paste("PolyPhen_", oseqDat$PolyPhen[li], sep="")
	nvarcol = paste(oseqDat$Consequence[li], npol, sep=";")
	oseqDat$Consequence[li] = nvarcol
	return(oseqDat)
}

CheckIntegratePolyPhen<-function(seqDat, tracker){
	print(exists("seqDat"))
	#first check to see if there is a polyphen column
	if(length(grep(pattern="PolyPhen", colnames(seqDat), ignore.case=T))){
		cat("\nPolyPhen column found.. \n")
		tracker[["PolyPhen-annotated variants"]] = sum(seqDat$PolyPhen!="\\N")
		#if there's a polyphen column, adjust the consequences
		sd = splitScoresOut(seqdat=seqDat)
		print("Split PolyPhen scores out")
		sd2 = appendPPScoresOHSU(oseqDat=sd)
		print(dim(sd2))
		print("scores appended to consequence types")
		seqDat=sd2
	}
	return(list(seqDat=seqDat, tracker=tracker))
}

FilterdbSNP<-function(ohsuSeqDat, tracker, s){
	

	mess  = 	paste("Out of",
								 nrow(ohsuSeqDat),
								 "variants found in the cohort, the number with dbSNP records is",
								 sum( ohsuSeqDat$in_dbsnp==1))
	tracker[["Number of variants in chohort with dbSNP records:"]] = mess
	cat("\n",mess,"\n")
	
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	while(TRUE){
		s = setting(s=s,prompt="Would you like to filter out dbSNP Values? (y/n)")
		if(s$.text%in%c("y","n")) break
	}
	filtLogicVector = rep(TRUE, nrow(ohsuSeqDat))
	if(s$.text=="y"){
		filtLogicVector = ohsuSeqDat$in_dbsnp==0
	}
	return(list(tracker=tracker, filtLogicVector=filtLogicVector, s=s))
}


#for processing coverage data; not used anymore
splitMultiGeneRows<-function(syms){	
	crows = grepl(pattern=",", x=syms)
	cgenes = syms[crows]
	sgenes = c()
	
	for(d in cgenes){
		
		cur = strsplit(x=d, split=",")[[1]]
		sgenes = c(sgenes, cur)
	}
	
	lessrows = syms[!crows]
	
	out = union(lessrows, sgenes)
	return(out)
}

coverageFilePrompt<-function(defaultfile = "./OHSUseqDat/coords_with_names_4_16_2010.txt"){
	#prompts user to select file for data input
	#provides default file option
	#returns file name
	fsel = readline(paste("\nTo select a file of sequence capture coverage data, Enter s\n",
												"To load the default data from \n",defaultfile,",\njust press enter \n",sep=""))
	if(fsel=="s"){
		pfile = file.choose()
	}else{
		pfile = defaultfile
	}
	cat("\nLoading coverage data from:\n",pfile,"\n")
	return(pfile)	
}

stackedGeneBar_OHSUseq<-function(seqDat){
	#stackedGeneBar
	#takes: tcga_som: initial input .maf table, after cleaning of repeat rows and gene symbols
	#outputs: stacked bar plot
	ufilt = seqDat
	gs = summarize_by(col=seqDat[,"Symbol"], display=F)
	top20 = head(gs[order(gs[,2],decreasing=T),],20)
	gsnames = top20$types
	
	toprows = ufilt[ufilt$Symbol%in%gsnames,]
	
	slimrows = toprows[,c("Symbol","Consequence")]
	
	#set up the output matrix
	mainCons = unique(slimrows$Consequence)

	mainCons = unique(mainCons)
	mmat = NULL
	mmat = matrix(data=0,ncol=length(mainCons), 
								nrow=length(gsnames), 
								dimnames=list(gsnames, mainCons))

	for(n in gsnames){
		#get the consequences for the gene names
		cons = slimrows[grepl(pattern=n, x=slimrows$Symbol),"Consequence"]
		msum = summarize_by(col=cons, display=F)
		mmat[n,msum$types] = msum$counts
	}
	mdf  = cbind.data.frame(gene = rownames(mmat),mmat)

	par(las=2) # make label text perpendicular to axis
	par(mar=c(5,7,1,1)) # increase y-axis margin.
	oldmar =c(5,4,4,2)
	oldlas = 0
	
	scol = 2

	colcols = getColorSequence(cnames=colnames(mmat))

	barplot(t(mmat), 
# 					title.adj = 0,
# 					title.outer=T,
					legend = colnames(mmat),
					xlab="Number of variants found in gene", 
					cex.sub=1.2,
					cex.main=1.2,
					main="", 
					col=colors()[colcols],
					horiz=TRUE,
					cex.names=1.2, 
					cex.lab=1.2, 
					args.legend=list(x="topright",cex=.9, inset=c(.01,.001)))
	print("printing legend")
# 	legend("topright", 
# 				 legend=colnames(mmat))
	title(main="\nVariant types for the top 20 most mutated genes", 
				outer=F, 
				adj=0)
# 	legend.text=T,
	# 					legend.args=list(cex=1.2),
# 	legend = colnames(mmat),
	par(las=oldlas)
	par(mar = oldmar)
}

getColorSequence<-function(cnames, colorPalleteFile=NULL){
	#first check cnames, 
	#			if colors already assigned to cnames, 
	#					set those colors already assigned
	#					remove the already assigned from the cnames and the set of colors
	#					save the new set of assignments
	#return the set of colors in the order of the cnames
	#out will be what is returned
	out = rep(0,times = length(cnames))
	names(out) = cnames
	
	if(is.null(colorPalleteFile)){
		cpalfile = checkFileCopyDefault(fname="./reference_data/uniqueColorPalette.txt")
	}else{
		cpalfile = colorPalleteFile
	}
	cpal = read.table(file=cpalfile,
										stringsAsFactors=F,
										header=F,sep="\t")
	
	colnames(cpal) = c("colorNumber", "colorDescription")
	
	#make dictionary
	assignments = read.table(file=system.file("extdata/colorMatches.txt", package = "packageDir"),
													 stringsAsFactors=F,
													 header=F,
													 sep="\t")
	#clean things up
	assignments[,2]=gsub(pattern=" ", replacement="", x=assignments[,2])
	
	colordict = assignments[,1]
	names(colordict) = assignments[,2]
	
	#pull all assignments out of dictionary
	#find which indexes in out are described in the color lib
	switchVector = cnames%in%names(colordict)
	if(sum(switchVector)){
		#for those that are there, switch them
		out[switchVector] = colordict[names(out)[switchVector]]
	}
	#remove the already chosen from the colordict
	cpalNums = cpal[!cpal[,1]%in%out[switchVector],1]
	
	#removed the already-chosen from the cnames
	cnames2 = cnames[!cnames%in%names(out)[switchVector]]
	#this puts  into cnames
	#for the remaining, assign colors
	out[cnames2] = cpalNums[1:length(cnames2)]
	
	newAssignments = cbind.data.frame(out[cnames2],cnames2)
	rownames(newAssignments) = NULL
	names(assignments) = c("number", "variantType")
	names(newAssignments) = c("number", "variantType")
	outmat = rbind(assignments, newAssignments)
	write.table(x=outmat, sep="\t",
							file="./reference_data/colorMatches.txt", 
							quote=F, 
							row.names=F, 
							col.names=F)

	return(out)
}

# c("missense_variant;NMD_transcript_variant;PolyPhen_probably_damaging",               
# 	"missense_variant;PolyPhen_probably_damaging",                 
# 	"inframe_insertion",
# 	"frameshift_variant;splice_region_variant;intron_variant;feature_elongation",
# 	"frameshift_variant;splice_region_variant;feature_truncation",
# 	"frameshift_variant;feature_elongation",
# 	"frameshift_variant;feature_truncation",
# 	"missense_variant;NMD_transcript_variant;PolyPhen_unknown",          
# 	"missense_variant;NMD_transcript_variant;PolyPhen_possibly_damaging",                
# 	"missense_variant;PolyPhen_possibly_damaging", 
# 	"missense_variant",
# 	"missense_variant;NMD_transcript_variant;PolyPhen_unknown",
# 	"stop_gained;NMD_transcript_variant",
# 	"inframe_deletion",
# 	"splice_donor_variant;feature_truncation",
# 	"frameshift_variant;splice_region_variant;feature_elongation",
# 	"frameshift_variant;NMD_transcript_variant;feature_elongation",
# 	"inframe_insertion",
# 	"stop_gained",
# 	"splice_donor_variant;feature_truncation",
# 	"splice_acceptor_variant;nc_transcript_variant;feature_elongation",
# 	"missense_variant;splice_region_variant;PolyPhen_probably_damaging",
# 	"splice_region_variant;3_prime_UTR_variant;NMD_transcript_variant")

#FilterVariantType
#allows filtering by variant classification
#takes:					ohsuSeqDat: the somatic variants
#								tracker: the tracker object
#								verbose: logical, if input/output should be provided for the user
#								variant_selection: the names of the variant types to keep
#returns: 			list
#										$selected_variants: filtered variant records
#										$tracker: the tracker object
FilterVariantType<-function(ohsuSeqDat, tracker, s, variant_selection = NULL){
	line=""
	verbose=s$interactive
	if(is.null(variant_selection)){
		variant_selection = read.delim(header=F, 
																	 check.names=F,
																	 stringsAsFactors=F,
																	 "./input/AML_corrolated_overlap/defaultOHSUvariantTypes.txt", 
																	 sep=" ")
		variant_selection=as.vector(as.matrix(variant_selection)[1,])
		variant_selection = gsub(pattern=";$", replacement="", x=variant_selection)
	}
	#allCons = parseConsequenceLong(ohsuSeqDat$Consequence) #get all 
	allCons=ohsuSeqDat$Consequence
	stackedGeneBar_OHSUseq(seqDat=ohsuSeqDat)
	
	ohsuSeqDat_sum = summarize_by(col=allCons, 
																display=T, 
																barPlotTitle="Occurances of variant identifiers") #function provides user a summary of the somatic variant data
	cat("For descriptions of the above variants, please visit:\nhttp://uswest.ensembl.org/info/docs/variation/predicted_data.html\n")
	tracker[["variant types and counts"]] = ohsuSeqDat_sum
	#give user the option to filter it: 
	####################### readline() ####################### 
	
	userprompt = "\nEnter row numbers of variant types to be assumed as implying aberrational genes.\n"
	s[[userprompt]] = paste(variant_selection, collapse="; ")
	
	s = settingList(s=s, prompt=userprompt, set=ohsuSeqDat_sum[,1])
	
	selection = s$.text
	
	cat("\nLimiting analysis to lines containing these variant type identifiers:\n")
	print(selection)
	
	
	#selected_variants will contain the set of somatic variants which match the types selected for analysis
# 	glines = rep(F, times=nrow(ohsuSeqDat))
# 	for(vt in selection){
# 		curset = ohsuSeqDat$Consequence%in%vt#grepl(pattern=vt, x=ohsuSeqDat$Consequence)
# 		glines = glines | curset
# 	}
	
	glines = ohsuSeqDat$Consequence%in%selection
	
	selected_variants = ohsuSeqDat[glines,]
	tracker[["Retained variants with these tags as aberrations"]] = paste(selection, sep="; ", collapse="; ")
	tracker[["Ignored these variant tags"]] = paste(ohsuSeqDat_sum$types[!ohsuSeqDat_sum$types%in%selection],sep="; ",collapse="; ")
	cat("\nNow",nrow(selected_variants),"variants remain\n")
	return(list(selected_variants=glines, tracker=tracker, s=s))
}#FilterVariantType

parseConsequence<-function(cons){
	out = c()
	for(i in 1:length(cons)){
		cur = strsplit(x=cons[i], ";")[[1]]
		out = union(out, cur)
	}
	return(out)
}

parseConsequenceLong<-function(cons){
	outcount = 0
	for(i in 1:length(cons)){
		cur = strsplit(x=cons[i], ";")[[1]]
		outcount = outcount + length(cur)
	}
	print("Done counting the number of consequences")
	print(outcount)
	ph=1
	out = rep("", times=outcount)
	
	for(i in 1:length(cons)){
		cur = strsplit(x=cons[i], ";")[[1]]
		out[ph:(ph +length(cur) - 1)] = cur
		ph = ph + length(cur)
	}
	return(out)
}

formatIDs<-function(patIDs){
	
	#if -00 can be found, they're from the OHSU data
	if(length(grep("-00", patIDs[1]))){
		normids = c()
		for(id in patIDs){
			normids = c(normids, paste(strsplit(x=id, split="-00")[[1]], collapse=""))
		}
		return(normids)
	}
	#second, extract Tyner et al patient IDs
	
	if(length(grep("AML", patIDs))){
		#get AML ds patients
		amlpatIDs = patIDs[grep(pattern="AML", x=patIDs)]
		
		normdsid = c()
		for(id in amlpatIDs){
			normdsid = c(normdsid, strsplit(x=id, split="[A-Za-z]+")[[1]][2])
		}
		return(normdsid)
	}	
	print("Could not find OHSU or Tyner pattern patient IDs")
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


# seqfname = "~/tprog/main_131219/input/AML_corrolated_overlap/AMLSeqCapOverlapPatOnly.txt"
# seqcoverage = "~/tprog/main_131219/input/AML_corrolated_overlap/OHSUSeqcoverageForSequenceCapture.txt"
# paths_detail=path_detail
# study_name="amlseqcap"


processSequenceCaptureData<-function(settings, study){
	
	s=settings
	
	tracker = list()
	
	paths_detail=study@studyMetaData@paths
	study_name = study@studyMetaData@studyName

	s  = setting(s=s, prompt="Please select the file of sequence capture, gene variant data..")
	seqfname = s$.text
	
	s  = setting(s=s, prompt="Please select the file containing the set of gene identifiers representing the coverage of the sequence capture analysis..")
	seqcoverage = s$.text

	verbose=s$interactive #set this here to make sure something is in the $interactive slot

	
	covDatRes = loadCoverageData(path_detail=paths_detail, 
															 verbose=s$interactive, 
															 fname=seqcoverage) #"./input/AML_corrolated_overlap/OHSUSeqcoverageForSequenceCapture.txt")
	
	covDat = covDatRes$cov
	tracker[["Coverage data file name"]] =covDatRes$file
	
	covDatm = matrix(data=rep(T, length(covDat)), ncol=1, dimnames=list(covDat))
	
	# 	seqCaptureCoverageSummary = summaryTable4(paths_detail=paths_detail, 
	# 																						enrichment_tests=c(),
	# 																						dataSetName="Sequence Capture Coverage",
	# 																						targetname="sequenced",
	# 																						patientGeneMatrix=covDatm)
	
	ohsuseq = read.table(file=seqfname, 
											 sep="\t", 
											 header=T, 
											 stringsAsFactors=F)
	
	tracker[["Data file name"]] = seqfname
	tracker[["Number of rows found in the input data"]] = nrow(ohsuseq)
	ohsuseq = unique(ohsuseq)
	tracker[["Number of unique rows found in the input data"]] = nrow(ohsuseq)
	
	tracker[["Number of unique sample or patient IDs found in the input data"]] = length(unique(ohsuseq$alias))
	
	preFiltVarPerPatientPreFilt = summarize_by(col=ohsuseq$alias, display=F)
	tracker[["Variants per patient before filtering"]] = paste(names(summary(preFiltVarPerPatientPreFilt[,2])), summary(preFiltVarPerPatientPreFilt[,2]), sep=':', collapse=" ")
	
	while(T){
		s = setting(s=s,prompt="Have manual gene symbol corrections already been made? (y/n)")
		if(s$.text%in%c("y","n")) break
	}
	ohsuseq$Symbol = corsym(symbol_set=ohsuseq$Symbol, symref=paths_detail$symtable, verbose=(s$.text=="n"))
	
	########reduce the sequence data to only that allowed by the coverage set
	ohsuseq = ohsuseq[ohsuseq$Symbol%in%covDat,]
	tracker[["Number of unique rows found in the input data, after limiting to only those containing genes targeted by the sequence capture"]] = nrow(unique(ohsuseq))
	tracker[["Number of unique patient or sample IDs found in the input data,  after limiting to only those containing genes targeted by the sequence capture"]] = length(unique(ohsuseq$alias))
	
	################## readline() ################## 
	snpfiltres = FilterdbSNP(ohsuSeqDat=ohsuseq, s=s, tracker=tracker)
	snpfiltLogicVector = snpfiltres$filtLogicVector
	tracker = snpfiltres$tracker
	s = snpfiltres$s
	#use filtLogicVector to subset the original ohsuseq data
	ohsuseq = ohsuseq[snpfiltLogicVector,]
	
	postdbSNPPerPatientPreFilt = summarize_by(col=ohsuseq$alias, display=F)
	tracker[["Variants per patient after removing those in dbSNP"]] = paste(names(summary(postdbSNPPerPatientPreFilt$counts)), summary(postdbSNPPerPatientPreFilt$counts), sep=':', collapse=" ") 
	######## remove rows with no gene symbol
	ohsuseq = ohsuseq[ohsuseq$Symbol!="\\N",]
	
	######## check for PolyPhen annotation
	ohsuseqppres = CheckIntegratePolyPhen(seqDat=ohsuseq, tracker=tracker)
	tracker=ohsuseqppres$tracker
	ohsuseq = ohsuseqppres$seqDat
	preFiltPolyCount = grep(pattern="polyphen", x=ohsuseq$Consequence, ignore.case=T)
	
	################## readline() ################## 
	varfiltres = FilterVariantType(ohsuSeqDat=ohsuseq, s=s, tracker=tracker)
	varfiltLogicVector = varfiltres$selected_variants
	tracker = varfiltres$tracker
	s = varfiltres$s
	
	#use filtLogicVector and snpfiltLogicVector to subset the original ohsuseq data
	ohsuFiltered = ohsuseq[(varfiltLogicVector),]
	stackedGeneBar_OHSUseq(seqDat=ohsuFiltered)
	cat("\nFiltering complete.. \n")
	postFilterPerPatientVar = summarize_by(col=ohsuFiltered$alias, display=F)
	tracker[["Total unique variants, after all filtering"]] = nrow(ohsuFiltered)
	tracker[["Variants per patient, after all filtering"]] = paste(names(summary(postFilterPerPatientVar$counts)), summary(postFilterPerPatientVar$counts), sep=':', collapse=" ")
	postFiltPolyCount = grep(pattern="polyphen", x=ohsuFiltered$Consequence, ignore.case=T)
	
	tracker[["Variants filtered out with the help of PolyPhen"]] = length(preFiltPolyCount) - length(postFiltPolyCount)
	
	#make df to submit to toPGM
	#subset the columns of ohsuFiltered
	ohsuFilteredSubset = ohsuFiltered[,c("alias","Symbol")]
	colnames(ohsuFilteredSubset)<-c("Ids", "Symbol")
	cat("\nBuilding enrichment function input files\n")
	pgm = toPGM(sds=ohsuFilteredSubset)
	
	# 	tm = getTargetMatrix(tgenes=covDat, paths=paths_detail$paths)
	# 	seqCaptureSummary = summaryTable4(targetname="variant",
	# 																		verbose=verbose,
	# 																		individualEnrichment=T,
	# 																		target_matrix=tm,
	# 																		dataSetName=paste("Sequence capture data,",study_name),
	# 																		patientGeneMatrix=pgm, 
	# 																		paths_detail=paths_detail)
	
# 	tmpStudy = getStudyObject(study.name=study_name, path_detail=paths_detail)

	seqCaptureSummary = summaryTable(study=study,
																	 settings=s,
																	 coverage=rownames(covDatm),
																	 pgm=pgm,
																	 dataSetDescription=paste("Sequence capture data,",study_name),
																	 activeGeneDescription="variant",
																	 coverageDataSetDescription="Sequence Capture Coverage", 
																	 coverageGeneDescription="sequenced")
	
	seqCaptureSummary[["Data_work_up_notes"]] = tracker
	# 	seqCaptureSummary$coverage_summary = seqCaptureCoverageSummary
	
	return(seqCaptureSummary)
}#processSequenceCaptureData()


