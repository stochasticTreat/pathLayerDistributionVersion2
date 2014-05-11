#somatic_mutations_processing_v8.R

#runSomaticMutationsProcessing
#main interface for somatic mutations input arm
#implements armMain interface, thus this function is passed as armMain() inside the runArm() function
#	Takes: 	settings: a settings object
#					study: a study object (defined in initiateDataStructures.R)
# Returns: somatic mutations summary list of lists
runSomaticMutationsProcessing<-function(settings, study){
	
	if(settings$interactive){
		print("Running somatic mutations processing in interactive mode")
	}else{
		print("Running somatic mutations processing from saved settings")
	}
	#initialize objects
	path_detail = FullPathObject(S=study)
	somatic_summary = list() #this is the main summary object that is provided in the main name space for the main function
	s = settings
	interactiveTMP = s$interactive #this allows that, if a user enters something the program does not understand, 
		#the program will cycle back through and re-prompt the user and the settings will not be screwed up.
	############################################################
	##get somatic mutation files and extract, map mutations
	############################################################
	cat("\nMain interface for processing somatic mutation data\n")
	
	###selection structure
	#load chasm data/produce data set to send to CHASM
	#load somatic mutation data
	#quit
	while(T){
		
		p1 = "To load somatic mutation data: enter 1\nTo process somatic sequencing data, limiting the coverage, enter 2\nTo exit somatic mutation interface: enter 3.\n"
		#################### readline() ####################################
# 		uin = readline(prompt=p1)
		s = setting(s=s, prompt=p1)
		uin = s$.text
		
		if(uin=="1"){
			somatic_summary = processSomaticData(paths_detail=path_detail, 
																					 study=study,
																					 s=s,
																					 study_name=study@studyMetaData@studyName)
			s$interactive = interactiveTMP #not sure why I did this... 
			return(somatic_summary)
		}else if(uin=="2"){
# 			if(!exists("study_name")){
# 				study_name = "running_only_somatic_data_input"
# 			}
# 			htmlFname = paste("./output/",
# 												study_name,
# 												"_somatic_mutation_summary.html",
# 												sep="")
# 			htmlSummary(sumset=somatic_summary, fname=htmlFname)
		}else if(uin=="2"){
			source('./processSomaticSeqDataWithCoverage.R')
			somatic_summary = processSomaticDataWithCoverage(paths_detail=path_detail, 
																											 verbose=T, 
																											 study_name=study_name)
		}else if(uin=="3"){
			break
		}
		s$interactive = TRUE
	}
}

top20Hists<-function(unfilteredData,fileroot){


# 	unfilteredData=results$somatic_mutation_aberration_summary$unfiltered_data
	
	#get the top 20 mutations
	stacked = unfilteredData[,c("pid","Hugo_Symbol")]
	names(stacked)<-c("Ids", "Symbol")
	prepgm = toPGM(sds=stacked)
	gsum = prepgm%*%rep(1, ncol(prepgm))
	
	#now order the gsum: 
	gsum = gsum[order(gsum, decreasing=T),,drop=F]
	
	top30 = gsum[1:min(30,nrow(gsum)),,drop=F]
	
	for(i in 1:30){
		cn = rownames(top30)[i]
# 		start_positions = unfilteredData$Start_Position[unfilteredData$Hugo_Symbol==cn]#find all the start positions with Symbol == cn
		print(cn)		
		curgenmat = unfilteredData[unfilteredData$Hugo_Symbol==cn,,drop=F]
		postype2 = curgenmat
		postype2$Variant_Classification<-as.factor(postype2$Variant_Classification)
		bw = range(curgenmat$Start_Position)[2]-range(curgenmat$Start_Position)[1]
		bw = bw/75

		qp = qplot(Start_Position, 
							 data = postype2,  
							 binwidth=bw,
							 fill = Variant_Classification, 
							 main=paste("Positions and distributions of variants in gene", cn))
		qp + theme(axis.text.x=element_text(angle=-90))
		dir.create("./output/image_tmp/",recursive=T, showWarnings=F)
		ggsave(filename=paste("./output/image_tmp/variant_postions_and_types_for_gene_",cn,".png",sep=""), plot=qp)
	}
}

# stackedGeneHist<-function(postype){
# 	
# 	barDf = matrix(data=0, nrow=nrow(postype), ncol=length(unique(postype$Variant_Classification)))
# 	xlims=range(postype$Start_Position)
# 	uclasses = unique(postype$Variant_Classification)
# 	
# 	barindex = 1
# 	for(ci in 1:length(uclasses){
# 		type = uclasses[ci]
# 		#for each variant type
# 		varrows = postype$Variant_Classification==type
# 		cat(" ",sum(varrows)," ")
# 		#get all the start positions
# 		starts = postype$Start_Position[varrows]
# 		#get the histogram
# 		chist = hist(x=starts,
# 								 breaks=100)
# 		
# 		cb = chist$breaks
# 		
# 		chist$vartype = type
# 		
# 	}
# }


matchTabRow<-function(trow,tab){
	cols = colnames(tab)
	out = rep(T, nrow(tab))
	for(c in 1:ncol(tab)){
		lv = tab[,c]%in%trow[,c]
		out = out&lv
	}
	return(out)
}

getColorSequence<-function(cnames, 
													 fname="./reference_data/MAFcolorMatches.txt"){
	#first check cnames, 
	#			if colors already assigned to cnames, 
	#					set those colors already assigned
	#					remove the already assigned from the cnames and the set of colors
	#					save the new set of assignments
	#return the set of colors in the order of the cnames
	#out will be what is returned
	
	fname = checkFileCopyDefault(fname=fname)
	
	out = rep(0,times = length(cnames))
	names(out) = cnames
	
	cpal = read.table(file=system.file("extdata/uniqueColorPalette.txt", package = "packageDir"),
										stringsAsFactors=F,
										header=F,sep="\t")
	colnames(cpal) = c("colorNumber", "colorDescription")
	
	#make dictionary
	assignments = read.table(file=fname,
													 stringsAsFactors=F,
													 header=F,
													 sep="\t")

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
	write.table(x=outmat, 
							sep="\t",
							file=fname, 
							quote=F, 
							row.names=F, 
							col.names=F)
	return(out)
}#getColorSequence

stackedGeneBar<-function(tcga_som){
	print("***************creating stacked Gene Barplot***************")
	#stackedGeneBar
	#takes: tcga_som: initial input .maf table, after cleaning of repeat rows and gene symbols
	
	ufilt = tcga_som

	gs = summarize_by(col=tcga_som[,"Hugo_Symbol"], display=F)
	top20 = head(gs[order(gs[,2],decreasing=T),],20)
	gsnames = top20$types
	
	print("filtering")
	toprows = ufilt[ufilt$Hugo_Symbol%in%gsnames,]
	
	slimrows = toprows[,c("Hugo_Symbol","Variant_Classification")]

	#set up the output matrix
	mtypes = unique(slimrows$Variant_Classification)
	mmat = matrix(data=0,ncol=length(mtypes), nrow=length(gsnames), dimnames=list(gsnames, mtypes))

	for(n in gsnames){
		print(n)
		#get the rows for those names
		msum = summarize_by(col=slimrows[slimrows$Hugo_Symbol==n,"Variant_Classification"], display=F)
		mmat[n,msum$types] = msum$counts
	}
	
	mdf  = cbind.data.frame(gene = rownames(mmat),mmat)
	oldmar = par()$mar
	oldlas = par()$las
	par(las=2) # make label text perpendicular to axis
	par(mar=c(5,8,4,2)) # increase y-axis margin.	
	scol = 2
	colcols = getColorSequence(cnames=colnames(mmat), fname="./reference_data/MAFcolorMatches.txt")
	barplot(t(mmat), legend = colnames(mmat),xlab="Number of mutations found in gene",
					main="Mutation types for the top 20 most mutated genes", 
					col=colors()[colcols],
					horiz=TRUE, 
					cex.names=0.8)
	par(las=oldlas)
	par(mar = oldmar)
	#barplot(mmat, horiz=T)
}

#remove_dbSNP
#removes mutations with dbSNP records
#allows output and user input
#returns:  list$som_select : records with dbSNP values optionally removed
#							 $tracker			: the work-up tracker

remove_dbSNP<-function(tcga_som, tracker, s){
	dbsnp_set  = tcga_som[,"Dbsnp_Rs"]!=""
	tracker[["Number of variants found in dbSNP"]] = sum(dbsnp_set)
	if(sum(dbsnp_set)){
		cat("\n", as.integer(sum(dbsnp_set)), "of the somatic mutations were found to have dbSNP records.")
		cat("\nThese records are of types:")
		dbsnpstatus= unique(tcga_som[dbsnp_set,"Dbsnp_Val_Status"])
		statsum = summarize_by(col=tcga_som[,"Dbsnp_Val_Status"], display=F)
		colnames(statsum)<-c("dbSNP val status types","Number of records with type")
		statsum=as.matrix(statsum, quote=T)
		print(statsum)

		s = setting(s=s, prompt="Would you like to filter out those somatic mutations with dbSNP values? (please enter y or n) ")
		
		removedbSNP = s$.text
		if(sum(removedbSNP=="y")){
			cat("\nMutations with dbSNP records will be filtered out\n")
			tcga_som = tcga_som[!dbsnp_set,]
			cat("\nAfter filtering out dbSNP mutations,",nrow(tcga_som),"unique mutations remain\n")
			tracker[["Number of mutations after filtering by dbSNP record"]] = nrow(tcga_som)
		}
	}
	return(list(tracker=tracker, som_select=tcga_som, s=s))
}#remove_dbSNP

#filterMutationType
#allows filtering by mutation classification
#takes:					tcga_som: the somatic mutations
#								tracker: the tracker object
#								verbose: logical, if input/output should be provided for the user
#returns: 			list
#										$som_select: filtered mutation records
#										$tracker: the tracker object
filterMutationType<-function(tcga_som, tracker, s){

	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	print("Before removal of genes marked as \"UNKNOWN\"")
	stackedGeneBar(tcga_som)
	tracker[["Before removal of UNKNOWN genes, distribution of mutation types in top 20 most mutated genes"]]=save.plot("stackedGeneBarPreUnknownRemoval")
	
	tcga_som_no_unknown = tcga_som[tcga_som$Hugo_Symbol!="UNKNOWN",]
	print("After removal of genes marked as \"UNKNOWN\"")
	stackedGeneBar(tcga_som_no_unknown)#make another stacked gene bar after the "UNKNOWN" are removed
	tracker[["After removal of UNKNOWN genes, distribution of mutation types in top 20 most mutated genes."]]=save.plot("stackedGeneBarPostUnknownRemoval")

	tcga_som_sum = summarize_by(col=tcga_som[["Variant_Classification"]], left_margin_factor=1.3,
															display=T, 
															barPlotTitle="Counts of different types of somatic mutations") #function provides user a summary of the somatic mutation data
	tracker[["Distribution of variant types for all genes"]] = save.plot("Distribution of variant types for all genes")
	tracker[["Mutation types and counts"]] = tcga_som_sum
	#give user the option to filter it: 
	####################### readline() ####################### 
	cat("\n\n")
	prompt = "\nPlease enter the row numbers of the variant types you would like to analyze (sepparated by a space).\n"
	s = settingList(s=s, prompt=prompt, set=tcga_som_sum)
	
	selection = s$.text
	
	cat("\nLimiting analysis to these somatic mutation types:\n")
	print(selection)
	#som_select will contain the set of somatic mutations which match the types selected for analysis
	som_select = tcga_som[tcga_som[,"Variant_Classification"]%in%selection,]
	tracker[["Retained these mutation types"]] = paste(selection, sep="; ", collapse="; ")
	tracker[["Filtered out these mutation types"]] = paste(tcga_som_sum$types[!tcga_som_sum$types%in%selection],sep="; ",collapse="; ")

	cat("\nNow",nrow(som_select),"mutations remain\n")
	tracker[["Number of mutations after filtering by mutation type"]] = nrow(som_select)
	return(list(som_select=som_select, tracker=tracker, s=s))
}#filterMutationType

#makePatientGeneMatrix
#takes: 	somatic mutation records
#returns: patient gene matrix: 
#															rows: genes
#															columns: patients
#															values: logical, indicating if gene in patient is mutated
makePatientGeneMatrix<-function(som_select){
	unique_genes = unique(som_select[,"Hugo_Symbol"])
	pid = unique(som_select[,"pid"])
	cat("\n Building patient gene matrix")
	###### make patient gene matrix: patients = columns; genes = rows

	somatic_pgm = matrix(0, ncol=length(pid), nrow=length(unique_genes))

	rownames(somatic_pgm)<-unique_genes
	colnames(somatic_pgm)<-pid
	for(i in 1:nrow(somatic_pgm)){
		#find what genes that patient has mutated, then stick them in to the matrix
		curgene = as.character(unique_genes[i])
		curpats = as.character(som_select[som_select[,"Hugo_Symbol"] == curgene,"pid"])
		somatic_pgm[curgene,curpats] = 1
	}
	return(somatic_pgm)
}#makePatientGeneMatrix


FilterDuplicates<-function(tcga_som_raw, s, tracker, paths_detail){
	
	########################################################
	############### 									 Filtering duplicates
	
	cat("\n\nFiltering and processing mutation data... . \n")
	#find names for all somatic mutation data
	cat("Removing rows whose values for all columns are identical to another row.\n")
	
	tcga_som = tcga_som_raw[!duplicated(tcga_som_raw),]
	cat("\n",nrow(tcga_som),"rows remain\n")
	tracker[["Number of records after removing duplicate rows"]] = nrow(tcga_som)	
	cat("\nRemoving rows for which the values in the following colulmns are identical to those in other rows:\n")
	filtcolumns =c( "Hugo_Symbol", "Chrom","Ncbi_Build", "Start_Position", "Reference_Allele", "Tumor_Sample_Barcode" )
	filtcolumns2 =c( "Hugo_Symbol", "Chrom","Ncbi_Build", "Start_Position", "Reference_Allele", "Tumor_Sample_Barcode", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
	print(filtcolumns2)
	
	utest1 = unique(tcga_som[,filtcolumns])
	utest2 = unique(tcga_som[,filtcolumns2])
	cat("\n",nrow(utest2), "rows remains\n")
	
	tcga_som = tcga_som_raw[!duplicated(tcga_som_raw[,filtcolumns2]),]
	tracker[["Number of records after removing duplicate mutations"]] = nrow(tcga_som)
	tracker[[" Note: removal of duplicate mutations assumes these columns can serve as a unique key to identify the mutation."]]  = paste(filtcolumns2,sep=" ",collapse = " ")
	cat("\nReplicated rows in original data were likely technical replicates from sequencing on different platforms.\n")
	dupcol = c("Hugo_Symbol", "Chrom", "Start_Position", "Tumor_Sample_Barcode" )
	#dup = duplicated(utest1[,c("Hugo_Symbol", "Chrom", "Start_Position", "Tumor_Sample_Barcode" )])
	dup = duplicated(utest1[,dupcol])
	if(sum(dup)){
		print(tcga_som[dup,])
		cat("\n\n!!Please check your data, based on the values in these columns:\n")
		print(dup)
		cat("\nThere appears to be an issue where the base in question on the\n",
				"reference allele is different between the two matched normal samples.\n",
				"This may lead to spurious results.")
		blnk = readline("Press enter to continue.")
	}
	if(nrow(utest1)!=nrow(utest2)){
		#################### readline() ####################################
# 		blnk = readline("!!PLEASE NOTE: it appears at least one mutation is inconsistent\nbetween technical replicates (sequencing on different platforms).\nPlease press enter to continue with analysis.")
	
		s = setting(s=s, prompt="!!PLEASE NOTE: it appears at least one mutation is inconsistent\nbetween technical replicates (sequencing on different platforms).\nPlease press enter to continue with analysis.")
		
	}
	########################################################
	###### Correct gene symbols
	########################################################
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	s = setting(s=s, prompt="Have manual gene symbol corrections already been conducted? (y/n)")
	customCorrections = s$.text=="n"
	if(s$.text == "n") s[["Have manual gene symbol corrections already been conducted? (y/n)"]] = NULL
	if(customCorrections){
		s = setting(s=s, prompt="Run manual symbol corrections? (y/n)")
		
		customCorrections = s$.text=="y"
	}
	tcga_som[,"Hugo_Symbol"] = corsym(tcga_som[,c("Hugo_Symbol", "Chrom")],
																		symref=paths_detail$symtable, 
																		verbose=customCorrections)
	#####check again for duplicates
	minimalUniqueKey = c("Chrom", "Start_Position", "Reference_Allele", "Tumor_Sample_Barcode", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
	
	#check for problems
	#####################################################################################
	############### 	re-filtering for duplicates to clear up gene name annotation issues
	dup_index_forward = duplicated(tcga_som[,minimalUniqueKey])
	dup_index_backward = duplicated(tcga_som[,minimalUniqueKey], fromLast=T)
	alldupli = dup_index_backward | dup_index_forward
	
	tcga_som_notdup = tcga_som[!alldupli,]
	tcga_som_dup = tcga_som[alldupli,]
	
	#now scan the duplicated to remove any associated with an old sequencer
	orderedDuplicates = tcga_som_dup[order(tcga_som_dup$Start_Position),]
	
	sequencerTypes = unique(tcga_som[,"Sequencer"])
	if(length(sequencerTypes)>1){
		cat("\n\n")
		sequencerTypes = data.frame(sequencerTypes, stringsAsFactors=F)
		print(sequencerTypes)
		#################### readline() ####################################
		prompt = "It appears more than one sequencer type was used and that this has caused\nsome individual mutations to be annotated with outdated gene symbols.\nPlease select the number of newest model, with the best annotation from the list above: "
		s = settingList(s=s, prompt=prompt, set=sequencerTypes)
		
		seqselection = s$.text
		
		#go through the duplicated and remove the records using the old sequencer
		badindexes = c()
		for(i in 1:nrow(tcga_som_dup)){
			cur = tcga_som_dup[i,]#get one out
			#find its duplicate
			di = which(matchTabRow(trow=cur[,minimalUniqueKey],
														 tab=tcga_som_dup[,minimalUniqueKey]))
			#determine if one doesn't use the preferrable sequencer and throw it out 
			if(seqselection%in%tcga_som_dup$Sequencer[di]){#if one of the duplicates used the better sequencer
				badindexes = union(badindexes,di[!tcga_som_dup$Sequencer[di]%in%seqselection]) #find the index(es) of those that used the worse sequencer
			}
		}	
		tcga_som_dup_rem = tcga_som_dup[!1:nrow(tcga_som_dup)%in%badindexes,]
		#now merge the filtered duplicates with the main set:
		tcga_som = rbind(tcga_som_notdup,tcga_som_dup_rem)
	}
	
	#now check that tcga_som now has the same number of unique rows as with the minimal key
	errorCheckUniqueCount = nrow(unique(tcga_som[,minimalUniqueKey]))
	if(errorCheckUniqueCount!=nrow(tcga_som)){
		#################### readline() ####################################
		prompt = "!Allert, after repeated attempts, data can not be coorsed to remove unique records,\nplease check data to avoid spurioius results."
		s = setting(s=s, prompt=prompt)
	}

	return(list(tcga_som=tcga_som, tracker=tracker, s=s))
}#FilterDuplicates()


#maf = tcga_som_raw
addDbSNPToVariantClassification<-function(maf){
	#first find all the rows with dbSNP values
	dbSNPi = maf$Dbsnp_Rs!=""
	#for those rows with dbSNP vals, append a note to the variant type
	maf$Variant_Classification[dbSNPi] = paste(maf$Variant_Classification[dbSNPi], "In_dbSNP", sep="_")
	return(maf)
}

#processSomaticData
#takes:					paths_detail: a path list object
#								removedbSNP: logical, T if dbSNP values should be removed
#								fname:			 The name of the file containing the somatic mutation data to be processed
#								verbose:		 logical, if true gives extra output
#								study_name: 	the study name included, used as a file name prefix
processSomaticData<-function(study, 
														 study_name="no_study_name_provided",
														 paths_detail,
														 s=list()){

	fname_base = study_name
	tracker = list()
	# 	tracker[["Study name"]] = study_name
	
	########################################################
	############### Open the file and have a look ##########
	#################### readline() ########################
	########################################################
	s = setting(s=s, prompt="Select a .maf file containing the data set to be analyzed.")
	
	fname = s$.text
	
	fname = checkFileCopyDefault(fname=fname)
	
	cat("\nLoading file:", fname,"\n")
	tracker[["Data file name"]] = fname
	
	tcga_som_raw = read.delim(file=fname, header=T, sep="\t", stringsAsFactors=F, na.strings="-")
	
	cat("\n",nrow(tcga_som_raw),"rows read in from file.")
	cat("\nFixing patient ids...\n")
	pid = sapply(tcga_som_raw[,"Tumor_Sample_Barcode"], extract_pid)
	tcga_som_raw = cbind.data.frame(pid, tcga_som_raw, stringsAsFactors=F)#append extracted pids as a sepparated column
	
	cat("\n",length(unique(pid)),"unique patient sample id's found in input file.\n")
	tracker[["Unique patient IDs found in file"]] = length(unique(pid))
	tracker[["Rows of data read in from file"]] = nrow(tcga_som_raw)
	
	hgbuild = unique(tcga_som_raw$Ncbi_Build)
	cat("\nThe build of the human genome used to annotate the sequencing data:", hgbuild)
	tracker[["The build of the human genome used to annotate the sequencing data:"]] = hgbuild
	
	#################### readline() ####################################
	if(length(hgbuild)>1){
		readline("!!Warning, it appears the mutation records in the input sequencing\ndata were annotated using more than one unique build of the human geneome.\nThis may lead to spurious results!")
	}
	
	##########   process, clean the table
	
	#################### readline() ####################################
	filtRes = FilterDuplicates(tcga_som_raw=tcga_som_raw, 
														 tracker=tracker,
														 s=s,
														 paths_detail=paths_detail)
	tcga_som = filtRes$tcga_som
	tracker = filtRes$tracker
	s=filtRes$s
	
	s = setting(s=s, prompt="Would you like to append dbSNP status to variant classifications? (y/n)")
	if(s$.text=="y"){
		tcga_som = addDbSNPToVariantClassification(maf=tcga_som)
	}
	
	###### check number of genes that then have official HUGO symbols
	if("Status"%in%paths_detail$symtable){
		approvedHugoSymbols = paths_detail$symtable$Approved.Symbol[paths_detail$symtable$Status == "Approved"]
	}else{
		approvedHugoSymbols = unique(paths_detail$symtable$Approved.Symbol)
	}
	
	notApprovedHugoVector = tcga_som[!tcga_som[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols),"Hugo_Symbol"]
	numNotHugo = sum(!tcga_som[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols))
	cat("\n",numNotHugo," unique mutations do not have official HUGO symbols associated\n",sep="")
	tracker[["After symbol correction, the number of mutations without approved HUGO symbols:"]] = numNotHugo
	approvedHugoSymbols = paths_detail$symtable$Approved.Symbol
	tracker[["List of all symbols from the input that were not official HUGO symbols:"]] = matrix(data=tcga_som[!tcga_som[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols),"Hugo_Symbol"], ncol=1)
	
	#########################################################
	###############  Summarize patient somatic mutation stats
  subset = unique(pid)
	cat("\nPatient sample identifers loaded:\n")
 	print(subset)
	####extract the patient subset if a subset was established

	cat("\nIn this set of patients", as.character(nrow(tcga_som)), "unique somatic mutations were found.\n")
	
	####################### readline() ####################### 
	###### allow PolyPhen adjustment of data
	polyres = addPolyPhenResults(mafData=tcga_som, tracker=tracker, s=s)
	tcga_som = polyres$mafData
	tracker = polyres$tracker
	s = polyres$s
	
	cat(".")
	preFilteringGeneSummary = summarize_by(col=tcga_som[,"Hugo_Symbol"], display=F)
	#preFilteringWithHists = top20Hists(preFilteringGeneSummary)
	tracker[["Mutations per gene before filtering"]] = preFilteringGeneSummary
	
	preFiltCountByPatient = summarize_by(col=tcga_som$pid, display=F)
	pnames = preFiltCountByPatient$types
	preFiltCountByPatient = matrix(preFiltCountByPatient$counts, dimnames=list(pnames))
	
	tracker[["Mutations per patient, before any filtering"]] = as.matrix(preFiltCountByPatient)
	tracker[["Summary of mutations per patient before filtering"]] = paste(summary(preFiltCountByPatient), collapse="", sep="")

	###########################################################
	##########     Mutation type filtering
	
	####################### readline() ####################### 
	filtered = filterMutationType(tcga_som=tcga_som, tracker=tracker, s=s)
	som_select = filtered$som_select
	tracker = filtered$tracker
	s=filtered$s
	
	###############################################################
	###################   dbSNP filtering

	SNPless = remove_dbSNP(tcga_som=som_select, tracker=tracker, s=s)
	som_select = SNPless$som_select
	tracker = SNPless$tracker
	s=SNPless$s
	
	###### check number of genes that then have official HUGO symbols
	approvedHugoSymbols = paths_detail$symtable$Approved.Symbol[paths_detail$symtable$Status == "Approved"]
	notApprovedHugoVector = som_select[!som_select[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols),"Hugo_Symbol"]
	numNotHugo = sum(!som_select[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols))
	cat("\nAfter filtering, there are ",numNotHugo," unique mutations that do not have official HUGO symbols.\n",sep="")
	tracker[["Number of mutations without approved HUGO symbols (note: this was checked after filtering and symbol correction)"]] = numNotHugo
	
	##############################
	###### Filtering out mutations with HUGO symbols given as "UNKNOWN"
	
	som_select = som_select[som_select[,"Hugo_Symbol"]!="UNKNOWN",]
	tracker[["Number of unique mutations after removing records with gene identifier given as \"UNKNOWN\""]] = nrow(som_select)
	
	tracker[["Number of patients with mutations not removed by the filtering steps:"]] = length(unique(som_select$Tumor_Sample_Barcode))
	############################################################
	###### build output           ##############################
	############################################################
	uniquePatientSymbolsMutated= nrow(unique(som_select[,c("pid","Hugo_Symbol")]))
	uniquePatientSymbolsMutated2= nrow(unique(som_select[,c("Tumor_Sample_Barcode","Hugo_Symbol")]))
	uniquePatientMutations = nrow(som_select)
	tracker[["Number of genes across cohort that are mutated more than once in the same patient"]] = uniquePatientMutations - uniquePatientSymbolsMutated
	tracker[["Genes are considered to be in a mutated or normal state, thus the final number of unique, mutated genes passed to the enrichment analysis is:"]] = uniquePatientSymbolsMutated
	#########################################
	# create the patient gene matrix. 
	#########################################
	somatic_pgm = makePatientGeneMatrix(som_select)

	countByPatient = t(rep(T,nrow(somatic_pgm))%*%somatic_pgm)
	tracker[["Mutation counts by patient, after all filtering, for data set used in pathway analysis"]] = countByPatient
	tracker[["Summary statistics for mutations per patient after all filtering (data set used in pathway analysis)"]] = paste(summary(countByPatient), collapse="", sep="")
	tracker[["Number of mutations in final, cleaned data set used for pathway enrichment"]] = sum(somatic_pgm)
	table_out_fname = paste("./output/",study_name,"_summary_of_somatic_mut_by_pathway.txt",collapse="",sep="")

	somatic_summary = summaryTable(study=study,
																 pgm=somatic_pgm, 
																 originalDataMatrix=tcga_som, 
																 activeGeneDescription="mutated", 
																 dataSetDescription="somatic mutation data", 
																 settings=s)
	
	prefilt = tracker[['Mutations per patient, before any filtering']]
	postfilt = somatic_summary$patientsums
	twoHistOnePlot(dataset1=prefilt, dataset2=postfilt, 
								 x_label="Number of mutations", 
								 y_label="Number of patients",
								 legend_titles=c("Before filtering", "After filtering"),
								 main_title="Distributions of mutations in patients before and after filtering")

	tracker[["Distributions of mutations before and after fitering"]]=save.plot("DistMutBeforeAferFilt")
	somatic_summary[["Data_work_up_notes"]] = tracker
	somatic_summary[["preFilteringGeneSummary"]] = preFilteringGeneSummary
	#	print(names(somatic_summary$settings))
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	return(somatic_summary)
}#processSomaticData()


