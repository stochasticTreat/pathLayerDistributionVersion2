#data in: patientGeneLevels (patient gene matrix with continuous values for each gene)



#main function for runing drug screen
#implements armMain interface, thus this function is passed as armMain() inside the runArm() function
RunDrugScreen<-function(settings, study){
	
	#initialize objects
	path_detail = getPaths(study)
	drug_screen_summary = list() #this is the main summary object that is provided in the main name space for the main function
	s = settings
	interactiveTMP = s$interactive #this allows that, if a user enters something the program does not understand, 
																	#the program will cycle back through and re-prompt the user and the settings will not be screwed up.
	while(T){
		s = setting(s=s, prompt=paste("\n To analyze drug screen panel coverage (for a panel that has or has not been run), enter p\n",
																	"To process drug screen result set enter d\n",
																	"To save an HTML summary of the results enter h\n",
																	"To exit drug screen interface, enter q\n"))
		s$interactive=interactiveTMP
		mainsel = s$.text
		if(mainsel=="p"){
			drug_screen_summary = getPanelCoverage(path_detail=path_detail, s=s, study=study)
			return(drug_screen_summary)
		}else if(mainsel=="d"){
			#setwd("~/tprog/drug_screen_arm")
			drug_screen_summary = runPanelAnalysis(path_detail=path_detail, s=s, study=study)
			return(drug_screen_summary)
		}else if(mainsel=="q"){
			if(!length(drug_screen_summary)){return(NULL) 
			}else{return(drug_screen_summary)}
		}
		s$interactive=T
	}#while
	s$interactive = interactiveTMP
}

#isStackedFormat
#checks if one of the columns has the string "score", and that there are only 3 columns, 
#indicating data is in stacked format
#
#takes: fname: file name of drug screen results file
#returns: logical, indicating if data is in stacked format
isStackedFormat<-function(fname){
	cat("\nNote: if data is provided in stacked format, it must contain 3 columns, with the third containing the substring \"score\"\n")
	dstab = read.table(file=fname, header=T, comment.char="", sep="\t", stringsAsFactors=F,nrows=4)
	if(sum(grepl(pattern="score",x=colnames(dstab),ignore.case=T))>0 & ncol(dstab)==3) return(T)
	return(F)
}


#getTargetGenesFromDrugs
#unpacks drug targets from drug-gene target matrix
#examines which drugs are present on panel to determine which genes are targeted by panel
#takes: drugTargetsFileName: path to file containing drug target matrix: bipartate graph with drug names on one axis 
#				panelFileName : path to panel file (must have column <Drugs> containing drug names)
#returns: 
getTargetGenesFromDrugs<-function(drugTargetsFileName="./input/drug_targets_2.txt", panelFileName="./input/exp006_extract_for_db.txt",verbose=F){
	
	cat("\nCurrently this procedure only provides a list of genes which are targeted,",
			"\nbut does not give the ihibitory magnitude.\n")
	
	drug_targets <- read.delim(drugTargetsFileName, header=F, quote="",na.strings="-")
	drugs = drug_targets[1,]
	genes = as.character(drug_targets[,1])
	sysnames = sapply(genes, function(x){strsplit(x, "\\(")[[1]][1]}  )
	recs = (nrow(drug_targets)-1) * (ncol(drug_targets)-1)
	
	drug = rep("blank", recs)
	tgene = rep("blank", recs)
	sysname = rep("blank", recs)
	degree = rep(-1, recs)
	
	c = 0
	for(i in 2:nrow(drug_targets))#for each gene
	{
		for(j in 2:ncol(drug_targets))#for each drug
		{
			c = c + 1
			drug[c] = as.character(drug_targets[1,j])
			tgene[c] = genes[i]
			sysname[c] = sysnames[i]
			degree[c] = as.integer(drug_targets[i,j]) -1
		}
	}
	drug_gene_xref = data.frame(cbind(drug=toupper(drug), tgene=toupper(tgene), sysname=toupper(sysname), degree=degree), stringsAsFactors=F)
	mds = read.delim(panelFileName, header=T, sep="\t", stringsAsFactors=F,na.strings="-")
	paneldrugs = toupper(mds[,"Drug"]) #there are 158 unique names here,but some of them are replicates
	
	#now see which rows of out have drug names that are in dkmdrugs
	drugsinboth = drug_gene_xref[drug_gene_xref$drug%in%paneldrugs,]
	cat("\n",length(unique(drugsinboth[,"drug"])), " drugs from drug panel file match the drugs in the drug target matrix.\n",sep="")
	if(length(unique(drugsinboth[,"drug"]))!=length(unique(paneldrugs))){
		cat("\nThese are the drugs from the panel that were not found in the drug-target matrix:\n")
		notFoundSet = setdiff(paneldrugs, drug_gene_xref$drug)
		print(notFoundSet)
		cat("Some of these may be replicates (typically ending in the pattern _<replicate number>). ",
				"\nthose that are not replicates may have variations in notation or may need to be added",
				"\nto the target matrix file.")
		fnameDrugsNotFound="./output/drug_names_not_found_in_target_matrix.txt"
		cat("\n\nThe list of drugs not found in the target matrix was saved to", fnameDrugsNotFound,"\n")
		write.table(x=notFoundSet, file=fnameDrugsNotFound,quote=F, sep="\t")
	}
	#now get unique genes from dkmrows
	ugenes = unique(drugsinboth[,"sysname"])
	if(verbose){readline("\nPress enter to continue")}
	return(ugenes)
}

#cleans symbols, removing -phosphorylated and items in parenthesis
#takes: data frame with column titled "gene"
#returns: the input data frame with symbols in "gene" column cleaned
cleanSyms<-function(tab, s, symtab=NULL){
	cat("in cleanSyms()\n")
	if(is.null(symtab)){
		symtab = getHugoSymbols()
	}
	#remove -domain information
	allSymsPlusPhos = sapply(X=tab$gene,function(x){strsplit(x=x,split="\\(")[[1]][1]})
	#remove phosphorylation info
	allSyms = sapply(X=allSymsPlusPhos,function(x){strsplit(x=x,split="-phosphorylated")[[1]][1]})
	
	cat("\nIn the input file,",length(allSyms),"unique gene symbols were found")
	s = setting(s=s, prompt="Have manual gene symbol corrections already been made? (y/n) ")
	
	allSymsCor = corsym(symbol_set=allSyms, 
											symref=symtab, 
											verbose=s$.text=="n")
	tab$gene = allSymsCor
	return(list(txt=tab, s=s))
}

#takes drug screen file in stacked format and returns a patien scores matrix
DrugScreenFromStackedFormat<-function(s,fname=NULL){
	cat("in DrugScreenFromStackedFormat()\n")
	if(is.null(fname)) fname = file.choose()
	dstab = read.table(file=fname,header=T,comment.char="", sep="\t", stringsAsFactors=F)
	colnames(dstab) = c("SampleID","gene","ScoreValue")
	res = DrugMultiScoreReduce(dstab=dstab, s=s)
	pgm = res$txt
	s=res$s
	
	print(dim(pgm))
	print("returning from DrugScreenFromStackedFormat")
	return(list(txt=pgm, s=s))
}

#DrugMultiScoreReduce
#
# takes drug screen data in stacked format and produces a patient gene score matrix of drug screen scores. 
# If patient-gene combos are found to have multiple scores associated, the highest score is used
#
#
#takes: dstab: drug screen data in stacked format
#								must have columns with these names: SampleID, gene, ScoreValue    !!!
#
DrugMultiScoreReduce<-function(dstab, s){
	cat("in DrugMultiScoreReduce()\n")
	#remove blank rows: 
	blanki = which(dstab$gene=="")
	dstab = dstab[!1:nrow(dstab)%in%blanki,]
	#first get the unique set of patients:
	head(dstab)
	upat = unique(dstab$SampleID)
	
	#now get the table with the cleaned up gene names:
	res = cleanSyms(dstab, s)
	cdstab = res$txt
	s=res$s
	
	upat = unique(cdstab$SampleID)
	usym = unique(cdstab$gene)
	
	pgm = matrix(data=NA, nrow=length(usym), ncol=length(upat), dimnames=list(usym, upat))
	
	for(p in upat){#for each patient
		psub = cdstab[cdstab$SampleID==p,]#get the set of rows for that patient
		for(sym in usym){#for each symbol annotated to that patient
			rawscores = psub[psub$gene==sym,"ScoreValue"]#get the set of rows in that patient with that symbol
			scores = as.numeric(rawscores)
			if(sum(is.na(scores))){
				cat("\nThe drug screen score value(s) \"",rawscores[is.na(scores)],"\" for the  gene \"",sym,"\" in patient \"",p,"\"\ncould not be coerced into a numeric value.\n",sep="")
				cat("Here are all the scores found for the above mentioned patient/gene combo:\n")
				cat(rawscores)
				stop("Error, not all scores could be coerced into numeric values (see above).\nPlease check the drug file and re-run program.")
			}
			if(length(scores)>1){
				if(sum(scores>0)>1){
					cat("Multiple scorecs found for patient", p, ", gene", sym, ". These are the scores:\n")
					print(scores)
				}#in the HNSCC data from 6 paired experiements: 
				#there are 18 circumstances where a patient is found to have multiple scores, one of which above zero
				#there are 7 circumstances where an individual patient has more than one score above zero
			} 
			pgm[sym, p] = max(scores)
		}
	}
	print(dim(pgm))
	print("returning from DrugMultiScoreReduce")
	return(list(txt=pgm, s=s))
}


#getBasicCoverage
#exhumes the base gene symbols from each line of a file
#removes any items in parenthesis and any with "-phosphorylated"
#calls the corsym function to attempt to correct the gene symbols found
#
#takes:		tab: table with column named "gene", containing HUGO gene symbols
extractAndCleanGeneSymbols<-function(tab, study, s, symtab=NULL){
	
	#remove blank rows: 
	blanki = which(tab$gene=="")
	tab2 = tab[!1:nrow(tab)%in%blanki,]
	#remove -domain information
	allSymsPlusPhos = sapply(X=tab2$gene,function(x){strsplit(x=x,split="\\(")[[1]][1]})
	#remove phosphorylation info
	allSyms = sapply(X=allSymsPlusPhos,function(x){strsplit(x=x,split="-phosphorylated")[[1]][1]})
	
	allSyms = unique(allSyms)
	cat("\nIn the input file,",
			length(allSyms),
			"unique gene symbols were found")
	s = setting(s=s, prompt="Have manual gene symbol corrections already been made? (y/n)")
	allSymsCor = corsym(symbol_set=allSyms, symref=getHugoSymbols(paths_detail=study),
											verbose=s$.text=="n")
	allSymsCor = unique(allSymsCor)
	return(list(txt = allSymsCor, s=s))
}#extractAndCleanGeneSymbols

extractCoverageSymbols<-function(s, study, fname=NULL, verbose=T){
	
	#open up file
	if(is.null(fname)){
		fname = "./input/Drug_screen_data_6_paired_experiments.txt"
	}
	dstab = read.table(file=fname,header=T, comment.char="", sep="\t", stringsAsFactors=F)
	#	dstab$TargetName.20.
	#try to find colum
	if(!"gene"%in%colnames(dstab)){#check if column 'gene' exists
		gcindex = NULL
		if("TargetName.20."%in%colnames(dstab)){
			headcol = "TargetName.20."
		}else{
			#else show user the head of the document, prompt then to select the column with gene symbols
			while(T){
				cat("\nThese are the columns of data available:\n")
				print(head(dstab))
				s = setting(s=s, prompt="Please type in the name of the column with the gene symbols: ")
				headcol = s$.text
				if(headcol%in%colnames(dstab)){
					break
				}else{
					cat("\nSorry, the column\"",headcol,"\"could not be found in the table column names.\nPlease try again.\n")
				}
			}
		}
		gcindex = which(colnames(dstab)==headcol)
		colnames(dstab)[gcindex] = "gene"
		#rename the gene symbol column 'gene'
	}
	bc = extractAndCleanGeneSymbols(tab=dstab, study=study, s=s)
	
	return(bc)
}#extractCoverageSymbols


#getPanelCoverage
#get's coverage report for panels that have or have not been run yet. 
#takes 		path_detail: path_detail list to be passed to able4
#					alreadyRun: provide "n" to analyze drug list from a panel not yet run
#											provide "a" to analyze target set from a panel already run (requires the patient gene scores file)
#returns: path summary list from able4
getPanelCoverage<-function(path_detail, study, alreadyRun=NULL, s=list()){
	
	tracker=list()
	
	if(is.null(alreadyRun)){
		s = setting(s=s, "Enter \"g\" to examine coverage using a set of gene names.\nEnter \"d\" to examine coverage using drug names, along with a drug target matrix: ")
		alreadyRun = s$.text
	}
	
	if(alreadyRun == "d"){
		########################################################################
		##################             examine coverage for a panel not yet run, from a drug list
		########################################################################
		s=setting(s=s,prompt="\nPlease enter file name for drug target matrix , or press enter to use the default, ./input/drug_targets_2.txt\n")
		drugTargetsMatrix = s$.text

		s=setting(s=s,prompt="\nPlease enter drug panel file name, or press enter to use the default, ./input/exp006_extract_for_db.txt\n")
		panelFile=s$.text

		tracker[["Drug screen coverage data file used"]] = panelFile
		tracker[["time stamp of drug screen data file"]] = as.character(file.info(panelFile)$mtime)
		#### get targeted genes
		uncorrected = getTargetGenesFromDrugs(drugTargetsFileName=drugTargetsMatrix, panelFileName=panelFile)
		
		#prep the list for coverage analysis:
		
	}else if(alreadyRun=="g"){
		########################################################################
		##################         examine coverage using a file with gene names:
		########################################################################
		#just open the drugscreen patient scores matrix file, extract the target names and clean them. 
		default = "./input/Drug_screen_data_6_paired_experiments.txt"
		
		s = setting(s=s, prompt="Please select a file")
		smSource = s$.text

		tracker[["Drug screen coverage data file used"]] = smSource
		tracker[["Time stamp for drug screen data file"]] = as.character(as.character(file.info(smSource)$mtime))
		tracker[["Time stamp for drug screen data file"]] = as.character(tracker[["Time stamp for drug screen data file"]])
		
	}
	
	res = extractCoverageSymbols(fname=smSource, s=s, study=study)
	targetlist=res$txt
	s=res$s
	
	print("Coverage analysis: ")
	
	panelCoveragePaths = summaryTable(study=study, 
																		coverageGeneDescription="drug_targeted", 
																		coverageDataSetDescription="Drug screen coverage",
																		settings=s, 
																		coverage=targetlist)
	
# 	panelCoveragePaths = summaryTable4(paths_detail=path_detail,
# 																				verbose=T,
# 																				targetname="drug_targeted",
# 																				patientGeneMatrix=targetlistmat,
# 																				dataSetName="Drug screen coverage",
# 																				enrichment_tests=c())

	panelCoveragePaths[["Data_work_up_notes"]] = tracker
	panelCoveragePaths[["settings"]] = s
	#format summaryTable4 output
	return(panelCoveragePaths)
}#getPanelCoverage

#getPatientSubset
#takes: 		patient_gene_matrix: rows = genes, columns = patient ids, cell values = value for gene in patient
#					verbose flag: for interactive mode/accepting user input
#					subset_id: must be provided if verbose is F, this is the string used by grep to search the patient ids
#returns:	patient_gene_matrix for only the subset of patients selected
getPatientSubset<-function(patient_gene_matrix, s, verbose=T, subset_id=NULL){
	if(!verbose&is.null(subset_id)){
		print("If the verbose argument of getPatientSubset(), the subset_id argument must be provided ")
		return(NULL)
	}
	if(verbose) cat("\nThese are the IDs for patient gene score sets that are loaded:\n")
	if(verbose) print(colnames(patient_gene_matrix))
	sub=NULL
	while(T)
	{
		cancer_subtype = subset_id
		
		if(verbose){
			s = setting(s=s,prompt=paste0("Enter a patient id or a substring or regular expression matching the cancer subtype you would like to analyze: ",
																		"\n(ex: \"AML\" for all AML patients with IDs such as \"AMLadult07335\")\n",
																		"(enter * to select all patient ids): "))
			cancer_subtype = s$.text
		} 
		colindex = grep(cancer_subtype, colnames(patient_gene_matrix), ignore.case=T) #get those col.s with the input substring
		sub = as.matrix(patient_gene_matrix[,colindex], nrow=nrow(patient_gene_matrix))#sub is then the extracted patient columns with the matching input string in their col title
		colnames(sub) = colnames(patient_gene_matrix)[colindex]
		rownames(sub) = rownames(patient_gene_matrix)
		if(verbose) print(colnames(sub)) # displays the patient ids so the user can check that they're the correct ones.
		line="y"
		if(verbose){
			s = setting(s=s,prompt="Is this the correct set of patient records? ( y or n ): ")
			line = s$.text
		} 
		if(line == "y"){break}
	}
	return(list(s=s, text=sub))
}

#allow user to select gene score cutoff
#takes: pgm: patient gene matrix, verbose flag, optional cutoff score, must be provided if verbose flag is FALSE
#returns: logic pgm

#'@title Set the cutoff for drug screen score
#'@description Anything above cutoff will be considered to be drug-sensitive.
#'@param pgm A patient gene matrix: bipartate graph with patient ids as columns, genes as rows and cell values the drug screen score. 
#'@param s A settings list object. 
#'@param verbose A logical flag indicating if additional information should be displayed. 
#'@param cutoffScore Optional, the score, above which genes will be considered sensitive. If not provided, interactive prompt will be shown to the user. 
#'@return \code{list} with slots pgm, a logical patient gene matrix and \code{s} the settings list object. 
#'@import grid
#'@import gridExtra
#'@import ggplot2
setCutoff<-function(pgm, s, verbose=T, cutoffScore=NULL){
	#require(grid)
	#require(gridExtra)
	
	if(!verbose&is.null(cutoffScore)){
		print("If the verbose argument of setCutoff() is F, the cutoffScore\n argument must be provided with a numberical value")
		return(NULL)
	}

	print(class(pgm))
	print(mode(pgm))
	print(pgm[1:2,1:2])
	if(verbose) cat("\nThis is the distribution of drug sensitivity scores:")
	
	# 	oldpar <- par(no.readonly=T)
	# #	oldpar <- par()
	# 	while(T){
	# 		tr = try(expr={
	# 	# 			par(mfrow=c(2,1))
	# 	# 			hist(pgm, 
	# 	# 					 ylab="Num. genes throughout cohort",
	# 	# 					 xlab="Gene score", 
	# 	# 					 main="Distribution of drug sensitivity scores for all genes\nexamined in cohort")
			p1 = simpleGGHist(dataSet=as.vector(pgm), showPlot=F,
									 xlab="Gene score", 
									 ylab="Num. genes throughout cohort", 
									 mainTitle="Distribution of drug sensitivity scores for all genes\nexamined in cohort")
			
			tmp = as.vector(pgm)
			tmp = tmp[tmp>0]
			
			p2 = simpleGGHist(dataSet=tmp, showPlot=F,
												xlab="Gene score", 
												ylab="Num. genes throughout cohort", 
												mainTitle="Distribution of drug sensitivity scores in all genes\nwith scores > 0, in all selected patients")
			
	# 			hist(tmp, 
	# 					 ylab="Num. genes throughout cohort",
	# 					 xlab="Gene score", 
	# 					 main="Distribution of drug sensitivity scores in all genes\nwith scores > 0, in all selected patients")
	# 			par(mfrow=c(1,1))
			grid.arrange(p1, p2)

	# 		}, silent=T)
	# 		
	# 		if(!is.error(tr)) break
	# 		print(tr)
	# # 		blnktmp = readline("\nSorry, there was an error displaying the plot.\nPlease increase the size of the 'plots' display window\nthen press enter to try again.")
	# 	}

	if(verbose){
		s = setting(s=s, prompt="Please select a drug screen score cutoff value:")
		cutoffScore = as.numeric(s$.text)
	} 
	pgmout = pgm > cutoffScore
	return(list(text=pgmout, s=s))
}

setDrugScreenCutoff<-function(patient_gene_values, s){
	
	logicPGM=NULL
	while(T){
		res = setCutoff(pgm=patient_gene_values, s=s)
		logicPGM = res$text
		s=res$s
		if(sum(logicPGM)<1){
			cat("\n")
			tmpreadline = readline("There were no drug-sensitive genes found in the cohort.\nTo continue with the analysis, please press enter, then try\na different drug-score cutoff and/or patient cohort selection.")
		}else{
			break
		}
	}

	s$.text = logicPGM
	
	return(s)
}

runPanelAnalysis<-function(path_detail, study, s=list()){
	
	####################################################
	#### load data
	s = setting(s=s, prompt="Please select a file with a drug screen results data set\n")
	drug_scores_fname = s$.text

	drug_scores_fname = checkFileCopyDefault(fname=drug_scores_fname)
	#check if there are only 3 columns --> this would indicate it's a stacked format file
	if(isStackedFormat(drug_scores_fname)){
		print("in stacked format")
		res = DrugScreenFromStackedFormat(fname=drug_scores_fname, s=s)
		patient_gene_levels = res$txt
		s = res$s
		targets = rownames(patient_gene_levels)
	}else{
		print("not in stacked format")
		patient_gene_levels = openPGM(fname=drug_scores_fname)
		approved_hugo=path_detail$symtable
		s = setting(s=s, prompt="Have manual symbol corrections been performed yet for the current data set? (y/n)")
		targets= corsym(symbol_set=row.names(patient_gene_levels),
										symref=path_detail$symtable, 
										verbose=s$.text=="n")
		
		row.names(patient_gene_levels)<-as.character(targets)

		cat("\nChecking that gene symbols in the patient gene scores matrix match approved HUGO symbols. . .\n")
	}
	
	print(dim(patient_gene_levels))

	tracker = list()
	tracker[["Drug screen results data file used"]] = drug_scores_fname
	tracker[["time stamp for drug screen data file"]] = file.info(drug_scores_fname)$mtime
	####################################################
	#####clean out the targets from the target matrix
	
	####################################################
	#output table#
	####################################################
	tout = NULL #this will contain the pathway summary data frame
	####################################################
	
	covset = rownames(patient_gene_levels)

	s = setDrugScreenCutoff(s=s, patient_gene_values=patient_gene_levels)
	logicPGM = s$.text
	print(dim(logicPGM))
	print("Starting path summary")
	stres = summaryTable(study=study, 
							 settings=s, 
							 coverage=targets, 
							 coverageGeneDescription="drug_targeted", 
							 coverageDataSetDescription="Drug screen coverage",
							 pgm=logicPGM, 
							 originalDataMatrix=patient_gene_levels, 
							 activeGeneDescription="drug_sensitive", 
							 dataSetDescription="Drug screen output data")
	
# 	stres$coverage_summary$pathsummary = stres$coverage_summary$pathsummary[,1:5]

	stres$targets = targets

	stres$Data_work_up_notes = tracker

	return(stres)
}


