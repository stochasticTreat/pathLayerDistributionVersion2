
coverageFilePrompt<-function(defaultfile = "./OHSUseqdat/coords_with_names_4_16_2010.txt"){
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

loadCoverageData<-function(path_detail, fname=NULL, verbose=F){
	#function to load seq.capture  gene coverage data
	fname = coverageFilePrompt(defaultfile=fname)
	tab = read.table(file=fname,header=F,stringsAsFactors=F)
	#check that usym are approved hugo symbols
	usym = corsym(symbol_set=tab[,1], symref=path_detail$symtable, verbose=verbose)
	usym = unique(usym)
	return(list(file=fname, cov=usym))
}#loadCoverageData


#processSomaticData
#takes:					paths_detail: a path list object
#								removedbSNP: logical, T if dbSNP values should be removed
#								fname:			 The name of the file containing the somatic mutation data to be processed
#								p_subset:		 The subset of patients to analyze. 
#								mutation_selection: the types of mutations to be retained after the filtering step; 
#																		must match those in the Variant_Classification column in the somatic mut. data file; TCGA format
#								verbose:		 logical, if true gives extra output
#								study_name: 	the study name included, used as a file name prefix
processSomaticDataWithCoverage<-function(study_name="no_study_name_provided", 
																				 settings=list(),
																				 study,
														 paths_detail,
														 removedbSNP=F,
														 fname = "./input/HNSCC_broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2many_patient_130418.maf", 
														 verbose=T, 
														 p_subset=NULL, 
														 mutation_selection = c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")){
	print("inside processSomaticDataWithCoverage()")
	fname_base = study_name
	tracker = list()
	tracker[["Study name"]] = study_name
	
	########################################################
	############### 		       Open the file and have a look
	
	#################### readline() ####################################
	if(verbose){
		if("s"==readline(paste("\n\nTo select a nucleotide mutations data set to analyze enter s\nPress enter to load the test data set from\n",fname)))
		{
			fname = file.choose()
		}
	}
	cat("\nLoading file:", fname,"\n")
	tracker[["Data file name"]] = fname
	
	tcga_som_raw = read.delim(file=fname, header=T, sep="\t", stringsAsFactors=F, na.strings="-")
	
	cat("\n",nrow(tcga_som_raw),"rows read in from file.")
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
														 paths_detail=paths_detail)
	tcga_som = filtRes$tcga_som
	tracker = filtRes$tracker
	print("filtered duplicates..")
	###### check number of genes that then have official HUGO symbols
	approvedHugoSymbols = paths_detail$symtable$Approved.Symbol[paths_detail$symtable$Status == "Approved"]
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
	
	# 	if(!is.null(p_subset)){
	# 		subset = pid[pid%in%p_subset]
	# 		cat("***The",length(subset),"patients in the ", patient_set_name, "cohort will be examined.***\n")
	# 		cat("\nTo examine all patients, press escape to cancel, then select all patients from main program interface,", 
	# 				" or reboot R or clear workspace and re-run the somatic mutations processing interface.")
	# 		adjusted_rows = tcga_som[["pid"]]%in%subset
	# 		tcga_som = tcga_som[adjusted_rows,]
	# 	}
	
	cat("\nIn this set of patients", as.character(nrow(tcga_som)), "unique somatic mutations were found.\n")
	
	###### allow PolyPhen adjustment of data
	polyres = addPolyPhenResults(mafData=tcga_som, tracker=tracker)
	tcga_som = polyres$mafData
	tracker = polyres$tracker
	
	cat(".")
	preFilteringGeneSummary = summarize_by(col=tcga_som[,"Hugo_Symbol"], display=F)
	tracker[["Mutations per gene before filtering"]] = preFilteringGeneSummary
	
	preFiltCountByPatient = summarize_by(col=tcga_som$pid, display=F)
	pnames = preFiltCountByPatient$types
	preFiltCountByPatient = matrix(preFiltCountByPatient$counts, dimnames=list(pnames))
	
	tracker[["Mutations per patient, before any filtering"]] = as.matrix(preFiltCountByPatient)
	tracker[["Summary of mutations per patient before filtering"]] = paste(summary(preFiltCountByPatient), collapse="", sep="")
	
	###########################################################
	##########     Mutation type filtering
	
	####################### readline() ####################### 
	filtered = filterMutationType(tcga_som, tracker, verbose)
	
	som_select = filtered$som_select
	tracker = filtered$tracker
	###############################################################
	###################   dbSNP filtering
	
	SNPless = remove_dbSNP(som_select, removedbSNP, tracker)
	som_select = SNPless$som_select
	tracker = SNPless$tracker
	
	###### check number of genes that then have official HUGO symbols
	approvedHugoSymbols = paths_detail$symtable$Approved.Symbol[paths_detail$symtable$Status == "Approved"]
	notApprovedHugoVector = som_select[!som_select[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols),"Hugo_Symbol"]
	numNotHugo = sum(!som_select[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols))
	cat("\nAfter filtering, ",numNotHugo," unique mutations do not have official HUGO symbols.\n",sep="")
	tracker[["Number of mutations without approved HUGO symbols (note: this was checked after filtering and symbol correction)"]] = numNotHugo
	
	write.table(x=som_select[!som_select[,"Hugo_Symbol"]%in%toupper(approvedHugoSymbols),"Hugo_Symbol"],file="./testoutput_num_wo_HUGO.txt")
	##############################
	###### Filtering out mutations with HUGO symbols given as "UNKNOWN"
	
	som_select = som_select[som_select[,"Hugo_Symbol"]!="UNKNOWN",]
	tracker[["Number of unique mutations after removing records with gene identifier given as \"UNKNOWN\""]] = nrow(som_select)
	
	tm=NULL
	covmat = NULL
	
	limitCov = readline("Would you like to limit the genomic coverage of this analysis? (y/n)")
	if(limitCov=="y"){

		covSet = loadCoverageData(path_detail=paths_detail, verbose=T)
		tracker[["Limited the genomic coverage of the input data"]] = "TRUE"
		tracker[["File name of the set of genes the coverage was limited to"]] = covSet$file
		tracker[["Number of genes the coverage was limited to"]] = length(covSet$cov)
		tm = getTargetMatrix(tgenes=covSet$cov, paths=paths_detail$paths)
		
		#reduce the current data to only that in the target set: 
		som_select = som_select[som_select$Hugo_Symbol%in%covSet$cov,]
		
		covpgm = PGMFromVector(genevector=covSet$cov)
		
		covmat = summaryTable(study=study)
		
		# 		covmat = summaryTable4(paths_detail=paths_detail,
		# 																		verbose=F,
		# 																		targetname="sequenced",
		# 																		dataSetName="somatic mutation coverage data",
		# 																		patientGeneMatrix=covpgm)
		
	}
	
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
	
	# 	somatic_summary = summaryTable4(paths_detail=paths_detail,
	# 																	target_matrix=tm,
	# 																	verbose=T,
	# 																	targetname="mutated",
	# 																	dataSetName="somatic mutation data",
	# 																	patientGeneMatrix=somatic_pgm, 
	# 																	outputFileName=table_out_fname)
	som_select = som_select[som_select$Hugo_Symbol%in%covSet$cov,]
	
	covpgm = PGMFromVector(genevector=covSet$cov)
	
	# 		covmat = summaryTable4(paths_detail=paths_detail,
	# 																		verbose=F,
	# 																		targetname="sequenced",
	# 																		dataSetName="somatic mutation coverage data",
	# 																		patientGeneMatrix=covpgm)

	somatic_summary = summaryTable(study=study,
																 settings=settings,
																 coverage=covSet$cov,
																 coverageDataSetDescription="somatic mutation coverage data",
																 coverageGeneDescription="sequenced",
																 pgm=somatic_pgm,
																 dataSetDescription="somatic mutation data",
																 activeGeneDescription="mutated")
	
	#returns: list(pathsummary=ordered_table_out, summarystats=sumStatTab, genelist=patientGeneMatrix)
	somatic_summary[["coverage_summary"]] = covmat
	prefilt = tracker[['Mutations per patient, before any filtering']]
	postfilt = somatic_summary$patientsums
	twoHistOnePlot(dataset1=prefilt, dataset2=postfilt, 
								 x_label="Number of mutations", 
								 y_label="Number of patients",
								 legend_titles=c("Before filtering", "After filtering"),
								 main_title="Distributions of mutations in patients before and after filtering")
	
	somatic_summary[["Data_work_up_notes"]] = tracker
	somatic_summary[["unfiltered_data"]] = tcga_som
	#somatic_summary$genesummary = merge(x=preFilteringGeneSummary, y=somatic_summary$genesummary,all=T)
	#colnames(somatic_summary$genesummary)<-c("Before_filtering", "After_filtering")
	somatic_summary[["preFilteringGeneSummary"]] = preFilteringGeneSummary
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	return(somatic_summary)
}#processSomaticData()