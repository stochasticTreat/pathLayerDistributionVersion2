#generic aberration analysis

PrintOpeningMessage<-function(){
	cat("\nThis is the generic aberration interface\n")
	cat("Through this, you can provide a user-defined\n",
			"patient-gene-matrix of any arbitrary data type\n")
	cat("The patient-gene matrix should consist of a\n",
			"matrix, where rows correspond to gene names\n",
			"and columns correspond to patients or samples\n")
	cat("Rows and columns should be named with gene-\n",
			"identifiers and patient identifiers, respectively\n.")
	cat("Values should true/false or 1/0 indicating genes\n",
			"That are aberrational or normal.\n")
}


SelectAndLoadPGM<-function(s,fname=NULL){
	if(is.null(fname)){
		s = setting(s=s,prompt="Please select the patient gene matrix file")
		fname = s$.text
	}
	pgm = openPGM(fname=fname)	
	s$.text = pgm
	return(s)
}


PGMFromVectorWithSelection<-function(s){
	s = setting(s=s,prompt="Please select the file containing your gene list.\nThis file should contain a single column of gene names.\n")
	genefname = s$.text
	genevector = read.table(file=genefname)
	pgm = PGMFromVector(genevector)
	s$.text = pgm
}

getPGM<-function(s){
	while(T){
		s = setting(s=s, prompt="Enter l to do an anlysis using a list of genes.\nEnter p to do an analysis using a patient-gene-matrix")
		uin=s$.text
		if(uin=="l"){
			cat("Vector..\n")
			s = PGMFromVectorWithSelection(s=s)
			return(s)
		}else if(uin=="p"){
			cat("PGM...\n")
			s = SelectAndLoadPGM(s=s)
			return(s)
		}
		cat("\nSorry, that input was not understood, please try again\n")
	}
}




RunGenericEnrichment<-function(settings, study){
	
	s = settings
	
	s = setting(s=s,prompt="Please enter the data type from which your input gene set comes: ")
	aberration_data_type = s$.text
	
	s = setting(s=s,prompt="Please enter the name of the data set that is being entered: ")
	dataSetName  = s$.text
	
	s = setting(s=s,prompt="Please provide a one-word description of the state of the genes that are considered active in the input data set (ex: mutated): ")
	targetName  = s$.text
	
	s = setting(s=s,prompt="Did the gene analysis only examine a subset of genes out of the genome? (y/n)")
	subset  = s$.text
	
	while(T){
		s = setting(s=s, prompt="Is the data from functional analysis (ex: drug screen or RNAi)\nor from aberration analysis (ex: somatic mutation or copy nubmer changes)?\n(please enter f or a)")
		funab = s$.text
		if(funab%in%c("f","a")) break
		message("Sorry, that input was not understood, please try again.")
	}

	target_list=NULL
	
	if(subset=="y"){
		s=setting(s=s,prompt="Please select the file of gene symbols that the analysis was limited to.")
		target_list = read.table(file=s$.text, sep="\t", header=F, stringsAsFactors=F, quote="")
		# target_matrix = getTargetMatrix(tgenes=target_matrix_list[,1], paths=study@studyMetaData@paths$paths)
	}
	
	# 	s = setting(s=s, prompt="Please select the file containing the list of genes you would like to analyze.")
	# 	genevector_input = read.table(file=s$.text, sep="\t", header=F, stringsAsFactors=F, quote="")
	# 	genevector=genevector_input
	# 	pgm = NULL
	
	aberration_data_type = gsub(pattern=" ", replacement="_", 
															x=paste(aberration_data_type, ifelse(test=(funab=="f"),
																							yes="functional",
																							no="aberration"), 
																			  "summary"))
	
	print("Running main generic aberraion analysis function:")
	res = RunGenericEnrichment_main(aberration_data_type=aberration_data_type,
																	study=study,
																	dataSetName=dataSetName,
																	pgm=NULL,
																	genevector=NULL,
																	targetName=targetName,
																	target_list=target_list,
																	path_detail=study@studyMetaData@paths, 
																	s=s)

	return(res)
}

numericThreshold<-function(pgm,targetName,s){
	
	hist(pgm, 
			 main=paste("Distribution of",targetName,"genes"), 
			 xlab=)
	
	prompt=paste("Please describe the numeric cutoff you would like to use to signify a gene that is an aberrational state.",
							 "Please include one of the operators >,<,<=,>= and a number.",
							 "(ex: >=2 to select all values greather than or equal to 2) ",
							 sep="\n")
	while(T){
		s = setting(s=s,prompt=prompt)
		sel=s$.text

		pgmOut = switch(gsub(pattern="[.0-9]*","",x=sel),
								 ">=" = pgm>=as.numeric(gsub(pattern="[<>=]",replacement="",x=sel)), 
								 "<" = pgm<as.numeric(gsub(pattern="[<>=]",replacement="",x=sel)),
								 ">" = pgm>as.numeric(gsub(pattern="[<>=]",replacement="",x=sel)),
								 "<=" = pgm<=as.numeric(gsub(pattern="[<>=]",replacement="",x=sel)),
								 "=" = pgm==as.numeric(gsub(pattern="[<>=]",replacement="",x=sel)),
								 badInput(pgm, sel))
		if(!is.null(pgmOut)) break
	}

	s$.text = pgmOut
	return(s)
}

badInput<-function(pgm,sel){
	cat("Input \'",gsub(pattern="[.0-9]*","",x=sel),"\' not understood\n")
	return(NULL)
}

textThreshold<-function(pgm,targetName, s){
	sumdat = summarize_by(col=matrix(data=pgm,ncol=1),
							 display=T,
							 barPlotTitle=paste("Counts of different types of", targetName, "genes."))
	s = settingList(s=s,
									prompt=paste("Please select the row numbers corresponding to the types of",targetName,"genes you would like to consider to be aberrational."),
									set=sumdat)
	sel=s$.text
	pgm = pgm%in%sel
	return(pgm)
}

applyThreshold<-function(pgm,targetName,s){
	cat("\nApplying threshold...\n")
	s$.text=pgm
	s = switch(typeof(pgm), 
				 "integer"=numericThreshold(pgm=pgm, targetName=targetName,s=s), 
				 "character"=textThreshold(pgm=pgm, targetName=targetName,s=s), 
				 "logical"=s)
	
# 	s$.text=pgm
	
	return(s)
}

RunGenericEnrichment_main<-function(path_detail,
																		study,
																		pgm=NULL, 
																		genevector=NULL,
																		dataSetName="generic aberrational gene data", 
																		aberration_data_type ="generic_aberration_data",
																		targetName="aberrational", 
																		target_list=NULL, 
																		s=s){
	
	if(is.null(pgm)&is.null(genevector)){
		cat("\nLoading patient gene matrix.. \n")
		s = getPGM(s=s)
		pgm=s$.text
	}else if(!is.null(genevector)){
		cat("\nLoading vector of genes.. \n")
		pgm = PGMFromVector(genevector)
	}
	
	s = setting(s=s, prompt="Have manual gene symbol corrections already been conducted? (y/n)")
	customCorrections = s$.text=="n"
	if(s$.text == "n") s[["Have manual gene symbol corrections already been conducted? (y/n)"]] = NULL
	print(".")
	if(customCorrections){
		print(".")
		s = setting(s=s, prompt="Would you like to run manual symbol corrections?\n(If you answer no, automatic corrections will be run) (y/n)")
		customCorrections = s$.text=="y"
	}

	rownames(pgm)<-corsym(symbol_set=rownames(pgm), 
												symref=path_detail$symtable, 
												verbose=customCorrections)
	if(!is.null(target_list)){
		if(!is.atomic(target_list)) target_list = target_list[,1]
		
		target_list <-corsym(symbol_set=target_list, 
												 symref=path_detail$symtable, 
												 verbose=customCorrections)
	}
	originalPGM = pgm
	s = applyThreshold(s=s, pgm=pgm, targetName=targetName)
	pgm = s$.text

	gsum = summaryTable(study=study, 
											coverage=target_list,
											settings=s, 
											originalDataMatrix=originalPGM,
											pgm=pgm,
											dataSetDescription=dataSetName,
											activeGeneDescription=targetName)

	
	typeName=aberration_data_type
	
	glist = list(summary=gsum, resTypeName=typeName, settings=s)
	
	return(glist)
}