
sectionTitles<-function( reslotnames ){
	sectiondict = list()
	sectiondict$settings = "Settings"
	sectiondict$Data_work_up_notes = "Data work up notes"
	sectiondict$summarystats = "Summary Statistics"
	sectiondict$pathsummary = "Pathway analysis results"
	sectiondict$patientsums = "Patient totals"
	sectiondict$genesummary = "Gene totals"
	sectiondict$path_summary_each_patient = "Individual patient summaries"
	sectiondict$active_genes_not_in_paths = "Affected genes not found in pathways"
	sectiondict$patientList = "Patient ID list"
	sectiondict$coverage_summary = "Summary of platform coverage"
	sectiondict$targets = "Gene identifiers found in input dataset"
	sectiondict$active_genes_ea_path = "Affected genes in each affected pathway"
	sectiondict = sectiondict[ intersect(reslotnames, names(sectiondict)) ]
	
	return(sectiondict)
}

sectionDescriptionDictionary<-function( reslots ){
	
	sectiondesc = list()
	length(sectiondesc) = length(reslots)
	names(sectiondesc)<-reslots
	sectiondesc$settings = "The settings used in the analysis of the current data type."
	sectiondesc$Data_work_up_notes = "Intermediate values found while preparing the data for path analysis."
	sectiondesc$summarystats = "Summary statistics for the cellular pathways, analysis-platform and cohort, as well as relevent platform and pathway coverage."
	sectiondesc$pathsummary = "Results from pathway analysis."
	sectiondesc$patientsums = "Number of affected genes for each patient."
	sectiondesc$genesummary = "Number of cohort patients in which each gene is found to be affected"
	sectiondesc$active_genes_ea_path = "Affected genes found in each affected pathway"
	return(sectiondesc)
}


checkRowNames<-function(tab){
	if(!nrow(tab)){
		message("Warning, table appears to have no rows...")
		return(tab)
	} 
	#check if the row names are the characters "1" through "10"
	rowNamesNotRowNumbers = sum(rownames(tab) != as.character(c(1:nrow(tab))))>0 
	if(rowNamesNotRowNumbers){
		#check if the row names are equal to the values in the first column
		firstColumnNotRowNames = sum(rownames(tab) == tab[,1])==0
		if(firstColumnNotRowNames){
			tab = cbind(rownames(tab),tab)
			if(!is.null(colnames(tab))) colnames(tab)[1]<-""
		}
	}
	return(tab)
}

checkColumnNames<-function(tab){
	colnames(tab)<-gsub(pattern="_", replacement=" ", x=colnames(tab), fixed=TRUE)
	return(tab)
}

customSort<-function( tab, sigcol ){
	
	desc = !sigcol%in%colnames(tab)
	#check for the sigcol
	sigColi = grep(pattern=sigcol, x=colnames(tab))
	sortCol = max(sigColi, ncol(tab)) #doing this handles the case that hyperg_p_w_FDR is not present
	#the actual sort
	tab = tab[order(tab[,sortCol], decreasing=desc),,drop=F]
	
	return(tab)
}

trimTable<-function( tab, numRows=50, siglim=0.05, sigcol=c("hyperg_p_w_FDR") ){

	tab = customSort(tab=tab, sigcol=sigcol)

	#if there are less than 20 rows, just show them all
	if(nrow(tab)<60){ #else just give the whole table
		
		return(tab)
		
	}else if(sigcol%in%colnames(tab)){ #if more than 20 rows and sigcol is there, show all sig, min 20
		
		nshow = max(numRows, sum(tab[,sigcol]<siglim))
		return(tab[1:min(nshow, nrow(tab)),,drop=F])
	}
	#if nothing else, only return rows 1-20
	return(tab[1:min(numRows, nrow(tab)),,drop=F])
}


addTable<-function( s, tab, sectionDescription, sn, sortAndTrim=TRUE ){

	#check to make sure there is actually data in the table
	if( sum(dim(tab)==c(0,0)) ){
		s <- addTo( s, newParagraph("The table is empty") )
		return(s)
	} 
	mtab = tab
	if(sortAndTrim){
		mtab  = trimTable( tab=tab, siglim=0.05, sigcol=c("hyperg_p_w_FDR"))
	}

	mtab = checkRowNames(tab=mtab)
	mtab = checkColumnNames(tab=mtab)
	
	desc = ifelse(test=is.null( sectionDescription[[sn]] ), 
								yes=sn,
								no=sectionDescription[[sn]] )[[1]] 
	
	if(nrow(mtab)!=nrow(tab)) desc = paste(desc, ". Showing",nrow(mtab), "out of", nrow(tab),"rows; see 'GET FULL TABLE' for the full set of rows.")
	
	f <- newTable( table=mtab, 
								 desc, 
								 significantDigits=5) # w/ caption
	
	if(nrow(mtab)!=nrow(tab)){
		
		cat(".external table needed.")
		tfname = paste0("./",sn,".txt")
		#check if the table is there already, if not, put it there
		if(!file.exists(tfname)) write.table(x=tab, file=tfname)
		f <- setTableFile( f, tfname )
		
	}

	s <- addTo( s, f )
	
	return(s)
}

customSlot<-function( s, sn, sectionDescription, curel ){
	
	if(sn == "settings"){
		cat("Class of settings:",class(curel))
		if( is.list(curel)&!is.data.frame(curel) ) curel = settingsAsDataFrame( settingsData=curel )
		
		curel = checkRowNames( tab=curel )
		curel = checkColumnNames( tab=curel )
		
		f <- newTable( table=curel, 
									 ifelse(test=is.null( sectionDescription ), 
									 			 yes=sn,
									 			 no=sectionDescription)[[1]] 
									 ) # w/ caption
		print("..")
		s <- addTo( s, f )
		print("...")
		return(s)
	}else if(sn == "Data_work_up_notes"){
		
		if( is.list(curel)&!is.data.frame(curel) ) curel = listToDf( lst=curel )
		
		curel = checkRowNames( tab=curel )
		curel = checkColumnNames( tab=curel )
		f <- newTable( table=curel, 
									 ifelse(test=is.null( sectionDescription ), 
									 			 yes=sn,
									 			 no=sectionDescription)[[1]] 
									 ) # w/ caption
		s <- addTo( s, f )
		
		return(s)
		
	}
	return(s)
	
}


Data_work_up_notes_nozzle<-function( r, dwun ){
	
	cat("Class dwun 1:", class(dwun),"...\n")
	
	if( is.list(dwun)&!is.data.frame(dwun) ){
		print("converting object class")
		dwun = listToDf( lst=dwun, namesFirstColumn=TRUE )
	} 
	cat("Class dwun 2:", class(dwun),"...\n")
	
	if(is.null(dwun)) return(r)
	if(is.list(dwun)){
		if(length(dwun)==0) return(r)
	}else	if(nrow(dwun)==0) return(r)
	
	s <-newSection("Data work up notes")
	table2 <- newTable( dwun,
											"Table of data work up notes" );
	
	for ( i in 1:nrow(dwun) )
	{
		cur = dwun[i,2]
		if ( grepl(pattern=".txt$|.png$", x=cur) )
		{	
			
			if(grepl(pattern=".png", x=cur)) resultsContent <- newFigure( basename(cur) )
			if(grepl(pattern=".txt", x=cur)) resultsContent <- newTable( "Click link to download table", file=basename(cur) )
			result1 <- addTo( newResult( "Click to see" ),
												addTo( newSection( dwun[i,1] ),  resultsContent ) );
		}
		else
		{
			result1 <- newResult( "" ) ;
		}
		table2 <- addTo( table2, result1, row=i, column=2 );
	}
	s <- addTo( s, table2 )
	r <- addTo( r, s )
	return(r)
}

resultsToNozzle<-function(resSet, resSetName, fname){
	
	if(!require( "Nozzle.R1" )){
		install.packages("Nozzle.R1")
		library("Nozzle.R1")
	}
	sorder = c('settings',
						 'Data_work_up_notes',
						 'summarystats',
						 'pathsummary',
						 "active_genes_ea_path",
						 'patientsums', 
						 'genesummary', 
						 'path_summary_each_patient', 
						 'active_genes_not_in_paths', 
						 'patientList' 
						 # 						 'coverage_summary', 
						 # 						 'targets'
	)
	#extract out the sections 
	dfileSections = setdiff(names(resSet), sorder)
	cat("Not including these sections:\n",dfileSections, "\n")
	missingSections = setdiff(sorder, names(resSet))
	usedSections  = intersect(sorder, names(resSet))
	
	sectionTitle  = sectionTitles(reslotnames=usedSections)
	sectionDescription = sectionDescriptionDictionary(reslots=usedSections)
	#build the report using the reportSections
	r <- newCustomReport( gsub(x=resSetName, pattern="_", replacement=" ", fixed=TRUE) )
	
	handledSlots = c("summarystats")
	cat("\nOutputting sections to nozzle:\n")
	for(sn in usedSections){
		cat("..",sn, "..")

		curel = resSet[[sn]]
		
		if( sn%in%handledSlots ){
			s <- newSection( sectionTitle[[sn]] )
			
			if(sn=="active_genes_not_in_paths"){
				colnames(curel)<-c("Pathway name","Affected genes in pathway")
			}
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn, 
										 sortAndTrim=FALSE )
			r <-addTo( r, s )
		}else if( sn=="path_summary_each_patient" ){
			if( length(resSet[[sn]]) ){
				s <- newSection( sectionTitle[[sn]] )
				print("starting patient summaries")
				s <- addAllPatientSums( sec=s, 
																 psums=resSet[[sn]], 
																 resSetName=resSetName, 
																 rootReportName=fname )
				print("adding patient sums to section")
				r <- addTo( r, s )
			}
		}else if( sn=="Data_work_up_notes" ){
			
			r <- Data_work_up_notes_nozzle( r=r, dwun=curel )
			
		}else if( sn == "settings" ){
			s <- newSection( sectionTitle[[sn]] )
			s <- customSlot( s=s, 
											 sn=sn, 
											 sectionDescription="settings", 
											 curel=curel )
			r <- addTo( r, s )
		}else if( class(curel)%in%c("data.frame","matrix") ){
			s <- newSection( sectionTitle[[sn]] )
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn )
			r <-addTo( r, s )
		}else if (class(curel)=="character"){
			s <- newSection( sectionTitle[[sn]] )
			
			if(length(curel)>1){
				cat(curel)
				lst = newList(isNumbered=TRUE)
				for(li in curel) lst <- addTo( lst, newParagraph(li) )
				s <- addTo( s, lst )
			}else{
				p <- newParagraph( curel[1] )
				s <- addTo( s, p )
			}
			r <-addTo( r, s )
			
		}else{
			s <- newSection( sectionTitle[[sn]] )
			p <- newParagraph( paste("Data type not yet supported: ",class(curel)) )
			s <- addTo( s, p )
			r <-addTo( r, s )
		}
		
	}
	#writeReport( r, filename="nozzleTest" );
	return(r)
}

resToReport<-function(resSet, resSetName, fname){
	nres = resultsToNozzle(resSet=resSet, resSetName=resSetName, fname=fname)
	writeReport( nres, filename=fname )
}

# getStudyFolderName<-function(stud){
# 	folderName = stud@studyMetaData@studyName
# 	#check if it's a test
# 	if(!grepl(pattern="^test", x=studName)){
# 		folderName = paste0("study_", folderName)
# 	}
# 	return(folderName)
# }

#adapter for normal results sets to be transformed into a nozzle report
nozzlesToFileStructure<-function(study, armName, nozzleTitle=""){
	
	stfold = study@studyMetaData@RootFile
	subfold = paste0(stfold, "/",armName, ".nozzleReport")
	resSetName = nozzleTitle
	if(is.null(nozzleTitle)) resSetName = paste(gsub(pattern="_", replacement=" ", fixed=TRUE, x=armName))
	
	resToReport(res=STUDY@results[[armName]], 
							resSetName="Functional drug screen data for HNSCC", 
							fname=subfold)
	
}

#allPatSums = study@results$somatic_mutation_aberration_summary$path_summary_each_patient
#patSum = allPatSums[1]
#'@title Make a Nozzle report from an individual patient's results
#'@description Make a Nozzle report from an individual patient's results
#'@param patSum The list element, including item name of the patient's summary (note this must include the patient name, thus what is given by \code{list[1]} would be correct, not what is subscripted like \code{list[[1]]})
#'@param rootReportName The file path to the root nozzle report ex: "./output/study_mutSigAnalysisOnly/results/somatic_mutation_aberration_summary.nozzleReport"
#'@param resSetName 
#'@export
#'@return The path to the file that the patient summary was saved to.
patientSummaryToNozzleReport<-function(patSum, rootReportName, resSetName){
	
	#take a patient summary
	#set the correct target directory
	#get the name of the patient
	patName = names(patSum)
	#get the name of the root folder
	# 	rootFolder = paste0(gsub(pattern=".nozzleReport[,/.\a-zA-Z0-9]*$", replacement="",x=rootReportName, perl=TRUE), "/")
	rootFolder = dirname(rootReportName)
	# 	print(rootFolder)
	linkText = paste0("./path_summary_each_patient/",patName, "/",patName,".nozzleReport.html")
	patFolder = paste0(rootFolder,"/path_summary_each_patient/",patName)
	dir.create(path=patFolder, recursive=TRUE, showWarnings=FALSE)
	fname = paste0(patFolder,"/",patName,".nozzleReport")
	print("fname inside patientsummarytonozzlereport():")
	print(fname)
	#adjust the title
	resSetTitle = paste(gsub(pattern="_",replacement=" ", x=resSetName),"for patient",patName)
	#send the summary to the regular resToReport
	resToReport(resSet=patSum[[patName]], 
							resSetName=resSetTitle, 
							fname=fname)
	#return the file name that the nozzle was saved to

	return(linkText)
}

addAllPatientSums<-function( sec, psums, rootReportName, resSetName ){
	
	patnames = names(psums)
	
	for(pat in patnames){
		print(pat)
		#make the patient summary
		psumFname = patientSummaryToNozzleReport( patSum=psums[pat], 
																						 rootReportName=rootReportName, 
																						 resSetName=resSetName )
		#add a new paragraph to the report with a link to the patient summary

		p <- newParagraph(asLink(url=psumFname), paste("Patient",pat,"summary"))

		sec <- addTo( sec, p )
		
	}
	
	return(sec)
}

makeSelectNozzleReport<-function(study){
	
	
	sections = names(study@results)
	while(T){
		print(matrix(data=sections, ncol=1, dimnames=list(1:length(sections), "Results ready for output")))
		uin = readline(prompt="Which data set(s) would you like to build (a) nozzle report for?\n(please enter the numbers sepparated by a space)\n")
		uin  = as.numeric(strsplit(x=uin, split=" ")[[1]])
		if( !sum(!uin%in%1:length(sections)) ) break
		message("There was an error, please try your input again.")
	}
	
	sections = sections[uin]
	
	for(sect in sections){
		
		fileName = paste0(study@studyMetaData@RootFile,"/results/",sect,"/nozzleSummary.nozzleReport")
		resToReport(resSet=study@results[[sect]], 
								resSetName=sect, 
								fname=fileName)
		
	}
	
	
	
}
