
overlapSectionTitles<-function( reslotnames ){
	
	sectiondict = list()
	sectiondict$imageSlots = "Data visualizations"
	sectiondict$"Aberration enriched, not drug targeted" = "Dark pathway summary"
	sectiondict$"Pathway overlaps of genes in aberration enriched, not drug targeted paths" = "Genes overlapping between dark pathways"
	sectiondict$"Aberrationally enriched, containing drug targets" = "Significantly aberrational and functionally targeted"
	sectiondict$"Drug targeted, not aberrationally enriched" = "Drug targeted, not significantly aberrational"
	sectiondict$"Paths containing drug-sensitive genes" = "Paths containing sensitive targets"
	sectiondict$"Enriched for aberration and enriched for sensitive drug targets" = "Significantly aberrational and significantly sensitive"
	sectiondict$"Aberration enriched, containing sensitive targets" = "Significantly aberrational pathways with sensitive targets"
	sectiondict$"settings" = "Settings used for overlap analysis"
	sectiondict$overlap_analysis_each_patient ="Individual patientient overlap analyses"
	sectiondict = sectiondict[ intersect(reslotnames, names(sectiondict)) ]
	
	return(sectiondict)
	
}

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

overlapSectionDescriptionDictionary<-function( reslots ){
	
	sectiondesc = list()
	length(sectiondesc) = length(reslots)
	names(sectiondesc)<-reslots
	sectiondesc$imageSlots = "Visualizations describing the overlap between aberrational and functionally targeted pathways."
	sectiondesc$"Aberration enriched, not drug targeted" = "Pathways which were found to be aberrational but not targeted by the functional assay (ie. by the drug screen or siRNA)."
	sectiondesc$"Pathway overlaps of genes in aberration enriched, not drug targeted paths" = "Genes found in multiple dark pathawys"
	sectiondesc$"Aberrationally enriched, containing drug targets" = "Pathways found to be significantly aberrational and targeted by the functional assay."
	sectiondesc$"Drug targeted, not aberrationally enriched" = "Pathways which were/are targeted by the functional assay but not significantly aberrational."
	sectiondesc$"Paths containing drug-sensitive genes" = "Pathways containing sensitive targets"
	sectiondesc$"Enriched for aberration and enriched for sensitive drug targets" = "Pathways which are significantly aberrational and significantly sensitive (significance determined by enrichment or other test) "
	sectiondesc$"Aberration enriched, containing sensitive targets" = "Pathways which were found to be significantly aberrational and to contain sensitive targets."
	sectiondesc$"settings" = "Settings used for overlap analysis"
	sectiondesc$overlap_analysis_each_patient  = "Overlap analyses for individual patients"
	return(sectiondesc)
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
	sectiondesc$overlap_analysis_each_patient = "The overlap analysis was run on individual patients. Results can be found by following the links below."
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

checkColumnNames<-function(tab, sn){
	if( is.null(colnames(tab)) & ncol(tab)==1 ) colnames(tab) <- sn
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
		
		nshow = max(numRows, sum(tab[,sigcol]<siglim, na.rm=T))
		return(tab[1:min(nshow, nrow(tab)),,drop=F])
	}
	#if nothing else, only return rows 1-20
	return(tab[1:min(numRows, nrow(tab)),,drop=F])
}


addTable<-function( s, tab, sectionDescription, sn, sortAndTrim=TRUE, resSet=NULL ){

	#check to make sure there is actually data in the table
	if( sum(dim(tab)==c(0,0)) ){
		s <- addTo( s, newParagraph("The table is empty") )
		return(s)
	} 
	if( is.vector(tab) & class(tab)=="character" ){
		s <- addTo( s, newParagraph(tab) )
		return(s)
	} 
	mtab = tab
	if(sortAndTrim){
		mtab  = trimTable( tab=tab, siglim=0.05, sigcol=c("hyperg_p_w_FDR"))
	}

	mtab = checkRowNames( tab=mtab )
	mtab = checkColumnNames( tab=mtab, sn=sn )
	
	desc = ifelse(test=is.null( sectionDescription[[sn]] ), 
								yes=sn,
								no=sectionDescription[[sn]] )[[1]] 
	print("desc:")
	print(desc)
	if(nrow(mtab)!=nrow(tab)) desc = paste(desc, ". Showing",nrow(mtab), "out of", nrow(tab),"rows; see 'GET FULL TABLE' for the full set of rows.")
	
	f <- newTable( table=mtab, 
								 desc, 
								 significantDigits=5) # w/ caption
	
	f <- addPathDiagrams( tab=tab, sn=sn, ntab=f, resSet=resSet )
	
	if( nrow(mtab) != nrow(tab) ){
		
		cat(".external table needed.")
		tfname = paste0("./",sn,".txt")
		#check if the table is there already, if not, put it there
		if(!file.exists(tfname)) write.table(x=tab, file=tfname)
		f <- setTableFile( f, tfname )
		
	}

	s <- addTo( s, f )
	
	return(s)
}

pathSlots<-function(){
	
	pslots = c("pathsummary",
						 "Aberration enriched, not drug targeted",
						 "Aberrationally enriched, containing drug targets",
						 "Drug targeted, not aberrationally enriched",
						 "Paths containing drug-sensitive genes",
						 "Enriched for aberration and enriched for sensitive drug targets",
						 "Aberration enriched, containing sensitive targets"
						 )
	return(pslots)
}

addPathDiagrams<-function( tab, sn, ntab, resSet ){
	
	if( !sn%in%pathSlots() | !length(resSet$allPathImages) | !nrow(tab) ){#if the slot does't contain pathways, return the nozzle table element. 
		return(ntab)
	}
	
	pdict = resSet$allPathImages
	
	pidcol = max(1, grep( pattern="path id", x=colnames(ntab$table) ) )
	wantedPaths = ntab$table[,pidcol]
	wantedCleanPNames = gsub(pattern="[/;:*.~<>]", replacement="_", x=wantedPaths) #sanatize pathname
# 	wantedImageFileNames = paste(wantedCleanPNames,"_path_image2.png",sep="") #make the image file name
	
	for ( i in 1:nrow(ntab$table) )
	{
		cat("row number,",i," ")
		#check for path diagram
		pname = wantedCleanPNames[i]
		if( pname%in%names(pdict) ){
			cat("Found diagram for path '",pname,"'\n")
			cur = paste0("./allPathImages/", basename( pdict[[pname]] ) )
			resultsContent <- newFigure( cur )
			result1 <- addTo( newResult( "Pathway diagram" ),
												addTo( newSection( ntab$table[i,pidcol] ),  resultsContent ) );
		}else{
				result1 <- newResult( "" ) ;
		}
		ntab <- addTo( ntab, result1, row=i, column=pidcol )
	}
	return(ntab)
}

formatSettingsToNozzleTable<-function( curel, sn, sectionDescription ){
	
	if( is.list(curel)&!is.data.frame(curel) ) curel = settingsAsDataFrame( settingsData=curel )
	
	curel = checkRowNames( tab=curel )
	curel = checkColumnNames( tab=curel, sn=sn )
	
	f <- newTable( table=curel, 
								 ifelse(test=is.null( sectionDescription ), 
								 			 yes=sn,
								 			 no=sectionDescription)[[1]] 
	) # w/ caption
	return(f)
}

customSlot<-function( s, sn, sectionDescription, curel ){
	
	if(sn == "settings"){

		st <- formatSettingsToNozzleTable( curel=curel, sn=sn, sectionDescription=sectionDescription )

		s <- addTo( s, st )

		return(s)
	}else if(sn == "Data_work_up_notes"){
		
		if( is.list(curel)&!is.data.frame(curel) ) curel = listToDf( lst=curel )
		
		curel = checkRowNames( tab=curel )
		curel = checkColumnNames( tab=curel, sn=sn )
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
	
	if( is.null(dwun) ) return(r)
	if( is.list(dwun) ){
		if( length(dwun)==0 ) return(r)
	}else	if( nrow(dwun)==0 ) return(r)
	
	cat("Number of rows in the data work up notes:",nrow(dwun),"\n")
	
	s <-newSection("Data work up notes")
	table2 <- newTable( dwun,
											"Table of data work up notes" );
	
	for ( i in 1:nrow(dwun) )
	{
		cur = as.character(dwun[i,2])
		
		if ( grepl(pattern=".txt$|.png$", x=cur) )
		{	
			if( grepl(pattern=".png", x=cur) ) resultsContent <- newFigure( basename(cur) )
			if( grepl(pattern=".txt", x=cur) ) resultsContent <- newTable( "Click link to download table", file=basename(cur) )
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

armResultsToNozzle<-function(resSet, resSetName, fname){
	
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
	

	cat("\nOutputting sections to nozzle:\n")
	for(sn in usedSections){
		cat(".. adding",sn, "..")

		curel = resSet[[sn]]
		
		if( sn=="summarystats" ){
			s <- newSection( sectionTitle[[sn]] )
			
			if(sn=="active_genes_not_in_paths"){
				colnames(curel)<-c("Pathway name","Affected genes in pathway")
			}
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn, 
										 sortAndTrim=FALSE, 
										 resSet=resSet )
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
			cat("added settings..\n")
		}else if( class(curel)%in%c("data.frame","matrix") ){
			s <- newSection( sectionTitle[[sn]] )
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn, 
										 resSet=resSet )
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
		cat("*")
	}
	cat("\nAdded all\n")
	#writeReport( r, filename="nozzleTest" );
	return(r)
}

resToReport<-function(resSet, resSetName, fname){
	
	if( grepl( x=resSetName, pattern="overlap" ) ){
		cat("\nOverlap analysis found with name '",resSetName,"'\n")
		nres = overlapAnalysisToNozzle( resSet=resSet, 
																		resSetName=resSetName, 
																		fname=fname )
	}else{
		cat("\nRegular result set found with name '",resSetName,"'\n")
		nres = armResultsToNozzle(resSet=resSet, resSetName=resSetName, fname=fname)
	}

	writeReport( report=nres, filename=fname )
}


addImageSlots<-function( rs, sectionTitle ){
	
	s <- newSection( sectionTitle )
	rsnames = names(rs)
	rs = as.vector(rs, mode="character")
	names(rs)<-rsnames
	
	ovi = grep(pattern="overlap_venn_diagram", names(rs))
	abtarg = grep(pattern="aberrational_and_targeted", names(rs))
	oviFname = paste0("./imageSlots/", basename( rs[ovi] ) )
	abtargFname = paste0( "./imageSlots/", basename( rs[abtarg] ) )
	
	fig1 <- newFigure(file=oviFname, " The overlap between functionally targeted and aberrational pathways. ")
	fig2 <- newFigure(file=abtargFname, " Relation between number of functional targets and degree of aberration in pathways with aberrational genes and targeted genes. ")
	s <- addTo( s, fig1 )
	s <- addTo( s, fig2 )
	
	remainingImages <- setdiff(1:length(rs), c(ovi, abtarg))
	
	if(length(remainingImages)){
		for(im in remainingImages){
			curFname = paste0("./imageSlots/",basename(rs[[im]]))
			capt = gsub(pattern="^graphic",replacement="", x=basename(rs[[im]]))
			capt = gsub(pattern=".png$",replacement="", x=capt)
			capt = gsub(pattern="_",replacement=" ", x=capt)
			s <- addTo( s, newFigure( file=curFname ), capt )
		}
	}
	
	return(s)
	
}

overlapSettings<-function( curel, s ){
	
	for(cs in names(curel)){
		stitle = gsub(pattern="_", replacement=" ", x=cs)
		ss<-newSubSection( stitle )
		st <- formatSettingsToNozzleTable( curel=curel[[cs]], sn=stitle, sectionDescription=NULL )
		ss <- addTo( ss, st )
		s <- addTo( s, ss )
	}
	
	return(s)
}

overlapAnalysisToNozzle<-function( resSet, resSetName, fname ){
	
	if(!require( "Nozzle.R1" )){
		install.packages("Nozzle.R1")
		library("Nozzle.R1")
	}
	# 	[1] "list"
	# 	[1] "data.frame"
	# 	[1] "data.frame"
	# 	[1] "data.frame"
	# 	[1] "data.frame"
	# 	[1] "data.frame"
	# 	[1] "matrix"
	# 	[1] "matrix"
	# 	[1] "list"
	# 	[1] "list"
	# 	[1] "list"
	# 	[1] "list"
	# 	[1] "list"
	sorder = c(
# 		'imageSlots',
						 'Aberration enriched, not drug targeted',
						 'Pathway overlaps of genes in aberration enriched, not drug targeted paths',
						 'Aberrationally enriched, containing drug targets',
						 'Drug targeted, not aberrationally enriched',
						 'Paths containing drug-sensitive genes',
						 'Enriched for aberration and enriched for sensitive drug targets',						 
						 'Aberration enriched, containing sensitive targets',
#						 'combined_aberrations_summary',
#						 'functional_drug_screen_summary',
						 'overlap_analysis_each_patient',
						 'settings'
# 						 'allPathImages'
						)
	#extract out the sections 
	dfslots = c('Aberration enriched, not drug targeted',
							'Pathway overlaps of genes in aberration enriched, not drug targeted paths',
							'Aberrationally enriched, containing drug targets',
							'Drug targeted, not aberrationally enriched',
							'Paths containing drug-sensitive genes',
							'Enriched for aberration and enriched for sensitive drug targets',						 
							'Aberration enriched, containing sensitive targets')

	dfileSections = setdiff(names(resSet), sorder)
	cat("Not including these sections:\n",dfileSections, "\n")
	missingSections = setdiff(sorder, names(resSet))
	usedSections  = intersect(sorder, names(resSet))
	
	sectionTitle  = overlapSectionTitles(reslotnames=usedSections)
	sectionDescription = overlapSectionDescriptionDictionary(reslots=usedSections)
	#build the report using the reportSections
	r <- newCustomReport( gsub(x=resSetName, pattern="_", replacement=" ", fixed=TRUE) )
	
	cat("\nOutputting sections to nozzle:\n")
	for(sn in usedSections){
		cat("\nSection:",sn, "..\n")
		
		curel = resSet[[sn]]
		
		if( sn=="imageSlots" ){
			if( length(resSet[[sn]]) ){
				s <- addImageSlots( rs = resSet[[sn]], sectionTitle = sectionTitle[[sn]] )
				r <- addTo( r, s )
			}
		}else if( sn%in%dfslots ){
			
			s <- newSection( sectionTitle[[sn]] )
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn, 
										 resSet=resSet )
			r <-addTo( r, s )
			
		}else if( sn == "settings" ){
			
			s <- newSection( sectionTitle[[sn]] )
			s <- overlapSettings( curel, s )
			r <- addTo( r, s )
			
		}else if( class(curel)%in%c("data.frame","matrix") ){
			s <- newSection( sectionTitle[[sn]] )
			s <- addTable( s=s, 
										 tab=curel, 
										 sectionDescription=sectionDescription, 
										 sn=sn, 
										 resSet=resSet )
			r <-addTo( r, s )
		}else if( sn=="overlap_analysis_each_patient" ){
			if( length(resSet[[sn]]) ){
				s <- newSection( sectionTitle[[sn]] )
				print("starting patient overlap summaries")
				s <- addAllPatientOverlaps( sec=s, 
																psums=resSet[[sn]], 
																resSetName=resSetName, 
																rootReportName=fname )
				print("adding patient sums to section")
				r <- addTo( r, s )
			}
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
		}#for
	}
	cat("\nAll sections written...\n")
	#writeReport( r, filename="nozzleTest" );
	return(r)
}

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
#'@param resSetName The name of the results set being added to the nozzle page.
#'@export
#'@return The path to the file that the patient summary was saved to.
patientSummaryToNozzleReport<-function(patSum, rootReportName, resSetName){
	
	#take a patient summary
	#set the correct target directory
	#get the name of the patient
	patName = names(patSum)
	#get the name of the root folder
	# 	rootFolder = paste0(gsub(pattern=".nozzleReport[,/.\a-zA-Z0-9]*$", replacement="",x=rootReportName, perl=TRUE), "/")
	rootFolder = dirname(rootReportName) #set the folder
	# 	print(rootFolder)
	#build the links
	linkText = paste0("./path_summary_each_patient/",patName, "/",patName,".nozzleReport.html")
	patFolder = paste0(rootFolder,"/path_summary_each_patient/",patName)
	dir.create(path=patFolder, recursive=TRUE, showWarnings=FALSE)
	fname = paste0(patFolder,"/",patName,".nozzleReport")
	#adjust the title
	resSetTitle = paste(gsub(pattern="_",replacement=" ", x=resSetName),"for patient",patName)
	#send the summary to the regular resToReport
	resToReport(resSet=patSum[[patName]], 
							resSetName=resSetTitle, 
							fname=fname)
	#return the file name that the nozzle was saved to

	return(linkText)
}

patientOverlapToNozzleReport<-function(patSum, rootReportName, resSetName){
	
	#take a patient summary
	#set the correct target directory
	#get the name of the patient
	patName = names(patSum)
	#get the name of the root folder
	# 	rootFolder = paste0(gsub(pattern=".nozzleReport[,/.\a-zA-Z0-9]*$", replacement="",x=rootReportName, perl=TRUE), "/")
	rootFolder = dirname(rootReportName) #set the folder
	# 	print(rootFolder)
	#build the links
	linkText = paste0("./overlap_analysis_each_patient/",patName, "/overlap_analysis/",patName,".nozzleReport.html")
	patFolder = paste0(rootFolder,"/overlap_analysis_each_patient/",patName,"/overlap_analysis")
	dir.create(path=patFolder, recursive=TRUE, showWarnings=FALSE)
	fname = paste0(patFolder,"/",patName,".nozzleReport")
	#adjust the title
	resSetTitle = paste(gsub(pattern="_",replacement=" ", x=resSetName),"for patient",patName)
	#send the summary to the regular resToReport
	resToReport( resSet=patSum[[patName]]$overlap_analysis, 
							resSetName=resSetTitle, 
							fname=fname )
	#return the file name that the nozzle was saved to
	
	return(linkText)
	
}

addAllPatientOverlaps<-function( sec, psums, rootReportName, resSetName ){
	
	patnames = names(psums)
	
	for(pat in patnames){
		cat("\nPatient:")
		print(pat)
		#make the patient summary
		psumFname = patientOverlapToNozzleReport( patSum=psums[pat], 
																							rootReportName=rootReportName, 
																							resSetName=resSetName )
		#add a new paragraph to the report with a link to the patient summary
		
		p <- newParagraph( paste(" Patient",pat,"summary"), asLink(url=psumFname, " Click this link to see patient overlap summary. "))
		
		sec <- addTo( sec, p )
		
	}
	
	return(sec)
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

		p <- newParagraph( paste(" Patient",pat,"summary"), asLink(url=psumFname, " Click this link to see patient summary. "))

		sec <- addTo( sec, p )
		
	}
	
	return(sec)
}

makeSelectNozzleReport<-function(study){
	
	sections = names(study@results)
	while(T){
		print(matrix(data=sections, ncol=1, dimnames=list(1:length(sections), "Results ready for output")))
		uin = readline(prompt="Which data set(s) would you like to build (a) nozzle report(s) for?\n(please enter the numbers sepparated by a space)\n")
		uin  = as.numeric(strsplit(x=uin, split=" ")[[1]])
		if( !sum(!uin%in%1:length(sections)) ) break
		message("There was an error, please try your input again.")
	}
	
	sections = sections[uin]
	
	for(sect in sections){
		
		fileName = paste0(study@studyMetaData@RootFile,"/results/",sect,"/nozzleSummary.nozzleReport")
		resToReport( resSet=study@results[[sect]], 
								resSetName=sect, 
								fname=fileName )
		
	}
	
	
	
}
