#drugDbInterface.R
#Finds possible drugs to traget driver or sensitive and driver pathways


#open drug gene association file
#open drug info file

#provide facility for user to pick which drugs they are interested in


# rx="GABA"
rowsContaining<-function(df, rx){
	
	rowsWith = rep(FALSE, times=nrow(df))
	if(rx!=""){
		for(i in 1:ncol(df)){
			rowsWith = rowsWith|grepl(pattern=rx, x=df[,i])
		}
	}
	return(rowsWith)
}

makeStacked<-function(dfin){
	ltmp = list()
	#first make a list out of it
	for(i in 1:nrow(dfin)){
		ltmp[[dfin[i,1]]] = strsplit(x=dfin[i,2], split="; ")[[1]]
	}
	tabd = list_to_table(pth=ltmp)
	
	#now, make the output data frame and fill it
	dfout = data.frame(matrix(ncol=2,nrow=sum(tabd), dimnames=list(1:sum(tabd), c("geneID","drugID"))))
	
	curi = 1
	for(i in 1:nrow(tabd)){
		if(sum(tabd[i,])){
			lindex = curi + sum(tabd[i,]) - 1
			if(lindex == curi){
				dfout[curi,1] = rownames(tabd)[i]
				dfout[curi,2] = colnames(tabd)[tabd[i,]]
			}else{
				dfout[curi:lindex,1] = rownames(tabd)[i]
				dfout[curi:lindex,2] = colnames(tabd)[tabd[i,]]
			}
			curi=lindex + 1
		}
	}
	return(dfout)
}


#loads the drug target data files from drugbank.ca
#returns two column matrix with columns containing 1) hugo symbols <geneID>  and drug bank internal drug ID <drugID>
getDrugTargetData<-function(hugo, 
														fname="./reference_data/drugDB/drugbank/all_target_ids_all.csv"){
	if(!file.exists(fname)){
		warning(paste("The drug target database file,",
									fname,
									"could not be found.\n",
									"Data for this file can be downloaded from http://www.drugbank.ca/",
									"\nFrom the downloads tab/section, download the links for\n",
									"'Drug Target Identifiers', 'All Drugs'\n",
									"'and 'Drug Enzyme Identifiers', 'All Drugs'",
									"\nCurrent working directory:",getwd()))
		
		return(NULL)
	}
	ptm<-proc.time()
	cat("\nLoading drug data from: \n", fname)
	sep=","
	header=T
	quote="\""
	delim=","
	# 	fname="./drugDB/drugbank/all_target_ids_all.csv"

	targtab = read.table(file=fname, sep=sep, header=header, quote=quote, stringsAsFactors=F)
	htab = targtab[targtab$Species == "Homo sapiens",]
	cat("... loaded\n")
	#next, attempt-correction of gene Name column 
	tmpName = corsym(symbol_set=htab$Gene.Name,verbose=F, symref=hugo)
	htab$Gene.Name = tmpName
	
	notHugoIndex = which(!tmpName%in%hugo$Approved.Symbol)
	notHugo= htab[notHugoIndex,]
	
	hextract = hugo[hugo$HGNC.ID%in%notHugo$HGNC.ID,]
	rownames(hextract)<-hextract$HGNC.ID
	
	#now obtain the hugo symbols by their row names
	htab$Gene.Name[notHugoIndex] = hextract[notHugo$HGNC.ID,]$Approved.Symbol
	htab$Gene.Name[is.na(htab$Gene.Name)] = ""
	#extract the GSEA format lines now: 
	targetTable = htab[,c("Gene.Name", "Drug.IDs")]
	
	targetTable = makeStacked(dfin=targetTable)
	totalTime = proc.time()  - ptm
	cat("File load took", totalTime["elapsed"], "seconds.\n")
	return(targetTable)
}#getDrugTargetData

#opens and returns drug meta data file
getDrugData<-function(drugFname = "./reference_data/drugDB/drugbank/drug_links.csv"){
	if(!file.exists(drugFname)){
		warning(paste("The file,",
									drugFname,
									"could not be found.\nData for this file can be downloaded from http://www.drugbank.ca/", 
									"\nBy navigating to DrugBank's 'Downloads' page, 'External Links' tab,",
									"\nthen clicking the 'Download' button for 'Approved drugs'",
									"\n\nCurrent working directory:",
									getwd()))
		
		return(NULL)
	}
	drugTab = read.table(drugFname, sep=",", header=T, stringsAsFactors=F)
	return(drugTab)
}




addPanelMembership<-function(dmd, dtd, STUDY){
	
	cat("Adding drug target counts and numbers of panel members")
	targetTab = as.data.frame(matrix(data=0, nrow=nrow(dmd), 
																	 ncol=2, 
																	 dimnames=list(dmd$DrugBank.ID, c("Total_targets", "New_Targets"))), 
														stringsAsFactors=F)

	if(is.null(STUDY@results$functional_drug_screen_summary)){
		warning("Cannot find summary of drug screen data.\nHas drug screen data been analysed yet?")
		return(dmd)
	}
	
	for(dn in dmd$DrugBank.ID){
		#for current drug, check how many targets there are
		targetTab[dn,] = c(length(dtd$geneID[dtd$drugID==dn]), 
											 sum(!dtd$geneID[dtd$drugID==dn]%in%STUDY@results$functional_drug_screen_summary$targets))
	}
	
	dmd2 = cbind.data.frame(dmd, targetTab, stringsAsFactors=F)
	cat("..\n")
	return(dmd2)	
}


addBangForBuck <- function (pathsToTarget, STUDY, dtd) {
	#' @title addBangForBuck
		 #' @description adds bang for buck analysis which shows how many paths each genes exists in. 
		 #' @param pathsToTarget the set of pathways to be targeted
		 #' @param STUDY the study object
		 #' @param dtd drug target data table as returned by getDrugTargetData()
		 #' @return dtd table with "number of paths containing gene" column added
		 bfbGenes = BangForBuck(darkPaths=pathsToTarget, path_detail=STUDY@studyMetaData@paths)	
		 print("Building output tables...")
		 bfbAddTargets = merge(x=bfbGenes, all.x=T,
		 											y=dtd, 
		 											by.x="row.names", 
		 											by.y="geneID")
		 colnames(bfbAddTargets)<-c("Gene symbol", "Number of dark paths containing gene", "Names of paths containing gene","Drugbank ID")
		 return(bfbAddTargets)
}

test.getPathIdsToTarget<-function(){
	STUDY=getTestStudyObject()
	pids = getPathIdsToTarget(STUDY=STUDY)
	
}

#'@title getPathIdsToTarget
#'@param pathsToSearch the set of pathways to search for targets. If this is provided, function passes this input as the output
#'@param STUDY the study object
#'@return Vector of path ids/path names
getPathIdsToTarget <- function (STUDY, pathsToSearch = NULL) {
	#1) associate paths in overlap analysis with drugs
	# 1a: which paths? aberrational, not drug-targeted
	if(!is.null(pathsToSearch)){
		pathsToTarget=pathsToSearch
		print("pathsToSearch is not null")
	}else if(!is.null(STUDY@results$overlap_analysis)){		
		dark = STUDY@results$overlap_analysis$"Aberration enriched, not drug targeted"
		pathsToTarget = dark$path_id
		print("pathsToSearch is null")
	}else{
		print("pathsToSearch is null and the overlap_analysis is null.. ")
		branchName = selectBranchToFindDrugsFor(study=STUDY)
		pathsToTarget = STUDY@results[[branchName]]$pathsummary$path_id
	}
	
	return(pathsToTarget)
}

test.selectBranchToFindDrugsFor<-function(){
	
	tbranch = selectBranchToFindDrugsFor(study=STUDY)
	
}

selectBranchToFindDrugsFor<-function(study){
	
	snames = names(study@results)
	snames = snames[!snames%in%"overlap_analysis"]
	snames2 = gsub(pattern="_", replacement=" ", x=snames)
	print(matrix(data=snames2, 
							 ncol=1, 
							 dimnames=list(1:length(snames), "data set")))
	while(T){
		line = readline("Please enter the number corresponding to the data set you would like to use to find drug targets: ")
		if(!is.na(as.numeric(line))) break
	}
	
	sel = snames[as.numeric(line)]
	
	return(sel)
	
}


# fdat = read.table(file=fname, 
# 									allowEscapes=F,
# 									stringsAsFactors=F,
# 									quote="", 
# 									fill=F,
# 									strip.white=T,
# 									comment.char="",
# 									header=T, 
# 									sep="|")

readClinicalTrialsFile<-function(fname){
	cat("\nReading fie:", fname,"\n")
	if(!file.exists(fname)){
		warning(paste("The file,",
									fname,
									"could not be found.\nData for this file can be downloaded from https://clinicaltrials.gov/ct2/search/advanced",
									"At this page, select all desired trial phases, click search, click the download button\n",
									"and in the dialog box that comes up, select the maximum number of studies,\n",
									"select all fields,\n",
									"change the format to ",
									"\nCurrent working directory:",getwd()))
		
		return(NULL)
	}
	fdat = try(read.table(file=fname, 
												allowEscapes=F,
												stringsAsFactors=F,
												quote="", 
												fill=T,
												strip.white=T,
												comment.char="",
												header=T, 
												sep="|"), silent=T)
	if(!is.data.frame(fdat)){
		cat("\nNote: there was a issue reading in the data file:\n",fname,"\nHere is the error that was caught:\n\"",fdat,"\"\n")
		cat("Attempting to eliminate unquoted newline characters")
		fdat = readClinicalTrialsFile2(fname)
	}
	print(head(fdat))
	cat("Data dimensions:",  dim(fdat)[1], "rows", dim(fdat)[2], "columns\n")
	return(fdat)
}

readClinicalTrialsFile2<-function(fname="./2012 Pipe delimited text output/subset_clinical_study_noclob.txt"){
	cat("\nReading file:", fname,"\n")
	cat("File will be cleaned of unquoted new-line characters\nusing the assumption that every line starts with the string \"NCT\"\n")
	cat("\nDepending on the file size, this process may be slow...\n")
	f <- readChar(fname, nchars=file.info(fname)["size"], TRUE)
	cat("\nInitial file read in... \ncleaning... ")
	f2 = gsub(pattern="\n(?!NCT)",replacement=" ",x=f, perl=T)
	cat("\nFinal formatting..\n")
	fdat = read.table(text=f2,
										allowEscapes=F,
										stringsAsFactors=F,
										quote="", 
										fill=F,
										strip.white=T,
										comment.char="",
										header=T, 
										sep="|")
	write.table(x=fdat, file=fname,sep="|")
	print(head(fdat))
	cat("Data dimensions:",  dim(fdat)[1], "rows", dim(fdat)[2], "columns\n")
	return(fdat)
}

getClinicalRef_depricated<-function(){
	
	intervFname = "./2012 Pipe delimited text output/interventions.txt"
	interv = readClinicalTrialsFile(fname=intervFname)
	
	cc1 = interv[,c("INTERVENTION_NAME", "NCT_ID")]
	
	intervOtherNameFname = "./2012 Pipe delimited text output/intervention_other_names.txt"
	intervOtherName = readClinicalTrialsFile(fname=intervOtherNameFname)
	
	cc2 = intervOtherName[,c("OTHER_NAME", "NCT_ID")]
	
	intervBrowseFname = "./2012 Pipe delimited text output/intervention_browse.txt"
	intervBrowse = readClinicalTrialsFile(fname=intervBrowseFname)
	
	cc3 = intervBrowse[,c("MESH_TERM", "NCT_ID")]
	
	colnames(cc1)[1]<-"name"
	colnames(cc2)[1]<-"name"
	colnames(cc3)[1]<-"name"
	
	cc1 = rbind(cc1, cc2[!cc2$name%in%cc1$name,])
	
	cc1 = rbind(cc1, cc3[!cc3$name%in%cc1$name,])
	
	return(cc1)
}

loadClinicalTrailsData<-function(fname="./reference_data/drugDB/exported_study_fields_phase3_phase4.tsv"){
	cdat = read.table(file=fname,
										header=T, 
										sep="\t", 
										quote="", 
										comment.char="", 
										stringsAsFactors=F)
	cat("\nFile read in with", ncol(cdat), "columns and", nrow(cdat),"rows.\n")
	return(cdat)
}

appendClincalTrials<-function(dmd, fname="./reference_data/drugDB/exported_study_fields_phase3_phase4.tsv"){
	
	if(!file.exists(fname)){
		warning(paste("The file",fname,"could not be found.\n",
									"This file is needed to append clinical trial information.\n",
									"This file should come from clinicaltrials.gov \n",
									"and should contain columns:\n",
									"Interventions\n Phases\n"))
		return(dmd)
	}
	ptm<-proc.time()
	cat("\nLoading clinical trial data...\n")
	ctdat = loadClinicalTrailsData(fname=fname)
	
	udn = dmd$"Drug name"
	#for each drug name grep all the rows in the clinical trials data that has that drug name
	phases = rep("", times=length(udn))
	cat("Finding clinical trails for drugs...\n")
	print(length(udn))
	pb = txtProgressBar(min=1,max=length(udn), style=3)
	for(i in 1:length(udn)){
		dn = udn[i]
		
		rowsWithDrug = grep(pattern=dn, x=ctdat$Interventions, ignore.case=T)
		drugsTrialPhases = ctdat$Phases[rowsWithDrug] 
		phaseNums0 = gsub(pattern="Phase ", replacement="", x=drugsTrialPhases)
		phaseNums = as.numeric(unlist(strsplit(x=phaseNums0, split="[|]")))
		maxPhase = "No data"
		if(length(phaseNums)){
			maxPhase = max(phaseNums, na.rm=T)
		}
		phases[i] = paste("Phase",maxPhase)
		setTxtProgressBar(pb, i)
	}
	totalTime = proc.time() - ptm
	out = cbind(dmd, phases)
	colnames(out)[ncol(out)] <- "Max.Phase.Clinical.Trial"
	cat("\nFinished appending clinical trial phase (time elapsed:", totalTime[1],"seconds)")
	return(out)
}


toGSEAclinical<-function(dfin){
	
	uids = unique(dfin[,1])
	dfout = cbind.data.frame(uids, rep("", times=length(uids)), stringsAsFactors=F)
	colnames(dfout)<-c("ids", "values")
	rownames(dfout)<-uids
	
	for(i in uids){
		cind = which(dfin$name==i)
		dfout[i,2] = paste(dfin$NCT_ID[cind], collapse=", ")
	}
	return(dfout)
}


importDrugDbData<-function(STUDY){
	#open the data files
	
	dtd1 = getDrugTargetData(hugo=STUDY@studyMetaData@paths$symtable, 
													 fname="./reference_data/drugDB/drugbank/all_target_ids_all.csv")
	dtd2 = getDrugTargetData(fname="./reference_data/drugDB/drugbank/all_enzyme_ids_all.csv",
													 hugo=STUDY@studyMetaData@paths$symtable)
	
	#merge the two gene-drug association lists
	dtd1 = rbind(dtd1,dtd2[!dtd2$geneID%in%dtd1$geneID,])
	
	dtd = dtd1
	
	dtd = unique(dtd)

	return(dtd)
}

test.makeDrugSelectionWorksheet<-function(){
	
	stud = getTestStudyObject()
	res1 = makeDrugSelectionWorksheet(STUDY=stud)
	checkTrue("RXRA"%in%res1[res1$"Drug name"=="Bexarotene",c("Gene symbol")])
	
	limpaths  = stud@results$somatic_mutation_aberration_summary$pathsummary$path_id#= read.table(file="./testLimitPaths.txt", header=F, sep="\t", stringsAsFactors=F)
	res2 = makeDrugSelectionWorksheet(STUDY=stud, pathsToSearch=limpaths)
	
}

#'@title Uses genomic data from the provided Study object to produce a table of pertinent drug-gene associations.
#'@param STUDY A \code{Study} object
#'@param pathsToSearch A \code{vector} of pathway names. The set of cellular pathways that is to be searched for drug targets.
#'@return Spread sheet relating genes in paths to drugs and clinical trials. 
#'@export
makeDrugSelectionWorksheet<-function(STUDY, pathsToSearch=NULL){
	
	dtd0 = importDrugDbData(STUDY=STUDY)
	dmd0 = getDrugData()
	
	#		This many drugs don't have drug meta data
	#check how I'm getting the dmd and see if I can at least get the drug names. .. 
	# 	sum(!dtd0$drugID%in%dmd0$DrugBank.ID)
	# 	[1] 4720
	# 	readline("some rows in output do not have drugs associated")
	
	colnames(dmd0)[2]<-"Drug name"
	rownames(dmd0)<-dmd0$DrugBank.ID
	
	darkPaths = STUDY@results$overlap_analysis$'Aberration enriched, not drug targeted'$path_id
	liableGenes = getGenesFromPaths(pids=darkPaths, STUDY=STUDY)
	
	dtd=dtd0[dtd0$geneID%in%liableGenes,]
	dmd = dmd0[dmd0$DrugBank.ID%in%dtd$drugID,]
	
	ptm<-proc.time()
	#add clinical trial information
	# 	dmd = appendClinicalTrialIDs(dmeta=dmd)
	dmd = appendClincalTrials(dmd)
	totalTime = proc.time()  - ptm
	cat("Appending clinical trial data took", totalTime["elapsed"], "seconds.\n")
	# 	> colnames(dmd)
	# 	[1] "DrugBank.ID"          "Drug name"            "CAS.Number"           "Drug.Type"            "KEGG.Compound.ID"     "KEGG.Drug.ID"        
	# 	[7] "PubChem.Compound.ID"  "PubChem.Substance.ID" "ChEBI.ID"             "PharmGKB.ID"          "HET.ID"               "UniProt.ID"          
	# 	[13] "UniProt.Title"        "GenBank.ID"           "DPD.ID"               "RxList.Link"          "Pdrhealth.Link"       "Wikipedia.Link"      
	# 	[19] "Drugs.com.link"       "NDC.ID"               "ChemSpider.ID"        "BindingDB.ID"         "TTD.ID"               "clinical_trial_IDs" 
	
	dmd = addPanelMembership(dtd=dtd, dmd=dmd, STUDY=STUDY)
	
	#make the data frames
	pathsToTarget = getPathIdsToTarget(STUDY=STUDY, pathsToSearch=pathsToSearch)
	
	# 	genes = getGenesFromPaths(pids=pathsToTarget, STUDY=STUDY)
	
	ptm<-proc.time()
	bfbAddTargets = addBangForBuck(pathsToTarget=pathsToTarget, STUDY=STUDY, dtd=dtd)
	totalTime = proc.time() - ptm
	cat("Appending number of paths containing each gene took", totalTime["elapsed"], "seconds.\n")
	
	# 	> colnames(bfbAddTargets)
	# 	[1] "Gene symbol"                     
	# 	"Number of paths containing gene" 
	# 	"Names of paths containing gene"  
	# 	"Drugbank ID"  
	
	bfbTargDrugData = merge(x=bfbAddTargets, y=dmd, by.x="Drugbank ID", by.y="DrugBank.ID", all.x=T)
	
	bfbTargDrugData = addMutationCount(bfbTargDrugData=bfbTargDrugData, STUDY=STUDY)
	bfbTargDrugData = unique(bfbTargDrugData)
	bfbTargDrugData = addNumberOfPathsTargeted(STUDY=STUDY, dtd=dtd, bfbTargDrugData=bfbTargDrugData)
	bfbTargDrugData = adjustColumns(tab=bfbTargDrugData)
	
	bfbTargDrugData = bfbTargDrugData[!is.na(bfbTargDrugData$"Drug name"),] #remove rows with no drug name. .. . ... .. .. .. 
	
	print("Output tables complete, returning data...")
	return(bfbTargDrugData)
}

adjustColumns<-function(tab){
	print("adjusting columns...")
	corder= c("Gene symbol", 
						"Number of dark paths containing gene", 
						"Aberrations in gene, across cohort", 
						"Names of paths containing gene", 
						"Drug name",
						"Drug.Type",
						"Total_targets",
						"New_Targets",
						"clinical_trial_IDs")
	
	foundc = corder%in%colnames(tab)
	
	if(sum(!foundc)) message("Could not find these data columns: ",paste(corder[!foundc],collpase=" "))
	corder = corder[foundc]
	
	corder= c(corder, setdiff(colnames(tab), corder))
	#setdiff(corder,colnames(tab))
	tab = tab[,corder]
	
	colnames(tab) = gsub(pattern="[_.]", replacement=" ", x=colnames(tab))
	
	return(tab)
}

addNumberOfPathsTargeted<-function(STUDY, bfbTargDrugData, dtd, significanceColumn="hyperg_p_w_FDR"){
	cat("\nAddding number of paths targeted... \n")
	ugene = unique(bfbTargDrugData$"Gene symbol")
	udrug = unique(bfbTargDrugData$"Drugbank ID")
	darkPathCount = rep(0, times=length(udrug))
	names(darkPathCount)<-udrug
	abPathCount = rep(0, times=length(udrug))
	names(abPathCount)<-udrug
	
	#arent the drug selected because they target dark paths?
	
	for(d in udrug){
		if(!is.na(d)){

			geneSet = dtd$geneID[dtd$drugID==d]
			targpaths = whichPaths(STUDY=STUDY, 
														 geneList=geneSet,
														 pathList=STUDY@results$overlap_analysis$'Aberration enriched, not drug targeted'$path_id)
			darkPathCount[d]=length(targpaths)
			
			abPaths = STUDY@results$overlap_analysis$combined_aberrations_summary$pathsummary
			sigAbPaths = abPaths$path_id[abPaths[,significanceColumn]<0.05]
			targpaths2 = whichPaths(STUDY=STUDY, 
														 geneList=geneSet,
														 pathList=sigAbPaths)
			abPathCount[d]=length(targpaths2)
			
		}else{
			cat("NA value found in drug IDs.. skipping\n")
		}
			#if
	}#for
	
	bfbTargDrugData2 = cbind.data.frame(bfbTargDrugData, 
																			darkPathCount[bfbTargDrugData$"Drugbank ID"], 
																			abPathCount[bfbTargDrugData$"Drugbank ID"])
	colnames(bfbTargDrugData2) = c(colnames(bfbTargDrugData),
																 "Number of dark paths targeted by drug",
																 "Number of significantly aberrational paths targeted by drug")

	return(bfbTargDrugData2)
}

#'@title whichPaths
#'@description Indicates the set of pathways a particular set of genes belongs to
#'@param pathList the set of pathways to be selected from
#'@param STUDY the study object
#'@param geneList the list of genes
#'@return A vector of path IDs
whichPaths<-function(pathList, geneList, STUDY){
	
	#first reduce the path matrix to only those targeted
	genesInPaths = intersect(x=geneList, y=colnames(STUDY@studyMetaData@paths$paths))
	reduced = STUDY@studyMetaData@paths$paths[pathList,genesInPaths,drop=F]
	if(is.null(ncol(reduced))){
		print(geneList)
		readline("error, please inform developer (code:142401)")
		return(c())
	}
	#second take the dot product to get paths per gene
	ppg = reduced%*%rep(T, times=ncol(reduced))
	#get the logic matrix
	lppg = ppg>0
	rownames(lppg)<-rownames(ppg)
	#get the path names
  targetedPathNames = rownames(lppg)[lppg[,1]]
	return(targetedPathNames)
}


addMutationCount<-function(bfbTargDrugData, STUDY){
	cat("\nAdding mutation counts... \n")
	geneMuts = as.data.frame(matrix(data=0,
																	ncol=1, 
																	nrow=nrow(bfbTargDrugData)),
													 stringsAsFactors=F)
	colnames(geneMuts) = "Aberrations in gene, across cohort"
	allMutCounts = STUDY@results$somatic_mutation_aberration_summary$genesummary
	
	muti = bfbTargDrugData$"Gene symbol"%in%rownames(allMutCounts)
	
	geneMuts[muti,] = allMutCounts[bfbTargDrugData$"Gene symbol"[muti],]
	
	bfbTargDrugData2 = cbind.data.frame(bfbTargDrugData, geneMuts, stringsAsFactors=F)
	
	return(bfbTargDrugData2)
}



