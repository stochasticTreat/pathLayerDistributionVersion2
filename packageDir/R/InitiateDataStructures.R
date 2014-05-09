#initiate data structures


Path_Detail<-setRefClass(Class="Path_Detail", fields=list(name="character", 
																													file="character", 
																													info="character", 
																													date="character", 
																													source="character",
																													graphite="list",
																													gene_overlap_counts="matrix", 
																													full_path_length="matrix", 
																													symtable="data.frame", 
																													paths="matrix", 
																													symbol_type="character"))

PathSummaryRunner<-setRefClass(Class="PathSummaryRunner", 
															 fields=c('patientGeneMatrix',
															 				 'coverage_summary',
															 				 'coverage',
															 				 'indPatientSummarySettings',
															 				 'original_data_matrix', 
															 				 'min_gene_frequency', 
															 				 'min_gene_count', 
															 				 'individualEnrichment', 
															 				 '.verbose', 
															 				 '.significanceTests', 
															 				 'settings', 
															 				 'targetname', 
															 				 'dataSetName', 
															 				 'study', 
															 				 '.targetMatrix', 
															 				 'patientsum', 
															 				 'path_summary_each_patient', 
															 				 'gene_count_matrix', 
															 				 '.gene_frequency_matrix', 
															 				 '.gene_vector', 
															 				 'genomicnotpw', 
															 				 'active_genes_ea_path'))

setClass("StudyMetaData", representation(paths="Path_Detail",
																				 settings="list",
																				 studyName = "character",
																				 geneIdentifierType="character",	
																				 GeneIdentifierLookup = "data.frame",
																				 RootFile = "character"))

setClass("Study", representation(results="list", 
																 arms = "list",
																 studyMetaData = "StudyMetaData"))

# studyClass<<-setRefClass("Study", fields=c("results","arms","studyMetaData"))

setClass("DataArm", representation(description="character", 
																	 title="character", 
																	 scriptFile="character", 
																	 mainFunction="function"))

setClass("pathSummarySettings", 
				 representation(#original/input data:
					paths_detail = "list",
					patientGeneMatrix = 'matrix',
					original_data_matrix="data.frame",
					target_matrix='matrix',
					#settings/data descript
					targetname="character",
					min_gene_frequency="numeric",
					min_gene_count="numeric",
					dataSetName="character",
					individualEnrichment='logical',
					enrichment_tests="character",
					#miscelanious
					outputFileName="character",
					verbose='logical'))

studyFolder<-function(s){
	return(s@studyMetaData@RootFile)
}

Paths<-function(s){
	return(s@studyMetaData@paths$paths)
}

PathMetaData<-function(s){
	return(s@studyMetaData@paths)
}

studyName<-function(s){
	return(s@studyMetaData@studyName)
}

isSummarySet<-function(testList){
	coreSet=c("genesummary", "patientsums")
	return(sum(!coreSet%in%names(testList))==0)
}

test.isSummarySet<-function(){
	checkTrue(!isSummarySet(testList=list(genesummary="ooo")))
	checkTrue(isSummarySet(testList=list(genesummary="ooo",patientsums="dude")))
}

testGetPathListObject<-function(){
	
	bip = path_detail$paths
	
}


getSummaryObject<-function(pid){
	sl = list()
# 	cat(paste(names(STUDY@results$functional_drug_screen_summary), collapse=", "))
	slnames<-c("pathsummary",
						 "summarystats",
						 "patientGeneMatrix",
						 "patientsums",
						 "genesummary",
						 "active_genes_ea_path",
						 "path_summary_each_patient",
						 "active_genes_not_in_paths",
						 "patientList",
						 "original_data_matrix",
						 "coverage_summary",
						 "targets",
						 "Data_work_up_notes")#,
						 #"settings")
 	for(n in slnames) sl[[n]] = data.frame()
	
	sl$path_summary_each_patient = list()
	sl$patientGeneMatrix = matrix(ncol=length(pid), nrow=0, dimnames=list(NULL, pid))
	sl$genesummary = matrix(ncol=length(pid), nrow=0, dimnames=list(NULL, pid))
	return(sl)
}

dataType<-function(dat){
	if(is.list(dat)){
		if(is.data.frame(dat)) return("data.frame")
		return("list")
	}else if(is.matrix(dat)){
		return(paste("matrix", mode(dat)))
	}else if(is.vector(dat)){
		return(paste("vector", mode(dat)))
	}else if(is.object(dat)){
		if(isS4(dat)) return("S4 object")
		return("classed object")
	}
	return("unknown data type")
}

patients<-function(dset){
		return(colnames(dset$patientGeneMatrix))	
}

targets<-function(dset){
	return(rownames(dset$patientGeneMatrix))
}

patientSums<-function(dset){
	psum = rep(1, times=nrow(dset$patientGeneMatrix))%*%dset$patientGeneMatrix
	return(psum)
}

# 
# for(i in names(STUDY@results$functional_drug_screen_summary)){
# 	cat(i,dataType(STUDY@results$functional_drug_screen_summary[[i]]),"\n",sep="\t\t")
# } 

getDefaultSettings<-function(settingFileName = NULL){
	
	
# 	summaryTableSettings = read.table(file=settingFileName, header=T, sep="\t")
	summaryTableSettings = loadSettings(fname=settingFileName)
	return(list(defaultSummaryTable = summaryTableSettings))
	
}

getStudyObject<-function(	study.name="",
												 	settings=list(), 
													resf=list(), 
													GeneIdentifierLookup=data.frame(), 
													path_detail=Path_Detail$new(), 
													arms=list(), 
													geneIdentifierType="HUGO", 
													rootFolder="./output/"){
	
	if(is.null(resf)) resf = list()
	settings$defaultSummaryTable=getDefaultSettings()[[1]]
	
	if(is.null(geneIdentifierType)) geneIdentifierType = "HUGO"
	
	rootFolder = ifelse(test=is.null(rootFolder), 
				 no=rootFolder,
				 yes="./output/")
	
	rootFolder = ifelse(test=rootFolder=="./output/", 
											no=rootFolder,
											yes=paste(rootFolder, "study_", study.name, "/", sep="", collapse=""))
	
	MetaTmp = new("StudyMetaData", 
								settings=settings,
								studyName=study.name,
								RootFile=rootFolder, 
								paths=path_detail,
								geneIdentifierType=geneIdentifierType,
								GeneIdentifierLookup=GeneIdentifierLookup)
	
	study<-new("Study",
						 results=resf,
						 studyMetaData=MetaTmp)
	
	return(study)
}

##loadDataArm
#takes:	scriptFile:
#				mainFunction: the main function for the data arm, 
#												must implement data arm interface,
#															taking <STUDY> and <settings> as arguments.
#				arms: The list of arms objects
#				title: must be one word, no special characters, used
#				description: This is detailed description of what the data arm will do
#
#example usage:
# arms = loadDataArm(description="Process drug screen data",
# 									 title="functional_drug_screen_summary", 
# 									 scriptFile="./drug_screen_nuevo.R", 
# 									 mainFunction=RunDrugScreen, 
# 									 arms=arms)
loadDataArm<-function(mainFunction, 
											arms, 
											title, 
											description, 
											scriptFile="no script file provided, loading data arm main function from local environment"){
	if(file.exists(scriptFile)) source(scriptFile)
	tmp = new("DataArm", 
						title=title, 
						scriptFile=scriptFile,
						mainFunction=mainFunction, 
						description=description)
	arms[[title]] = tmp
	arms$dictionary[description] = title
	return(arms)
}

isEmpty<-function(val){
	return(length(val)==0)
}

armDescriptionList<-function(study){
	desc = names(study@arms$dictionary)
	return(desc)
}

Paths<-function(s){
	return(s@studyMetaData@paths$paths)
}

FullPathObject<-function(S){
	return(S@studyMetaData@paths)
}

PathMetaData<-function(s){
	return(s@studyMetaData@paths)
}

