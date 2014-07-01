#initiate data structures

#'@title Path_Detail class
#'@description Contains all information relevent to cellular pathways used in pathway analysis.
#' @section Fields:
#' \describe{
#' \item{\code{name}:}{Object of class \code{"character"} The name of the pathway repository used. }
#' \item{\code{file}:}{Object of class \code{"character"} The file from which the pathway repostory was loaded.}
#' \item{\code{date}:}{Object of class \code{"character"} The date the pathway repository was loaded.}
#' \item{\code{info}:}{Object of class \code{"character"} Any additional information reguarding the pathway repository.}
#' \item{\code{source}:}{Object of class \code{"character"} The source of the pathway repository. For example, Reactome, via the graphite bioconductor package )}
#' \item{\code{graphite}:}{Depricated. Only here for legacy purposes.} 
#' \item{\code{gene_overlap_counts}:}{Object of class \code{"matrix"} One column, with row names given as gene names. Values are the number of pathways each gene belongs to.} 
#' \item{\code{full_path_length}:}{Object of class \code{"matrix"} One column, pathway names given as row names. Values given as the number of genes annotated to each pathway.}
#' \item{\code{symtable}:}{Object of class \code{"data.frame"} The official, approved set of symbols or gene identifiers and any cross references.}
#' \item{\code{paths}:}{Object of class \code{"matrix"} The actual storage of pathways in bipartate graph format, with row names as the path names, column names as gene identifiers and values logical representing gene membership in pathawys.}
#' \item{\code{symbol_type}:}{Object of class \code{"character"} The type of gene identiers used, for example HUGO or Uniprot.}
#' }
#' @import methods 
#' @exportClass Path_Detail
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

#'@title PathSummaryRunner
#'@description This reference class object is used in the pathway analysis. It provides all available data for internal functions and user-defined functions to conduct pathway analysis. 
#' @section Fields:
#' \describe{
#' \item{\code{patientGeneMatrix}:}{
#'	 \code{Object of class "Matrix"} Bipartate graph with patient ids as columns, genes as rows and values, logical indicating if gene is 'active' (drug sensitive, aberrational, etc.) in patient. 
#'	 }
#' \item{\code{coverage_summary}:}{
#'	 \code{Object of class "list"} Same data as is found in [[result]] set, only filled if analysis platform in question is coverage-limited, ie, coverage slot of PathSummaryRunner is filled. 
#'	 }
#' \item{\code{coverage}:}{
#'	 \code{Object of class "Character"} A character vector containing symbols for the set of species (ie genes) covered by the analysis platform. If this is provided, pathway analysis will be limited to included only genes in the coverage set. This should be provided if an analysis platform is not considered to have full-genome coverage (ex: sequence capture data for 3000 unique genes.).
#'	 }
#' \item{\code{indPatientSummarySettings}:}{
#'	 \code{Object of class "list"} Settings object for analysis of individual patients (as opposed to a full cohort). 
#'	 }
#' \item{\code{original_data_matrix}:}{
#'	 \code{Object of class "matrix"} The original, untransformed data set. 
#'	 }
#' \item{\code{min_gene_frequency}:}{
#'	 \code{Object of class "numeric"} The minimum frequency across the cohort for a gene to be considered as active (aberrational, drug sensitive, etc..) in the pathway analysis. (default 0 imposes no threshold)
#'	 }
#' \item{\code{min_gene_count}:}{
#'	 \code{Object of class "numeric"} Analogous to min_gene_frequency -- the minimum number of patients with gene in active state for it to be considered active in the pathway analysis. (default 1 imposes no threshold)
#'	 }
#' \item{\code{individualEnrichment}:}{
#'	 \code{Object of class "logical"} A flag indicating if path analysis should be run for individual patients. (note, if individual patient analysis is run, min_gene_frequency and min_gene_count will automatically be set to 0 and 1, respectively, in the analysis of individual patients (though will be maintained in the analysis of any full cohort))
#'	 }
#' \item{\code{.verbose}:}{
#'	 \code{Object of class "logical"} Flag indicating if extra output should be shown. 
#'	 }
#' \item{\code{.significanceTests}:}{
#'	 \code{Object of class "list"} Each slot contains a function. Each function implements the path test interface, accepting a PathSummaryRunner as an input, and returning a data frame with path names or identifiers as row names and values indicating pathway's significance (ex: p-values). 
#'	 }
#' \item{\code{settings}:}{
#'	  \code{Object of class "list"} The settings list object. 
#'	 }
#' \item{\code{targetname}:}{
#'	  \code{Object of class "character"} One word description of the genes affected (ex: mutated)
#'	 }
#' \item{\code{dataSetName}:}{
#'	  \code{Object of class "character"} Description of the data set analyzed by the current arm. 
#'	 }
#' \item{\code{study}:}{
#'	  \code{Object of class "Study"} The study object for the current analysis. 
#'	 }
#' \item{\code{.targetMatrix}:}{
#'	  \code{Object of class "matrix"} Bipartate graph matrix. This is a version of the pathway bipartate graph reduced to only contain genes considered "active" (aberrational or functionally sensitive). This is the primary matrix used in computing pathway significance. 
#'	 }
#' \item{\code{patientsum}:}{
#'	  \code{Object of class "matrix"} Numeric values indicating the number of genes affected in each patient. Patient IDs are given as row names, values are the number of genes affected in each patient.
#'	 }
#' \item{\code{path_summary_each_patient}:}{
#'	  \code{Object of class "list"} The results sets for each patient. List names are set as patient identifiers. 
#'	 }
#' \item{\code{gene_count_matrix}:}{
#'	  \code{Object of class "matrix"} Values are \code{"numeric"} Row names are gene identifiers. Values are the number of patients with affected gene.
#'	 }
#' \item{\code{.gene_frequency_matrix}:}{
#'	 \code{Object of class "matrix"} The frequency across cohort that each gene is active.
#'	 }
#' \item{\code{.gene_vector}:}{
#'	 \code{Object of class "matrix"} Gene identifiers are given as row names. Values are logical, indicating if genes are considered 'active' in current cohort (or patient, for individual patient analysis).
#'	 }
#' \item{\code{genomicnotpw}:}{
#'	 \code{Object of class "matrix"} Values are \code{"character"} identifiers for affected genes not found in pathways.
#'	 } 
#' \item{\code{active_genes_ea_path}:}{
#'	 \code{Object of class "matrix"} Values are pasted together \code{"character"} vectors of identifiers for genes found active in each pathway. Row names are given as pathway identifiers. 
#'	 }
#' \item{\code{pathsummary}:}{
#'	 \code{Object of class "data.frame"} containing main coverage and enrichment summary data for each affected pathway.
#'	 }
#' }
#'@import methods
#'@exportClass PathSummaryRunner
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
															 				 'active_genes_ea_path', 
															 				 "pathsummary"))

#'@title StudyMetaData
#'@description Object containing crucial study meta data. 
#'@slot paths Path_Detail reference class object
#'@slot settings A list of settings list objects. One slot in list is generally provided for each arm. 
#'@slot studyName String. The name of the study and name of the file folder used to save the study.
#'@slot geneIdentifierType String. The name of the type of gene identifiers used in the study. 
#'@slot GeneIdentifierLookup depricated. A slot for a set of official, approved gene identifiers and mappings for unofficial gene identifiers. 
#'@slot RootFile character data. The path to folder where the study is saved. 
#'@export
setClass("StudyMetaData", representation(paths="Path_Detail",
																				 settings="list",
																				 studyName = "character",
																				 geneIdentifierType="character",	
																				 GeneIdentifierLookup = "data.frame",
																				 RootFile = "character"))

#'@title Study
#'@description Contains all data related to a single study, including results for each data input arm and the overlap analysis, the main functions for each data input arm, and the study meta data object.
#'@slot results list. Names provided as data input arm names. Each slot contains the results set from its respective data input arm.  
#'@slot arms A list with each slot containing a DataArm object, slot names given as the names of each data input arm
#'@slot studyMetaData the StudyMetaData object
#'@export
setClass("Study", representation(results="list", 
																 arms = "list",
																 studyMetaData = "StudyMetaData"))

#'@title DataArm
#'@description Contains the information necessarry to connect a data input arm to the interactive version of this program. 
#'@slot description character. The description of the data input arm provided in the main menu.
#'@slot title character. The title of the path analysis arm. This is the title by which the results from the analysis arm will be saved in the results list.
#'@slot scriptFile A string providing the full or relative file path of a script file containing all functions pertinent to running the data arm. This can be omitted if all functions are already loaded into the namespace, and the mainFunction slot is filled. 
#'@slot mainFunction function. The name of the main function for the arm in question.
#'@export
setClass("DataArm", representation(description="character", 
																	 title="character", 
																	 scriptFile="character", 
																	 mainFunction="function"))


#'@title loadDataArm
#'@description Connects data input arm to program, providing essential information and the main execution function. 
#'@param mainFunction The main execution function for the arm. If the script provided, the function does not need to be loaded into memory when the call is made to loadDataArm -- the script file will be run then the function will be called upon to be placed into the mainFunction slot of a DataArm object.
#'@param arms The list containing all DataArm objects for the study.
#'@param title character string. The one word description of the data arm. (ex: functional_drug_screen_summary) This will be used internally to refer to any results associated with the data arm, and externally to name folders when saving data. 
#'@param description The description of the data input arm to be used as a main menu option. (ex: "Process drug screen data")
#'@param scriptFile A file providing the script where the functions implementing the data import arm can be found. Note: This allows that a script can be provied which contains the definition of the function provided for the mainFunction argument. Thus, this script will be executed before the mainFunction is called or referened. 
#'@return The arms list with the additional arm added. 
#'@export
#'@examples
#'study = getStudyObject(study.name="testDataSets", 
#' 											 geneIdentifierType="HUGO")
#'arms = slot(object=study, name="arms")
#'arms = loadDataArm(description="Load drug screen data",
#'	title="functional_drug_screen_summary", 
#'	mainFunction=packageDir:::RunDrugScreen, 
#'	arms=arms)
loadDataArm<-function(mainFunction, 
											arms, 
											title, 
											description, 
											scriptFile="no script file provided, loading data arm main function from local environment"){
	cat("\nSetting arm named '",title,"'... ")
	if(file.exists(scriptFile)) source(scriptFile)
	tmp = new("DataArm", 
						title=title, 
						scriptFile=scriptFile,
						mainFunction=mainFunction, 
						description=description)
	arms[[title]] = tmp
	arms$dictionary[description] = title
	cat(" .. arm set.\n")
	return(arms)
}


studyFolder<-function(s){
	return(s@studyMetaData@RootFile)
}

Paths<-function(s){
	return(s@studyMetaData@paths$paths)
}

# PathMetaData<-function(s){
# 	return(s@studyMetaData@paths)
# }

studyName<-function(s){
	return(s@studyMetaData@studyName)
}

isSummarySet<-function(testList){
	coreSet=c("genesummary", "patientsums")
	return(sum(!coreSet%in%names(testList))==0)
}

testGetPathListObject<-function(){
	
	bip = path_detail$paths
	
}


getSummaryObject<-function(pid){
	sl = list()
	# cat(paste(names(STUDY@results$functional_drug_screen_summary), collapse=", "))
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

#'@title getDefaultSettings
#'@description Returns the default settings object for the summaryTable path analysis. 
#'@return A list with a slot named "defaultSummaryTable", containing the default settings for the summaryTable function.
#'@export
#'@examples
#'defSettings = getDefaultSettings()
#'print(names(defSettings$defaultSummaryTable))
getDefaultSettings<-function(){
	# summaryTableSettings = read.table(file=settingFileName, header=T, sep="\t")
	summaryTableSettings = loadSettings(fname=NULL)
	return(list(defaultSummaryTable = summaryTableSettings))
	
}

#'@title getStudyObject
#'@description Function to initilize a new study object. 
#'@param study.name The name of the study (this will be used as a folder name if study is saved, so this should be valid as a file name.)
#'@param settings A settings list object. 
#'@param resf a results list object.
#'@param GeneIdentifierLookup Depricated. A table of approved gene identifiers. 
#'@param path_detail A Path_Detail object. 
#'@param arms The list of arm objects. 
#'@param geneIdentifierType String, a one word name of the type of gene identifiers used (ex : HUGO, Uniprot)
#'@param rootFolder character string. The file path where the study should be saved. 
#'@export
#'@examples
#'study = getStudyObject() #this will give a study object with no pathways loaded (pathways must first be loaded for analyses to take place)
#'study2 = getStudyObject(study.name="testStudy", path_detail=getDefaultPaths())
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

isEmpty<-function(val){
	return(length(val)==0)
}

#'@title Get all arm descriptions
#'@description Retreives arm descriptions from the arms slot of a Study object. This thus contains descriptions of all the arms loaded and ready for use. 
#'@param study A \code{Study} object
#'@return \code{vector} of study arm descriptions. 
#'@export
armDescriptionList<-function(study){
	desc = names(study@arms$dictionary)
	return(desc)
}

Paths<-function(s){
	return(s@studyMetaData@paths$paths)
}



