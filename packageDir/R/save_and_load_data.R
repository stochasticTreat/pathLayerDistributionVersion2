
#moveImage
#images are placed in a temp folder when they are created, 
#because of this, the need to be moved when the study is saved
checkAndMoveImage<-function(curfile, cfolder){
	# 	curfile = STUDY@results$somatic_mutation_aberration_summary$Data_work_up_notes$"Distributions of mutations before and after fitering"
	#check if it's an image file
	if(is.vector(curfile)){
		if(mode(curfile)=="character"&length(curfile)==1){
			if(grepl(pattern=".png$", x=curfile, ignore.case=T)){
				#check if it's in the tmp folder
				if(grepl(pattern="imageTemp", x=curfile, ignore.case=T)){
					#then move the image
					#re-create the file name
					newFname = paste(cfolder, basename(curfile), sep='/')
					cat("\nCopying to file",newFname,"\n")
					dir.create(path=dirname(newFname), showWarnings=F, recursive=T)
					my.file.rename(from=curfile, to=newFname)
					curfile = newFname
				}
			}
		}
	}
	return(curfile)
}


my.file.rename <- function(from, to) {
	todir <- dirname(to)
	if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
	file.rename(from = from,  to = to)
}

#save_and_load_data.R
nestedListToStudy<-function(res){
	print("inside nestedListToStudy")
	#first check that there was a study meta data folder
	
	if("studyMetaData"%in%names(res)){
		print("opening study saved in new format")
		if("settings"%in%names(res$studyMetaData)){
			#resurect the settings tree
			sl = settingsDFtreeTolistTree(dftree=res$studyMetaData$settings)
		}else{
			sl=list()
		}
		print("Getting study object....")
		#get the study object
		so = getStudyObject(study.name=res$studyMetaData$studyName, 
												geneIdentifierType=res$studyMetaData$geneIdentifierType, 
												path_detail=Path_Detail$new(file=res$studyMetaData$path_file_name),
												settings=sl,
												resf=res$results, 
												rootFolder=res$rootFolder)
		print("got study object.")
	}else if(length(res)>1){
		print("opening study saved in old format")
		#it's a legacy study
		#initilize with out settings
		pathFileName = res$path_file_name
		restmp = res
		res$rootFolder=NULL
		res$path_file_name = NULL
		res$study_name = NULL
		so = getStudyObject(study.name=restmp$study_name,
												rootFolder=restmp$rootFolder,
												path_detail=list(file=restmp$path_file_name), 
												resf=res)
	}else{
		print("Constructing new study.. ..  ")
		#it's a new study
		so = getStudyObject(study.name=res$study_name)
	}
	#set paths
	return(so)
}

#saveStudy: 
#saves all pertinent data of a study object so that it can latter be recalled
#takes: study: study object
#				path: the path to the location where the study folder tree should be saved
#							if this argument is not supplied, the default, '.output' will be used
#returns: nothing
#'@title saveStudy
#'@description Function for saving study object to a human-readable format. 
#'@param study A Study object to be saved. 
#'@param path The file path where the study folder should be saved. 
#'@export
saveStudy<-function(study, path="./output/"){
	# 	setClass("StudyMetaData", representation(paths="list",
	# 																					 settings="list",
	# 																					 study.name = "character",
	# 																					 geneIdentifierType="character",	
	# 																					 GeneIdentifierLookup = "data.frame",
	# 																					 RootFile = "character"))
	# 	setClass("Study", representation(results="list", 
	# 																	 arms = "list",
	# 																	 studyMetaData = "StudyMetaData"))
	#initilize folder
	study_name = study@studyMetaData@studyName
	#append the word "study" to the front of the study name
	#if the study name starts with "test", dont append "study"
	if(!grepl(pattern="^test", x=study_name, ignore.case=T)) study_name = paste0("study_", study_name)
	
	if(study_name!=""){
		path=paste(path,study_name,sep="/")
	}
	if(!file.exists(path)){
		dir.create(path, showWarnings=F, recursive=T)
	}
	#save the studyMetaData object
	savestudyMetaData(study=study, path=path)
	#save the results object
	saveSummary(summ=study@results, study_name="results", path=path)
}

savestudyMetaData<-function(study, path){
	# 	setClass("studyMetaData", representation(paths="list",
	#																						 path_file_name = "character",
	# 																					 settings="list",
	# 																					 studyName = "character",
	# 																					 geneIdentifierType="character",
	# 																					 RootFile = "character"))
	md = study@studyMetaData
	
	#save the settings list
	settingsTrunk = saveSettings(set=md@settings)
	
	mdlist = list(studyName=md@studyName,
								geneIdentifierType = md@geneIdentifierType, 
								path_file_name=md@paths$file, 
								settings=settingsTrunk, 
								RootFile=path)
	saveSummary(summ=mdlist, study_name="studyMetaData", path=path)
}

is.singlevalue<-function(val){
	if(is.data.frame(val)|is.matrix(val)){
		print("data frame or matrix")
		return(FALSE)
	}
	
	if(is.vector(val)){
		if(length(val)>1){
			print("length more than one")
			return(FALSE)
		}else if(length(val)==1){
			print("length is one")
			return(TRUE)
		}else{
			print("Warning, zero length vector in data work up notes")
			return(TRUE)
		}
	}
	print("Some other data type")
	return(FALSE)
}

test.saveDataWorkUpNotes<-function(){
	set = STUDY@results$somatic_mutation_aberration_summary$Data_work_up_notes
	ts = saveDataWorkUpNotes(set=set, fpath="./testoutput/saveDataWorkUpNotes")
	
}

#'@title saveDataWorkUpNotes
#'@description Transforms the data work up notes, from a list, as the are produced by the arms, 
#'@description and makes them into a data fram and an automatically named set of files for any 
#'@description tables contained in the work up notes list.
#'@param set the data work up notes
#'@param fpath the file path that any tables should be saved to.
#'@return data.frame with two columns: the notes and the values of the notes. 
#'@export
saveDataWorkUpNotes<-function(set, fpath){
	if(is.data.frame(set)){
		print("data workup notes already formatted.")
		return(set)
	} 
	#initilize the matrix:
	outmat = matrix(data="", nrow=length(set), ncol=2, dimnames=list(NULL, c("Note","Value")))
	
	for(i in 1:length(set)){
		cname = names(set)[i]
		print(names(set)[i])
		set[[i]] = checkAndMoveImage(curfile=set[[i]], cfolder=fpath)
		if(!is.singlevalue(val=set[[i]])){
			#save the files
			#getFileName
			cur_fname = paste0(fpath,"/","dataWorkUpNotesTable_", i, ".txt")
			#save the file to the right place
			dir.create(path=dirname(cur_fname),showWarnings=F,recursive=T)
			write.table(file=cur_fname,
									x=set[[i]],
									sep="\t")
			#build the data frame output
			outmat[i,] = c(cname, cur_fname)
		}else{
			#build the data frame output
			outmat[i,] = c(cname, as.character(set[[i]]))
		}
	}
	return(outmat)
}

#saveSummary
#creates a nested file folder structure from a nested list structure
#tables, data frames and vectors are saved to their own files
saveSummary<-function(summ, study_name="", path = "./output"){
	#establish the directory path
	if(study_name!=""){
		path=paste(path,study_name,sep="/")
	}
	
	if(!file.exists(path)){
		dir.create(path, showWarnings=F, recursive=T)
	}
	
	cat("\nMain file path:", path,"\n")
	if(length(summ)){
		print(names(summ))
		for(i in 1:length(summ)){#for each item in the summary
			curdat = summ[[i]]
			curname = names(summ)[i]
			cat(" ", curname)
			if(curname =="settings"){ 
				#if settings, then reformat the settings object 
				#and pass back the refrormatted object (then a list)
				cat("\nreformatting settings object\n")
				curdat = saveSettings(set=curdat)
			}
			if(curname=="Data_work_up_notes"){
				print(" reformatting data work up notes\n")
				curdat = saveDataWorkUpNotes(set=curdat, fpath=path)
			}
			checkAndMoveImage(curfile=curdat, cfolder=path)
			
			fname=paste0(path,"/",curname,".txt")
			if(is.data.frame(curdat)|is.matrix(curdat)){
				
				cat(" matrix/data frame to:\t\t",curname, ".txt", "\n",sep="")
				write.table(file=fname,
										x=curdat,
										sep="\t")
			}else if(is.vector(curdat)&!is.list(curdat)){
				
				if(length(curdat)>1){
					cat(" vector with row names to:\t\t", curname, ".txt", "\n",sep="")
					write.table(file=fname,
											x=curdat,
											sep="\t")
				}else{
					cat(" Single value vector to:\t\t",curname, ".txt", "\n",sep="")
					write.table(file=fname,
											x=curdat,
											sep="\t")
				}
			}else if(is.list(curdat)){
				#if it's a list, make a directory out of that element's name and call saveSummary, 
				#passing the path as the path and the variable name as the study name
				if(length(curdat)){#if there's anything in the list
					ndir = paste(path,curname,sep="/")
					cat(" list to:\t\t",curname, ".txt", "\n",sep="")
					dir.create(ndir, showWarnings=F, recursive=T)
					saveSummary(summ=curdat,path=ndir)
				}else{
					cat("***There's nothing in the list\n")
				}
			}else if(!is.null(curdat)){
				cat("\n Error!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Do not know how to save the data in:\n",curname,"\n")
				cat(" This would have been saved to path:\n",path,"\n")
			}
		}
	}
	
	cat("\n</saveSummary()>\n")
}


emptySignature<-function(testdf, pid=""){
	if(is.vector(testdf)&!is.data.frame(testdf)) return(F)
	if(max(dim(testdf))>1) return(F)
	if(min(dim(testdf))==0) return(T)
	#signs of an empty object: 
	#the only value is equal to the pid
	print("the pid:")
	print(pid)
	print("the data frame:::")
	print(testdf)
	print(is.data.frame(testdf))
	
	print(is.matrix(testdf))
	# 	print(testdf[1,1] == pid)
	print("the contents at 1,1")
	print(testdf[1,1])
	if(pid == testdf[1,1]) return(T)
	print("testing 0")
	#the only value is equal to the column name
	if(!is.null(colnames(testdf))){
		print("testing 1")
		if(colnames(testdf) == testdf[1,1]) return(T)
		print("testing 2")
		if(autoRowNames(tab=testdf) & (hasAutoColumnNames(tab=testdf) & testdf[1,1]==colnames(testdf))) return(T)
	}
	print("xx")
	return(F)
}

# 
# ttab = read.table(file="./output/study_integration_test_vs_data_for_xiao_ming/results/functional_drug_screen_summary/path_summary_each_patient/21L-8-1/patientList.txt")
# testdf=ttab
# pid = "21L-8-1"
# emptySignature(testdf=testdf, pid=pid)	

autoRowNames<-function(tab){
	#are the column names automatically generated?
	if(is.null(rownames(tab))) return(T)
	rns = rownames(tab)
	autorowname = sum(rns == as.character(1:length(rns))) == length(rns) #check if the row names are auto generated
	return(autorowname)
}

hasAutoColumnNames<-function(tab){
	if(checkIfDataFame(var=tab)) return(autoMatrixColNames(tab=tab))
	if(autoVectorCol(tab=tab)) return(T)
	print("hasAutoColumnNames() is not sure what type of data is being analyzed")
	return(F)
}

autoMatrixColNames<-function(tab){
	cns = colnames(tab)
	autocolname = sum(cns == paste("V", 1:length(cns), sep="")) == length(cns) #check if the col names are auto generated
	return(autocolname)
}

autoVectorCol<-function(tab){
	if(ncol(tab)>1) return(F)
	if(!is.null(colnames(tab))){ #so it doesn't crash if it's a matrix with NULL for column names
		if(colnames(tab)!="x") return(F)
	}
	return(T)
}

checkIfDataFame<-function(var){
	if(sum(var == paste("V", 1:length(var), sep="")) == length(var)) return(T)
	return(T)
}

autoColNamesDataFrame<-function(tab){
	cns = colnames(tab)
	autocolname = sum(cns == paste("X", 1:length(cns), sep="")) == length(cns) 
	return(autocolname)
}

checkIfMatrix<-function(var, varname){
	if(sum(grepl(pattern="matrix",x=varname, ignore.case=T))) return(T)
	if(autoMatrixColNames(var)) return(T)
	return(F)
}

cleanVector<-function(tab, tabname){
	#was it a named vector
	if(autoRowNames(tab=tab)){
		rownames(tab) = NULL
		return(as.vector(tab[,1])) #the elements in the vector were not named
	}
	#the elements in the vector were named
	tmp = as.vector(tab[,1])
	names(tmp)<-rownames(tab)
	if(tabname=="studyName"){
		print(tab)
		print(tmp)
		readline("Named vector, Press enter to continue")
	} 
	return(tmp)
}

cleanMatrix<-function(tab){
	tab = as.matrix(tab)
	if(autoRowNames(tab=tab)) rownames(tab)<-NULL
	if(autoMatrixColNames(tab=tab)) colnames(tab)<-NULL
	if(is.null(colnames(tab)) & is.null(rownames(tab))){
		attributes(tab)$dimnames<-NULL
		# 		readline("pause in cleanMatrix()")
	} 
	return(tab)
}

cleanDataFrame<-function(tab){
	if(autoRowNames(tab=tab)) rownames(tab)<-NULL
	# 	if(is.null(colnames(tab)) & is.null(rownames(tab))) attributes(tab)$dimnames<-NULL
	return(as.data.frame(tab))
}

CheckConvertDataType<-function(tab, tabname){
	print(tabname)
	#0 check that something is actually there
	if(is.null(tab)) return(NULL)
	print(mode(tab))
	#1 was it a vector
	if(autoVectorCol(tab=tab)){
		print("checking vector")
		return(cleanVector(tab=tab, tabname=tabname))
	}
	#was it a matrix?
	#check matrix autosave
	if(checkIfMatrix(var=tab, varname=tabname)){
		print("checking matrix")
		return(cleanMatrix(tab=tab))
	} 
	#was it a data frame?
	#check data frame autosave
	if(checkIfDataFame(var=tab)){
		print("checking data frame")
		return(cleanDataFrame(tab=tab))
	}
	print("!!!!!!!!!!!!! no other options; not sure what it is...")
	return(tab)
}

summaryObjectFromList<-function(sl, pid){
	
	#first get the blank object
	bl = getSummaryObject(pid)
	
	for(n in names(sl)){ #for each slot
		#second, screen the list to make sure there are no blank slots; if there are, NULL them
		print("------------------current n:")
		print(n)
		print(names(sl))
		print(pid)
		#third, input each slot into the blank summary / bl object
		print("************launching emptySignature()*******************")
		if(!emptySignature(testdf=sl[[n]], pid=pid)){ #if it doesn't have an empty signature
			bl[[n]] = sl[[n]]
		}
	} 
	return(bl)
}

#loadSummary
#loads the summary from one or more analyses, creating a nested list structure
#from the nested directory structure. 
loadSummary<-function(study_name="",path="./output"){
	if(study_name!=""){#if only a path is provided, open everything there
		path = paste(path,study_name,sep="/")
	}
	if(!file.exists(path)){
		cat("\n\nFile cannot be found at", path,"\n")
		tmp = readline("Press enter to continue..")
	}
	flist = list.files(path)
	if(!length(flist)) return(NULL)
	all_sum=list()
	
	#make sure the file doesn't just have an html doc in it. 
	if(sum(!grepl(pattern=".html$",flist))<1) return(NULL)
	
	#if it's a summary file structure, intitiate a summary object
	# 	po = "testReSave111713IntegrationTest/results/functional_drug_screen_summary/coverage_summary/patientsums.txt"
	# 	> path
	# 	[1] "./output/testReSave111713IntegrationTest/results/functional_drug_screen_summary/"
	
	for(i in 1:length(flist)){
		
		bname = flist[i]
		
		cat("\nOpening :", bname, "\n")
		cat("path:", path,"\n")
		fname = strsplit(bname,"\\.txt")[[1]][1]
		
		if("png"==strsplit(x=fname,split="\\.")[[1]][length(strsplit(x=fname,split="\\.")[[1]])]){
			all_sum[[fname]][[strsplit(x=fname,split="_image")[[1]][1]]] = fname
		}else if(file.info(paste(path,bname,sep="/"))$isdir){#if it's a directory, make it a list
			
			cat("\n\t\tSub directory:",bname,"\n")
			all_sum[[fname]] = loadSummary(study_name="",path=paste(path,bname,sep="/"))
			if(isSummarySet(testList=all_sum[[fname]])) all_sum[[fname]] = summaryObjectFromList(sl=all_sum[[fname]], pid=fname)
			
		}else if(grepl(pattern=".txt$", x=bname)){#make sure it's not an html file
			print(".txt file")
			all_sum[[fname]] = try({all_sum[[fname]]= read.table(file=paste(path,bname,sep="/"),
																													 stringsAsFactors=F,
																													 na.strings="-")}, 
														 silent=T)
			if(!is.data.frame(all_sum[[fname]])&!is.matrix(all_sum[[fname]])&mode(all_sum[[fname]])=="character"){ #handle exceptions
				if(grepl(pattern="error", x=all_sum[[fname]], ignore.case=T)){
					print("found error, trying other approaches")
					print(fname)
					if(sum(grepl(pattern="line 1 did not have",x=all_sum[[fname]]))>0){ #error handler
						cat("\nLoad exception handled\n")
						all_sum[[fname]]= read.table(file=paste(path,bname,sep="/"),
																				 stringsAsFactors=F,
																				 na.strings="-",
																				 header=T,
																				 sep="\t")
					}else if(grepl(pattern="first five rows are empty", x=all_sum[[fname]], ignore.case=T)){
						cat("\nEmpty file found..skipping\n")
						print(all_sum[[fname]])
						all_sum[[fname]]=NULL
					}
				}
			}
			#adjust data type:
			all_sum[[fname]] = CheckConvertDataType(tab=all_sum[[fname]], tabname=fname)
		}
	}
	return(all_sum)
}

checkSetDirectoryStructure<-function(verbose=T){
	
	neededDirectories = c("input","output","reference_data")
	wd = getwd()
	notFoundDirs = neededDirectories[!neededDirectories%in%dir(wd)]
	if(length(notFoundDirs)){
		cat("\nNOTICE: \nCreating needed directory structure in current working directory with folders input, output and reference_data.\n")
		for(dn in notFoundDirs) dir.create(path=dn, showWarnings=F)
	}else if(verbose){
		cat("\nCorrect directory structure found.\n")
	}
}

#orchestrates establishment of a new study, or loading of an old study. 
#takes: studyFolderName, path_detail, root, study.name
#returns: study object
#'@title initiateStudy
#'@description Orchestrates establishment of a new study, or loading of an old study
#'@param studyFolderName The name of the study folder. 
#'@param path_detail A Path_Detail object.
#'@param root The path to the folder containing a study to be loaded
#'@param study.name The name of the study to be initiated
#'@return The initiated Study object. 
#'@export
#'@import HGNChelper
initiateStudy<-function(studyFolderName=NULL, 
												path_detail=NULL, 
												root="./output", 
												study.name=NULL){

	cat("\n---------------------Initilizing Study---------------------\n")
	res1 = list()
	studyNameLine="s"
	checkSetDirectoryStructure()
	if(!is.null(study.name)) studyFolderName = paste("study_", study.name, sep="")
	if(is.null(studyFolderName)) studyNameLine= readline("To load a saved study, Enter s\nEnter a study name to start a new study with that name\nTo start a new study with the date as the study name, press enter ")
	
	if(studyNameLine=="s"|!is.null(studyFolderName)){
		
		if(is.null(studyFolderName)){
			folderName = selectStudy(root=root)
		}else{
			folderName = studyFolderName
		}
		
		res1 = loadSummary(path=folderName)
		sname = basename(folderName)
		res1$study_name = gsub(pattern="^study_", replacement="", x=sname)#put the study_name in, removing the prefix
		res1$rootFolder = folderName
		
	}else if(studyNameLine==""){#new study
		# 		print("B")
		study_name = paste("Analysis from",as.character(Sys.time()))
		study_name = gsub(pattern=":",replacement=".",x=study_name)
		res1$study_name = study_name
	}else{#new study
		# 		print("C")
		study_name = gsub(pattern=":|/",replacement=".",x=studyNameLine)
		res1$study_name = study_name
	}	
	
	studytmp = nestedListToStudy(res=res1)
	
	cat("Loading pathways.. . \n")
	if(is.null(studytmp@studyMetaData@paths$file) & is.null(path_detail)){
		if(studyNameLine=="s"){
			tmpbasket = readline(paste("Minor inconsistency found, please check a summary in the study folder \n(",
																 study_name
																 ,"/results/<data_type>) for the paths used, then press enter to continue."))
		}
		path_detail = getPaths()
	}else if(is.null(path_detail)){
		print("Cellular pathway repository not assigned.. loading.. ")
		path_detail=getPaths( path_file=as.character(studytmp@studyMetaData@paths$file))
	}
	# 	print("Check3")
	studytmp@studyMetaData@paths = path_detail
	# 	print("Check4")
	# system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	# 	print("Check5")
	return(studytmp)
}#initiateStudy


test.selectStudy<-function(){
	
	selected = selectStudy()
	
}

#allows interactive selection of a study file
#	takes: <optional> root folder of study
#returns: nested list structre (can subsequently be processed into a study object)
#'@import tcltk
selectStudy<-function(root = "./output"){
	#require(tcltk)
	dc = dir(root)
	while(T){
		if( !length(dc) ){
			tmp = readline("Press any key to continue and select the study folder")
			folderName  = tk_choose.dir(caption="Please select the folder containing the study you would like to load.")
		}else{
			
			#pull the directories out
			dc = dc[file.info(paste(root,dc,sep="/"))$isdir]
			#pull those out with the correct heading
			dc = dc[grep(pattern="^study_|^test", x=dc, ignore.case=T)]
			#make the set of display names
			prettyNames = gsub(pattern="^study_", replacement="", x=dc)
			
			cat("\nHere are the studies currently available:\n")
			print(as.data.frame(matrix(prettyNames,
																 ncol=1,
																 dimnames=list(1:length(prettyNames),
																 							"Study Name"))))
			selstud = readline("Enter the number of the study to be loaded: ")
			folderName = dc[as.integer(selstud)]
			cat("\nLoading data tables from folder named:\n",folderName,"\n")
			folderName  = paste0(root,"/",folderName)
		}
		if( !sum(!c("studyMetaData","results")%in%dir(folderName)) ) break
		message("Sorry, a study could not be found in the folder selected",
				"\n(a study should contain results and studyMetaData folders)",
				"\nPlease try another selection.\n")
		dc=character(0)
		
	}

	return(folderName)
}






