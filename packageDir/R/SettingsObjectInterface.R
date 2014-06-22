#SettingsObjectInterface.R

#setting list standard: 
#slot names are the prompts, slot values are the values selected by the user
#slot names starting with a "." are treated as special values that the user is not prompted for
#setings list will have appended to it these items: 
#		.interactive: 	indicates if user will be prompted even if a value has already been set
#		.text: 					the currenly active value, for immediate use by program
#		.study:					general information about the study such as: 
#														.study$paths: 	the path_detail object
#														.study$results: the results object for the study

# #how to make a dataWorkUpTracker: 
# #must take 
# settings = list()
# settings[["how cute are sloths?"]] = "Very, very cute"
# outputs = list(output_one = c("not too much here", "sloths are cute"))
# s = list(settings=settings, outputs=outputs)
# 
# WorkUpTracker

settingsDFtreeTolistTree<-function(dftree){
	print("..in settingsDFtreeTolistTree()")
	
	for(branch in 1:length(dftree)){
		if(class(dftree[[branch]])=="list"){
			dftree[[branch]] = settingsDFtreeTolistTree(dftree[[branch]])
		}else{
			dftree[[branch]] = dfToList(dftree[[branch]])
		}
	}
	print("leaving settingsDFtreeTolistTree")
	return(dftree)
}

dfToList<-function(df){
	out=list()
	for(n in rownames(df)){
		out[n] = df[n,1]
	}
	return(out)
}

listToDf<-function(lst){
	print("in listToDf")
	vals = sapply(X=names(lst), FUN=function(x){lst[[x]]})
	out = data.frame(vals=as.matrix(vals,ncol=1), stringsAsFactors=F)
	out[,1] = as.character(out[,1])
	return(out)
}


#'@title Load settings from a file folder tree into a Study object. 
#'@description Allows loading of settigns into a \code{Study} object. If a directory path name is provided in the fname argument, settings will be loaded from that file. If no directory is provided, the user will be prompted to interactively select a settings folder. 
#'@param study A \code{Study} object. 
#'@param fname Optional, a \code{string} giving the absolute or relative paths of the settings to be loaded. 
#'@param root Optional, a \code{string} giving a path to be checked for \code{Study} objects from which to retreive settings. 
#'@return A \code{Study} object with settings loaded. 
#'@export
selectAndLoadSettings<-function(study, fname=NULL, root="./output/"){
	fullPath=""
	while(T){
		if(is.null(fname)){
			#find all possible studies
			studs = getStudies(root)
			if(is.null(studs)){
				tmp = readline("\nNo studies found, press any key to continue and select a settings folder\n")
				fullPath = tk_choose.dir()
			}else{
				#prompt user to select study
				print(as.data.frame(matrix(studs,ncol=1,dimnames=list(1:length(studs),"Study Name"))))
				selstud = readline("Enter the number of the study from which settings are to be loaded\n(or enter f to select a folder containing settings): ")
				if(selstud=="f"){
					fullPath = tk_choose.dir()
				}else{
					#find the settings in that study
					study_name = studs[as.integer(selstud)]
					cat("\nChecking for settings record.. \n")
					#load the settings
					#form file path
					fullPath=paste(root, study_name, "/", "studyMetaData/settings/", sep="")
				}
			}
		}else fullPath=fname
		#check if the settings file exists
		if(file.exists(fullPath)) break
		cat("\n**** Sorry, a settings file cannot be found at:\n", fullPath, 
				"\n\nPlease select another study in which to find settings.\n")
		fname = NULL
	}
	
	loadedSettings = loadAllSettings(fname=fullPath)
	
	#place the settings in the study object
	
	study@studyMetaData@settings = loadedSettings
	study@studyMetaData@settings$defaultSummaryTable = getDefaultSettings()[[1]]
	return(study)
}


loadAllSettings<-function(fname){
	stree = loadSummary(path=fname, study_name="")
	ltree = settingsDFtreeTolistTree(dftree=stree)
	return(ltree)
}



loadSettings<-function(fname=NULL){
	
	print("Inside loadSettings()  ... ")
	if(is.null(fname)){
		message("Loading the default summary table settings")
		# data("defaultSummaryTable140504", verbose=T)
		lres = load(file=system.file("extdata/defaultSummaryTable140504.rdata", package = "packageDir"), 
				 verbose=T)
		out = get(lres[1])
	}else{
		settingsData = read.table(file=fname, header=T, sep="\t", comment.char="", stringsAsFactors=FALSE)
		out = dfToList(settingsData)
	}

	return(out)
}

getStudies<-function(root="./output/"){
	dc = dir(root)
	studs = c()
	dc = dc[grep(pattern="^study_|^test", x=dc, ignore.case=T)]
	for(f in dc){
		if(file.info(paste(root,f,sep=""))$isdir){
			studs = c(studs,f)
		}
	}
	return(studs)
}

#allows selection of a study file
#	takes: <optional> root folder of study
#returns: nested list structre (can subsequently be processed into a study object)
selectAndLoadResults<-function(root = "./output/"){
	studs = getStudies(root)
	cat("\nHere are the studies currently available:\n")
	print(as.data.frame(matrix(studs,ncol=1,dimnames=list(1:length(studs),"Study Name"))))
	selstud = readline("Enter the number of the study to be loaded: ")
	study_name = studs[as.integer(selstud)]
	cat("\nLoading data tables.. . \n")
	results = loadSummary(study_name=study_name,path="./output/")
	results$study_name = gsub(pattern="^study_", replacement="", x=study_name)
	return(results)
}

#saveSettings
#takes settings and performs one of several options based on their state: 
#if set is a list of settings lists: the settings lists are transformed to data frames and the 
#set is returned.
#If set is a settings list: it is transformed to a data frame and returned. 
#If set is a list of data frames, it is returned. 
#	takes: 	set: a list, described above
# 					#depricated: fname: a file name; set will be saved here if this is provided
#returns: list of data frames
saveSettings<-function(set, fname=NULL){
	
	if(!length(set)) return(NULL)#if nothing is there
	trunk=NULL
	if(is.data.frame(set)|is.data.frame(set[[1]])){
		trunk = set
	}else if(is.vector(set[[1]])&!is.list(set[[1]])){
		trunk = settingsAsDataFrame(set)
	}else if(is.list(set[[1]])){
		trunk = set
		for(branch in 1:length(trunk)){
			print(names(trunk)[branch])
			# 			trunk[[branch]] = settingsAsDataFrame(settingsData=trunk[[branch]])
			trunk[[branch]] = saveSettings(set=trunk[[branch]])
		}
	}else{
		print("!!!!!! Error: unknown data submitted as setting list")
	} 
	#optional save
	if(!is.null(fname)){#if a file name is provided
		saveSummary(summ=trunk, study_name="settings", path=fname)
	}else{
		print("Returning settings list tree with data frame leaves.")
		print(typeof(trunk))
		print(is.data.frame(trunk))
		return(trunk)
	} 
}

listTypes<-function(l){
	ltypes = sapply(X=l, FUN=function(x){ class(x) })
	return(ltypes)
}

settingsAsDataFrame<-function(settingsData){
	#filter out the metadata -- slots starting with a period, such as ".text"
	settingsData = settingsData[!grepl(pattern="^[.]", x=names(settingsData))]
	#inside the program the settings are stored as a list, transform this to a data.frame
	settingsDF = listToDf(lst=settingsData)
	return(settingsDF)
}

#the settings function collects user-entered settings or distributes previously defined user-entered settings
#using the settings object approch has several major advantages: 
#1) it collects all choices in one place, so that they can be easily passed, recorded and reported.
#2) it allows program to be run with previously selected settings
#3) it standardizes user-setting interface
#inputs: s: list object containing any settings already established.
#					prompt: the message string sent to the user as a prompt, also used to retreive values
#								
#output: list() object containg the settings and several other items:
#										slot $.text contains the most recent input by the user. This slot is not saved when settings are saved to a file
setting<-function(s, prompt, requireInput=T){
	if(is.data.frame(s)) s = dfToList(df=s)
# 	cat("\nsetting()\n")
	if(is.null(s$interactive)) s$interactive = F #if interactive flag is set to TRUE, a default value will be suggested
# 	cat(".")
	if(grepl(pattern="file", x=prompt, ignore.case=T)){
		return(settingFile(s, prompt, requireInput))
	}else{
		return(settingText(s, prompt, requireInput))
	}
}


test.settingFile<-function(){
	
	#create a file in advance, see if it will properly, quietly find it
	fname = "./testFile"
	write.table(x="testfile", file=fname)
	stest=list()
	prompt1 = "select the testFile.."
	stest[[prompt1]] = fname
	stest$interactive=FALSE
	
	t1res =settingFile(s=stest, prompt=prompt1)
	
	checkEquals(target=fname, current=t1res$.text)
	file.remove(fname)
	#this one should prompt becaue the file does not exist
	t1res =settingFile(s=stest, prompt=prompt1)
}

settingFile<-function(s, prompt, requireInput=T){
	if(prompt%in%names(s)){#if something has already been recorded for the prompt
		if(file.exists(s[[prompt]])){
			if(s$interactive=="TRUE"|!is.null(s$changeFiles)){#set to the character string "TRUE" so that the settings can survive file saves
				cat("\n", prompt, sep="")
				res = readline(paste("Enter s to select a file.\n Press enter to use the default file:\n", s[[prompt]]))
				while(T){
					if(res == ""){ #in this case use the default value
						break
					}else if(res=="s"){#in this case prompt to select a file
						while(T){
							s[[prompt]] = try(expr=file.choose(), silent=T)
							if(grepl(x=class(s[[prompt]]), pattern="error", ignore.case=T)){
								if("y"==readline("To skip this step enter y: ")) break
							}else{
								break
							}
							cat("\nPlease try again..\n")
						}
						break
					}
				}
			}
			s$.text = s[[prompt]]
			return(s)
		}else{
			if(s[[prompt]]!=""|s$interactive ) cat("\nThe file '", s[[prompt]], "' was not found, please select another.\n")
		}
	}
	#see if it's a prompt to select a file
	cat("\n",prompt,"\n")
	#it's a prompt
	while(T){
		fname = try(expr=file.choose(), silent=T)
		if(grepl(x=class(fname), pattern="error", ignore.case=T)){
			if(!requireInput) if("y"==readline("To skip this file selection step enter y")){fname=c();break} 
		}else{
			break
		}
		if( readline("\nWould you like to select a file? (y/n)\n")=="n"){
			fname=""
			break
		}
	}
	s[[prompt]] = fname
	s$.text = fname
	return(s)
}

settingText<-function(s, prompt, requireInput=T){
	if(prompt%in%names(s)){ #the value has already been set
		if(s$interactive){
			cat("\n",prompt,"\n")
			if(!is.null(s[[prompt]])) cat("Default value:\n", s[[prompt]])
			res = readline("Enter an option or just press enter to use the defalt value: ")
			if(res!="") s[[prompt]] = res
		}
		s$.text = s[[prompt]]
		return(s)
	}
	#if the value hasn't been set, the user will be prompted, then their input will be returned
	s[[prompt]] = readline(prompt)
	s$.text = s[[prompt]]
	return(s)
}

test.settingList<-function(){
	#test.settingList()
	s = list(interactive=F,YOmanwassup=c("option 1"))
	s[["YOmanwassup(Full list of options.)"]] = paste(c("1", "option 1"), collapse=";")
	#check that if the options have changed, the user will be re-prompted
	print("You should be prompted next")
	s = settingList(s=s, prompt="YOmanwassup",set=matrix(c("option 1","option 2"), ncol=1))
	resA = s$.text
# 	s = settingList(s=s, prompt="YOmanwassup choose one of those: ",set=c("option 1","option 2"))
	print("If you were not just prompted, the test failed. You should not be prompted now:")
	cat("\n\n next test:\n\n")
	print("You should not be prompted this time:")
	s = settingList(s=s, prompt="YOmanwassup",set=matrix(c("option 1","option 2"), ncol=1))
	resB = s$.text
	print("If you were just prompted, the test failed.")
	checkEquals(target=resA, current=resB)
	cat("\n\n next test:\n\n")
	print("Now trying with more sensical data, you should be prompted:")
	tm = cbind.data.frame(types=c("Missense_Mutation","Nonsense_Mutation"), counts=c(4,7), stringsAsFactors=F)
	prompt = "\nPlease enter the row numbers of the variant types you would like to analyze (sepparated by a space).\n"
	
	if(exists("somsets")){
		s=list(interactive=T)
		s[[prompt]] = somsets[[prompt]]
	}

	s = settingList(s=s, prompt=prompt, set=tm)	
	print("If you were not just prompted, the test failed.") 
	print(s)
}


#'@title settingList
#'@description Allows a user to select from a list of options.
#'@description Will prompt user with a list of options and assure one of the options is selected
#'@param s the settings list
#'@param prompt The message to be prompted to the user.
#'@param set data frame matrix or vector: the options the user has to select amongst
#'@return The settings list object is returned with the list of user selections in the prompt slot and the full list of options in a slot named prompt, '<prompt>(Full list of options.)'
#'@export
settingList<-function(s, prompt, set){

# 	print(s)
# 	print(set
	if(!is.matrix(set)&!is.data.frame(set)){
		if(is.vector(set)){
			set = matrix(data=set, ncol=1, dimnames=list(1:length(set), "Option descriptions"))
		}else{
			print("Error, illegal object provided for the 'set' argument")
		}
	}
	if(is.null(colnames(set))) colnames(set)<-"Option descriptions"
	cat("\nThese are the available options:\n")
	rownames(set) = 1:nrow(set)
	print(set)
	
	fullSet = paste(prompt, "(Full list of options.)", sep="")#this is provided so that
	
# 	print(fullSet)
# 	print(names(s)==fullSet)
# 	readline("That is the fullSet variable and it's equalities to the names of the items in the setting list")
	#the current options can be compared to the options available when a previous setting object was made
	#if the options have changed, the program will force prompt of the user
	
	#for when a list of objects is ultimately wanted and the set of numbers from a user might either be unstable between runs, 
	#or might be otherwise meaningless
	
	if(prompt%in%names(s)){
		#break appart the prompt into options

		s[[prompt]] = strsplit(s[[prompt]], split="; ")[[1]]
		
		#check if the options have changed
		forceInput = F
		if( !is.null(s[[fullSet]]) ){#if there's something stored for the fullset
			
			s[[fullSet]] = strsplit(x=s[[fullSet]], split="; ")[[1]]#split what's there

			if( sum( (!s[[fullSet]]%in%set[,1])) | sum(!set[,1]%in%s[[fullSet]]) ){
				#if there are any in the full set that aren't in the set or if there are any in the set that arent in the full set
				#if the options in the input set have changed, the user must be re-prompted
				
				cat("\nAvailable options have changed ")
				#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
				forceInput = T
			} 
		}
		
		if(!is.null(s[[prompt]]))     cat("\nDefault value:\n", s[[prompt]], "\n")
		#	s[[fullSet]] = paste(set[,1], collapse="; ")
		if( s$interactive|forceInput ){
			#display prompt
			cat("\n",prompt,"\n")
			
			if( !is.null(s[[prompt]]) ){
				
				line = readline("Enter your selection numbers sepparated by spaces, or just press enter to use the defalt value: ")
				while(T){#make sure everything the user enters can be used by the program
					if(line=="") break
					test = as.integer(strsplit(x=line,split=" ")[[1]])
					if( sum( !test %in% 1:length(set[,1]) ) ){
						cat("Sorry, one of the options in your entry\"", line,"\"was not an available option.\nPlease try again.\n")
					}else{
						break
					} 
				}
			}
			##process user input:
			if(line==""){#this will happen if the user selects the default values
				s$.text = s[[prompt]]
				s[[prompt]] = paste(s$.text, collapse="; ")
				s[[fullSet]] = paste(set[,1], collapse="; ")
				return(s)
			}else{#the user entered their own values
				selection = as.integer(strsplit(line, " ")[[1]])
				selection = set[selection,1]
				#make sure there are no NAs
				s$.text = selection[!is.na(selection)]
				s[[prompt]] = paste(s$.text, collapse="; ")
				s[[fullSet]] = paste(set[,1], collapse="; ")
				return(s)
			}
		}else{#if the mode is not interactive and a new option in set did not force input:
			s$.text = s[[prompt]] #strsplit(s[[prompt]], split="; ")[[1]]
			s[[prompt]] = paste(s$.text, collapse="; ")
			s[[fullSet]] = paste(set[,1], collapse="; ")
			return(s)
		}
	}#if(prompt is in the settings already available)

	#if this executes, there is not yet a prompt for the selection, prompt the user:
	while(T){
		cat("\n",prompt,"\n")
		line = readline("Enter your selection numbers or just press enter to use the default value: ")
		if(line!=""){
			selection = as.integer(strsplit(line, " ")[[1]])
			selection = set[selection,1]
		}
		
		#make sure there are no NAs
		s$.text = selection[!is.na(selection)]
		s[[prompt]] = paste(s$.text, collapse="; ")
		if(s[[prompt]]!="") break
	}
	s[[prompt]] = paste(s$.text, collapse="; ")
	s[[fullSet]] =  paste(set[,1], collapse="; ")
	return(s)
}#settingsList()