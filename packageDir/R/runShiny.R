#runShiny.R

prepDrugSelect<-function(){
	if( !file.exists("shinyDrugSelect/ui.R") ){
		uif = system.file("shinyDrugSelect/",package = "packageDir")
		file.copy(from=uif, to="./", recursive=TRUE)
	}
	source("./shinyDrugSelect/ui.R")
	source("./shinyDrugSelect/server.R")
	return("./shinyDrugSelect/")
}

#'@title Run the drug selection worksheet.
#'@description Running this function will open up an HTML page allowing interactive selection and exploration of targeted drugs. Information is supplied to link genes to drugs and filter by clinical trial phase, gene mutation in a cohort and cellular pathways. 
#'@param STUDY A \code{Study} object into which experimental data has been loaded. 
#'@return The \code{data.frame} displayed in the drug selection worksheet, relating drug targets, drugs, clinical trial phase, and other data. 
#'@export
#'@import shiny
runDrugWorksheet<-function(STUDY=STUDY){
	print("inside runShinyMain")
	# library("shiny")
	# setwd("..")
	buildNewWorksheet = T
	#if bfbTargDrugData exists and has rows in it
	if(exists("bfbTargDrugData")){
		if(!is.null(bfbTargDrugData)){
			buildNewWorksheet = "Y" == toupper(readline("A drug selection worksheet was found. \nWould you like the program to update the worksheet with the current study data?\n(enter Y or N ) "))
		}
	}
		#ask if it should be rebuilt
	
	if(buildNewWorksheet){
		cat("Building new drug selection worksheet.\n")
		bfbTargDrugData = makeDrugSelectionWorksheet(STUDY=STUDY)
	}else{
		cat("\nUsing old drug selection worksheet\n")
		bfbTargDrugData = STUDY@results$drugSelectionWorksheet
	}
	
	bfbTargDrugData <<- bfbTargDrugData
	
	# 	else{
	# 		print("option2")
	# 		bfbTargDrugData = makeDrugSelectionWorksheet(STUDY=STUDY)
	# 		bfbTargDrugData <<- bfbTargDrugData
	# 	}
	shinyDir = try(expr={prepDrugSelect()}, silent=T)
	
	#dirname(system.file("shinyDrugSelect/server.R",package = "packageDir"))
	if(grepl(pattern="error", x=shinyDir, ignore.case=TRUE)){
		mess = paste("Error, could not find directory ./shinyDrugSelect/ in current working directory,",getwd(),
									".\nPlease make sure R has write access to this directory so that this program can copy data to it, then try again.")
		message(mess)
		warning(mess)
	}else{
		runApp(shinyDir)
	}
	# "./shinyDrugSelect/")
	# cat("bfbTargDrugData exists:",exists("bfbTargDrugData"),"\n")
	#	cat("bfbTargDrugData is null:",is.null(bfbTargDrugData),"\n")
	# print(getwd())
	
	return(bfbTargDrugData)
}



