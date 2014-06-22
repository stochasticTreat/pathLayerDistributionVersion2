#runShiny.R

prepDrugSelect<-function(){
	uif = system.file("shinyDrugSelect/ui.R",package = "packageDir")
	serverf = system.file("shinyDrugSelect/server.R",package = "packageDir")
	source(uif)
	source(serverf)
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
	
	bfbTargDrugData=NULL
	
	if(exists("bfbTargDrugData")|!is.null(STUDY@results$drugSelectionWorksheet)){
		print("option1")
		if(is.null(bfbTargDrugData)) bfbTargDrugData <<- makeDrugSelectionWorksheet(STUDY=STUDY)
	}else{
		print("option2")
		bfbTargDrugData <<- makeDrugSelectionWorksheet(STUDY=STUDY)
	}

	prepDrugSelect()
	
	shinyDir = dirname(system.file("shinyDrugSelect/server.R",package = "packageDir"))
	
	runApp(shinyDir)#"./shinyDrugSelect/")
	
	return(bfbTargDrugData)
}



