#runShiny.R

prepDrugSelect<-function(){
	uif = system.file("shinyDrugSelect/ui.R",package = "packageDir")
	serverf = system.file("shinyDrugSelect/server.R",package = "packageDir")
	source(uif)
	source(serverf)
}

runShinyMain<-function(STUDY=STUDY){
	print("inside runShinyMain")
	library("shiny")
	# setwd("..")
	
	bfbTargDrugData=NULL
	
	if(exists("bfbTargDrugData")){
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



