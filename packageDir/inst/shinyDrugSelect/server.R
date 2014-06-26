#server.R
require(shiny)

newbfbTargDrugData=NULL

shinyServer(function(input, output){
	# 	
	# 	output$exitHTML<-reactive({
	# 		
	# 	})
	# 	output$myImage<-renderImage({
	# 		list(src="test.png", 
	# 				 contentType='image/png', 
	# 				 width=400, 
	# 				 alt="this is only a test image")
	# 	})
	newbfbTargDrugDataFunction<-reactive({
		print("Adjusting pathways")
		inFile = input$file1
		if(!is.null(input$file1)){
			limPaths = read.csv(inFile$datapath, header=F, sep="\t", quote="", stringsAsFactors=F)
			print(getwd())
			setwd("..")
			print(getwd())
			bfbTargDrugData <<- makeDrugSelectionWorksheet(STUDY=STUDY, plimit=limPaths[,1])
			setwd("./shinyDrugSelect/")
			displayData = bfbTargDrugData
		}else{
			displayData = bfbTargDrugData
			bfbTargDrugData<<-bfbTargDrugData
		}
		# 		newbfbTargDrugData<<- displayData
		
		return(displayData)
		
	})

	output$targetGeneData <- renderDataTable({
		
		inFile<-input$file1
		
		tmp = newbfbTargDrugDataFunction()
		
		if(is.null(tmp)){
			displayData = bfbTargDrugData
		}else{
			displayData = tmp
		}
		showRows = rep(T, times=nrow(bfbTargDrugData))
		if(input$exclude!=""){
			showRows = !packageDir:::rowsContaining(df=displayData, rx=input$exclude)
		}
		
		displayData[showRows,input$show_vars,drop=FALSE]
	})
	
	observe({
		if(input$exit == 0) return()
		stopApp()
	})
	
	output$downloadData <- downloadHandler(
		filename = function(){
			paste0('drugTargetingSpreadsheet', Sys.Date(), '.txt')
		},
		content = function(file) {
			write.table(x=bfbTargDrugData[,input$show_vars,drop=FALSE], file, sep="\t", row.names=F, col.names=T)
		}
	)
	
	# 	output$dynamicRanges<-renderUI({
	# 		
	# 	})
})


# , sanitize.text.function = function(x) x