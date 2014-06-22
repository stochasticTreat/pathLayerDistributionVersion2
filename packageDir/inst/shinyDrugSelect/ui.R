#ui.R
#the UI for finding drugs
require(shiny)
shinyUI(pageWithSidebar(
	headerPanel(title='Drug selection worksheet'),
	sidebarPanel(
		h4("Showing data for drugs which target dark pathways"),
		br(),
		actionButton("exit", "Exit web interface"),
		br(),		
		downloadLink(outputId='downloadData', label='Download spreadsheet'),
		br(),
		fileInput('file1', 'Limit pathways to those named in a .txt file',
							accept=c('text/csv', 'text/comma-separated-values,text/plain', '.txt')), 
		
		checkboxGroupInput(inputId="show_vars",
											 selected=c("Gene symbol", 
											 					 "Drug name", 
											 					 "Drug Type", 
											 					 "Number of dark paths containing gene",
											 					 "Names of paths containing gene",
											 					 "clinical trial IDs", 
											 					 "Total targets", 
											 					 "New Targets", 
											 					 "Aberrations in gene, across cohort"),
											 choices=colnames(bfbTargDrugData),
											 label="Columns of data to display")
	),
	mainPanel(
		h5("The selection boxes below use \"regular expressions\", thus use the pipe symbol \"|\" to sepparate \"or'd\" terms (ex: use \"Phase-4|Phase-3\" to get both phase 4 and phase 3 clinical trials). \"[1-9]+\" can be used to get all numbers greater than 0."),
		textInput(inputId="exclude", 
							value="",
							label="Input terms you would like to exclude from rows. \nSeparate multiple terms with the pipe character (\"|\") ex: GABA|Vitamin C excludes rows with GABA or Vitamin C."),
							
		dataTableOutput('targetGeneData')
# 		imageOutput('myImage')
	)
))

