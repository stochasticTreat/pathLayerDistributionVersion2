#main interface
########################################################################
#initiate peripheral library scripts
source('./main_functions.R')
#load accessory functions
if(!exists("acc_loaded")){
	source("./acc_functions.R")
}
source('./drugDbInterface.R')
# source("http://bioconductor.org/biocLite.R")

#accessory functions
if(!require("RUnit")){
	install.packages("RUnit")
	library("RUnit")
}
if(!require("tcltk")){
	install.packages("tcltk")
	library("tcltk")
}

acc_loaded=T
source('./corsym.R')
source('./summaryTable5.R')
source('./pathway_functions.R')
source('./InitiateDataStructures.R')
source('./OverlapAnalysisFunctions.R')
source('./path_paint.R')
source('./save_and_load_data.R')
source('./SettingsObjectInterface.R')

# source("http://bioconductor.org/biocLite.R")
# biocLite("RCytoscape")
library("RCytoscape")
library("graphite")
