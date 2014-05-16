#main interface
########################################################################
#initiate peripheral library scripts

# #load accessory functions

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

if(!require(ggplot2)){
	install.packages("ggplot2")
	if(!require(ggplot2)){
		print("error, could not install ggplot2!!!!!!!!!!!!!!!!!")
	}
}

if(!require("graphite")){
	install.packages("graphite")
	library("graphite")
}
if(!require(graphite)){
	print("Cannot load graphite.. . \nthat package is needed for the interface to cytoscape to work")
}

if (!require(RCytoscape))
	stop("the RCytoscape package is missing")

if (!require(RCytoscape))
	stop("the RCytoscape package is missing")

if(!require("VennDiagram")){
	install.packages("VennDiagram")
	library("VennDiagram")
}

if(!require(calibrate)){
	cat("\nAttempting to install package calibrate which allows better plotting options.\n")
	
	if(!require(calibrate)){
		install.packages("calibrate")
		library(calibrate)
		cat("\nError: could not install package calibrate!!!\n")
	}
}

if (!require(RCytoscape)) 
	stop("the RCytoscape package is missing")

#1 initilize
require("RCurl")

if(!require("biomaRt")){
	source("http://bioconductor.org/biocLite.R")
	biocLite("biomaRt")
	library("biomaRt")
}

if(!require("xtable")){
	print("Trying to install xtable so that HTML output can be generated")
	install.packages("xtable")
	if(require("xtable")){
		print("xtable installed and loaded")
	} else {
		stop("could not install xtable")
	}
}

library("xtable")
if(!require("hwriterPlus")){
	print("Trying to install hwriterPlus so that HTML output can be generated")
	install.packages("hwriterPlus")
	if(require("hwriterPlus")){
		print("hwriterPlus installed and loaded")
	} else {
		stop("could not install hwriterPlus")
	}
}

library("tools")

library("graphite")
library("rBiopaxParser")


# biocLite("RCytoscape")
library("RCytoscape")
library("graphite")

