require("roxygen2")

# dir.create("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/executionDir")
setwd("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/executionDir")

detach("package:packageDir", unload=TRUE)
remove.packages("packageDir")
roxygenize(overwrite=T, 
					 package.dir="/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir/")

# roxygenize()

readline("please run R CMD build packageDir then press enter")
library("packageDir", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
# cd ~/tprog/distribution/pathLayerDistributionVersion2

# install.packages("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir_0.1.tar.gz", repos = NULL, type = "source")

# devtools::install("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir")
library("packageDir")
STUDY = allInteractiveMainFunction()

library("devtools")
load_all(pkg="packageDir", export_all=T, reset=T)

packageDir:::list_to_table()


abacavirIntegrationTest = STUDY



# PACKAGENAME="packageDir"
# data("defaultSummaryTable140504")
# system.file("extdata/defaultSummaryTable140504.rdata", package = "packageDir")
# system.file("extdata/mog_map.csv", package = "HGNChelper")
# remove.packages("HGNChelper")
# install.packages("~/tprog/distribution/HGNChelper_0.3.0.tar.gz", repos = NULL, type = "source")

#'@examples
#' \dontrun{
#' #for basic interactive usage:
#' runInteractivePathAnalysis() 
#' #for interactive usage with additional data input arm(s) provided
#' runInteractivePathAnalysis(additionalArms=function(stud){
#'																										stud@arms = loadDataArm(description="Load drug screen data to show addition of a study arm",
#'																																				title="functional_drug_screen_summary_testAdd", 
#'																																				mainFunction=RunDrugScreen, 
#'																																				arms=stud@arms)
#'																										return(stud)
#'																									})}

sourceAllInFolder<-function(folname="../packageDir/R/"){
	
	if(!grepl(pattern="/$", x=folname)) folname = paste0(folname,"/")
	
	fnames = dir(folname)
	fnames = fnames[grep(pattern="r", x=fnames, ignore.case=T)]
	fnames = paste0(folname, fnames)
	for(fn in fnames){
		cat("..",fn,"..")
		source(fn)
	}
	
}
sourceAllInFolder()

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
												 decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
		fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
		as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.dim)
	names(out) <- c("Type", "Size", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}
# shorthand
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
