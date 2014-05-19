require("roxygen2")

# dir.create("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/executionDir")
setwd("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/executionDir")

roxygenize(overwrite=T, package.dir="/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir/")
roxygenize()
detach("package:packageDir", unload=TRUE)
remove.packages("packageDir")
readline("please run R CMD build packageDir then press enter")
# cd ~/tprog/distribution/pathLayerDistributionVersion2
#

install.packages("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir_0.1.tar.gz", repos = NULL, type = "source")

# devtools::install("/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/packageDir")
library("packageDir")
STUDY = runInteractivePathAnalysis()

library("devtools")
load_all("packageDir")
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
#'																									

