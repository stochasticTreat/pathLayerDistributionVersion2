cat("\nThis script installs needed dependencies from bioconductor.\n")

neededPackages = c('reactome.db', 'RCytoscape', 
									 'graphite', 'HGNChelper', 
									 'VennDiagram', 'calibrate', 'biomaRt', 
									 'xtable', 'hwriter', 'hwriterPlus', 'graph', 'plyr', 'rBiopaxParser', 'ggplot2', 
									 'gridExtra', 'shiny')

source("http://bioconductor.org/biocLite.R")
biocLite("RCytoscape")

for(pkg in neededPackages) biocLite(pkg)
