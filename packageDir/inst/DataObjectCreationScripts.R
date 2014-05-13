pmd = "/Users/samhiggins2001_worldperks/tprog/distribution/pathLayerDistributionVersion2/reference_data/paths/pathMetaData.txt"
fullTab = read.table(file=pmd, header=T, sep="\t", stringsAsFactors=F, comment.char="")

#not using this approach; insted placing the map table in extdata/

read.csv(system.file("extdata/pathMetaData.txt", 
										 package = PACKAGENAME), as.is = TRUE)