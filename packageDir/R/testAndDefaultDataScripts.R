

#'@title Get an example patient gene matrix
#'@description Retreives from package data a patient gene matrix with corresponds with the Reactome pathway, 'Abacavir metabolism' being completely enriched in active (mutated, drug sensitive, etc.) genes.
#'@return \code{matrix} object with column names set as patient IDs, row names set as HUGO gene symbols, and \code{logical} cell values indicating if the corresponding gene is active in the corresponding patient.
#'@export
getTestPGM<-function(){
	#the pgm used is the abacavir metabolism pgm from the moc somatic mutation data
	#pgm=STUDY@results$somatic_mutation_aberration_summary$patientGeneMatrix
	#save(pgm, file="./testData/pgm.rda")
	defname = system.file("testData/pgm.rda", package = "packageDir")
	# 	load(file="./testData/pgm.rda", verbose=T)
	load(file=defname, verbose=T)
	return(pgm)
}