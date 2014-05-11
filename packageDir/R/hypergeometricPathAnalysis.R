
#hypergeometricPathEnrichment

#'@title hypergeometricPathEnrichment
#'@description provides hypergeometric test for analysis of pathway enrichment
#'@param paths_detail: an object of the Path_Detail class
#'@param pathSigRunner: a object of the PathSummaryRunner reference class
#'@return data.frame with two columns, one with hypergeometric p-values for each path, the other with B.H./FDR adjusted p-value.
#'@note Implements the path test interface with one input, the path summary list, and the output, a data.frame with pathway names as row name(s) and named column(s) giving the result of the significance test. 
#'@export
hypergeometricPathEnrichment <- function (pathSigRunner, paths_detail) {
	
	targetMatrix=pathSigRunner$.targetMatrix
	
	if(!nrow(targetMatrix)) return(cbind.data.frame(hyperg_p_value=1,hyperg_p_w_FDR=c(1)))
	cat(" Hypergeometric test .. ")
	hyperg_p_value = rep(0,nrow(targetMatrix))
	total_active = ncol(targetMatrix)
	total_nodes = ncol(paths_detail$paths)
	for(i in 1:nrow(targetMatrix)){
		curpath = rownames(targetMatrix)[i]
		active_in_path = sum(targetMatrix[curpath,])
		total_in_path = sum(paths_detail$paths[curpath,])
		hyperg_p_value[i] = 1-phyper(active_in_path, 
																 total_active, 
																 total_nodes-total_active, 
																 total_in_path)
	}
	dout = cbind.data.frame(hyperg_p_value, hyperg_p_w_FDR=p.adjust(hyperg_p_value, method=c("BH")))
	rownames(dout)<-rownames(targetMatrix)
	return(dout)
}

pathAnalysisFunctions = list(hypergeometric=hypergeometricPathEnrichment)

