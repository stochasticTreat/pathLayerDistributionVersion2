# #processOHSUdat
# 
# # source('./processOHSUseqDat_functions.R')
# 
# #processOHSUseqDat_functions.R
# # source('./acc_functions.R')
# # 
# # if(!exists("acc_loaded")){
# # 	source("acc_functions.R")
# # }
# 
# if(!exists("path_detail")){
# 	path_detail = getPaths()
# }
# 
# 
# if(!exists("study_name")){
# 	study_name = gsub(pattern=":",replacement=".",x=paste("single_run_of_Sequence_Capture_data",Sys.time(), sep="_"))
# }
# 
# if(!exists("Sequence_CaptureSummary")) Sequence_CaptureSummary = NULL
# 
# #selection structure
# cat("\n\nWelcome to the Sequence_Capture data interface\n\n")
# while(T){
# 	cat("\n\nMain options for Sequence_Capture data interface:\n\n")
# 	mainsel = readline(paste("To process and analyze a Sequence_Capture data set enter p\n",
# 													 "To make HTML summary of Sequence_Capture data enter s\n",
# 													 "To exit Sequence_Capture interface enter e\n",sep=""))
# 	if(mainsel=="p"){
# 
# 		opfile = filePrompt(defaultfile="./input/AML_corrolated_overlap/AMLSeqCapOverlapPatOnly.txt")
# 
# 		Sequence_CaptureSummary = processSequenceCaptureData(seqfname=opfile,verbose=F,
# 																									 paths_detail=path_detail, 
# 																									 study_name=study_name)
# 	}else if(mainsel=="s"){
# 		if(!is.null(Sequence_CaptureSummary)){
# 			htmlSummary(sumset=Sequence_CaptureSummary$enrichmentSummary,
# 									fname=paste("./output/Summary_Of_path_enrichment_Using_Sequence_Capture_data_for",study_name,Sys.Date(),".html",sep=""))
# 			htmlSummary(sumset=Sequence_CaptureSummary$coverageSummary,
# 									fname=paste("./output/Summary_Of_Sequence_Capture_Coverage_for",study_name,Sys.Date(),".html",sep=""))
# 		}else{cat("\nAn analysis of Sequence_Capture data hasn't been conducted since this program\nwas loaded, so a summary of that analysis can't really be saved.\n")}
# 	}else if(mainsel=="e"){
# 		break
# 	}else{
# 		print("Sorry, program did not understand that input. Please try again.")
# 	}
# 	
# }
# 
