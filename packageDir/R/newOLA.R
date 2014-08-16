overlapmerge<-function(settings, study){
	
	keycols = c("path_id","full_path_length")
	#get all path results
	allMerged = NULL
	armNames = names(study@results)[grepl(pattern="_summary", x=names(study@results))]
	pathLists = list()
	
	for(n in armNames){
		
		#extract the path summary
		cur = study@results[[n]]$pathsummary
		curTitle = gsub(pattern="_summary", replacement="", x=n)
		
		#add to the pathlists
		pathLists[[curTitle]] = cur$path_id
		#adjust the column names
		colnames(cur)<-setColumnPrefix(coln=colnames(cur), 
														 prefix=curTitle, 
														 excols=keycols)
		
		#conduct the merge on the main data set
		if(is.null(allMerged)){#if it's the first
			allMerged = cur
		}else{
			#merge them
			allMerged = merge(x=allMerged, y=cur, by=keycols, all=TRUE)
		}
		
		#check for/handle coverage data
		if(notNullOrEmpty(study@results[[n]]$coverage_summary)){
			curcovpaths = study@results[[n]]$coverage_summary$pathsummary
			curCoverageTitle = paste0(gsub(pattern="_summary", replacement="", x=n),"_coverage")
			
			#add to pathlists
			pathLists[[curCoverageTitle]] = curcovpaths$path_id
			
			colnames(curcovpaths)<-setColumnPrefix(coln=colnames(curcovpaths), 
																		 prefix=curCoverageTitle, 
																		 excols=keycols)
			allMerged = merge(x=allMerged, y=curcovpaths, by=keycols, all=TRUE)
		}
	}
	overlapVenn(plist=pathLists, vennTitle="Overlaps of paths from different sources")
	
}


overlapVenn<-function(plist, vennTitle){
	names(plist)<-gsub(pattern="_",replacement=" ", x=names(plist), fixed=TRUE)
	venn(plist, simplify=T)
	# 	mtext(c("overlaps"),side=3,line=2)
	title(vennTitle)
}

notNullOrEmpty<-function(checkVal){
	if(class(checkVal)=="NULL") return(FALSE)
	
	if(is.vector(checkVal)){
		return(length(checkVal)>0)
	}
	
	if(is.data.frame(checkVal)|is.matrix(checkVal)){
		return(nrow(checkVal)>0)
	}
	message("Unchecked type:",class(checkVal)," returning True")
	return(TRUE)
}

#sets the name prefix for a data frame
#param coln the column names
#param prefix the prefix to be added to the column names
#excols the column that will not be modified
#return the adjusted set of column names
setColumnPrefix<-function(coln, prefix, excols){
	#set prefix
	changei = !coln%in%excols
	coln[changei]<-paste0(prefix,"__",coln[changei])
	return(coln)
}