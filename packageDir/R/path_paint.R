#paint_path.R


# #this function takes patient gene matracies, 
# #then, via the graphite and RCytoscape packages, 
# #outputs a path to cytoscape
# attributeDataSet<-function(g, attributeName){
# 	
# 	#first get the names of the nodes
# 	
# 	#nodeData(<graph_object>, <node_name>, <attribute_name>) = value
# 	nodeData(g, "A", "lfc") = -3
# 	nodeData(g, "B", "lfc") = 0
# 	nodeData(g, "C", "lfc") = 3	
# }

is.error.message<-function(biopaxFile){
	if(grepl(pattern="error$",x=class(biopaxFile), ignore.case=T)){
		return(T)
	}
	return(F)
}

fromReactome<-function(study){
	res = grepl(pattern="reactome", x=study@studyMetaData@paths$source, ignore.case=T)
	return(res)
}

clean.biopaxRecords<-function(pwrecord){
	
	pwr = unique(pwrecord)
	#check if there are duplicate path names
	#modify the oldest to include "duplicated
	
}


#'@title Add network diagrams of affected pathways. 
#'@description Outputs to cytoscape diagrams of affected pathways . 
#'@param study A \code{Study} object. 
#'@param limitCol The column from the results that is used to limit the number of path diagrams displayed. 
#'@param limitVal The value in the limitCol used to limit the number of pathways. Net work diagrams will only be produced for pathways with values in limitCol which are smaller than limitVal. 
#'@return A \code{Study} object will be returned with file names of network diagrams added to the imageSlots slot of the appropriate results sets. 
#'@import RCytoscape
#'@import rBiopaxParser
#'@import graphite
#'@import graph
#'@export
addPathwayImagesWithSelection<-function(study, 
																				limitCol="hyperg_p_w_FDR", 
																				limitVal=0.05){
	study_name = study@studyMetaData@studyName
	#require(RCytoscape)
	#require(rBiopaxParser)
	results = study@results
	path_detail = study@studyMetaData@paths
	# 	if(is.null(path_detail$graphite)){
	# 		results = pathwaysFromBiopax(study)
	# 
	# # 		cat("\n!!! sorry, at this time pathway images are only available\nif pathways from the Graphite package are used.")
	# 		return(results)
	# 	}
	
	abi = grep("aberration|overlap", x=names(results))
	options = matrix(1:length(abi), dimnames=list(names(results[abi]),"option number"))
	print(options)
	
	sel = readline("Please enter the number corresponding to the data set you want to add pathway images for.\nEnter \"all\" to add path images for all aberration data (takes a long time).\n")
	
	# 	draw_network_prompt = readline("To display cytoscape network diagram of pathway overlap(s) press enter.\n(must have Cytoscape turned on and cytoscape plugin \"cytoscapeRPC\" installed and activated to display diagram) ")
	
	if(sel=="all"){
		resSetName = NULL
	}else if(!is.na(as.numeric(sel))){
		resSetName = rownames(options)[as.numeric(sel)]
	}else{
		readline("\nSorry, could not understand input.\nPress any key to continue.")
		return(results)
	}

		
	#pull out the needed pathway names
	
 if(pathsAreFromGraphite(study)){
		
		if(!length(path_detail$graphite)) path_detail$graphite = .loadGraphitePaths(study)

		if(resSetName=="overlap_analysis"){
			results = paintOverlap(results=results, study_name=study@studyMetaData@studyName, paths_detail=path_detail)
		}else{
			dir.create(path=paste("./output/", study_name,"/",resSetName,"/","imageSlots/", sep=""),recursive=T,showWarnings=F)
			results=addPathwayImages(study_name=study_name, 
															 results=results, 
															 path_detail=path_detail,
															 resSetName=resSetName,
															 sigtest=limitCol, 
															 siglimit=limitVal)
		}
		
	 }else	if(fromReactome(study)){
	 	
	 	pathnames = extractPathNames(resSetName=resSetName, study=study, limitCol=limitCol, limitVal=limitVal)
	 	notAvailableBiopax = pathwaysFromBiopax(study=study, pathNames=pathnames, resSetName=resSetName)
	 	
	 }
	
	#system('/usr/bin/afplay ./reference_data/Submarine.aiff')
	study@results = results
	return(study)
}#addPathwayImagesWithSelection()

pathsAreFromGraphite<-function(stud){
	return(sum(grepl(pattern="graphite", x=stud@studyMetaData@paths$source, ignore.case=T))>0)
}

#'@title load graphite-supplied pathways
#'@description Internal function used to ensure graphite supplied network diagrams are loaded and made available for cell network display.
#'@param stud The \code{Study} object
#'@return The pathway repository selected.
.loadGraphitePaths<-function(stud){

	repset = list(Reactome=reactome, NCI=nci, KEGG=kegg)
	
	selPaths = repset[[stud@studyMetaData@paths$name]]
	
	return(selPaths)
	
}


#'@title Get names of significant pathways.
#'@description Extracts path names given a results set name, study and limitCol
#'@param resSetName the name of a summary table results set
#'@param study a study object
#'@param limitCol the column the returned pathways should be limited by
#'@param limitVal the value in the limitCol the returned pathways should be limited to
#'@return vector, a list of pathways from study meeting the parameters described in the function argument
#'@export
extractPathNames <- function (resSetName, study, limitCol, limitVal) {

		resSet = study@results[[resSetName]]
		if(resSetName=="overlap_analysis"){
			pathNames = resSet$"Aberration enriched, not drug targeted"[resSet$"Aberration enriched, not drug targeted"[,limitCol]<limitVal,1]
			pathNames = c(resSet$"Aberration enriched, containing sensitive targets", pathNames)
		}else{
			pathNames = resSet$pathsummary[resSet$pathsummary[,limitCol]<limitVal,1]
		}
		plens = study@studyMetaData@paths$full_path_length[pathNames,]

		
		if(sum(plens>=600)){
			while(T){
				res = try({			hist(plens,breaks=40, 
													 xlab="Path length", 
													 ylab="Number of paths",
													 main="Distribution of path lengths\nfor selected subset of pathways")}, silent=T)
				if(!is.error(res)){
					break
				}
				tmp = readline("Please increase the plot window size, then press enter..")
			}
			cat("\nOne or more very large pathways (>600 genes) were found.\nDisplay of these pathways takes a long time.\n")
			skipLargePathways = "y" == readline("Would you like to skip display of large pathways?\n(please enter y or n)")
			if(skipLargePathways){
				plens = study@studyMetaData@paths$full_path_length[pathNames,]
				pleni = plens<600
				cat("\nOnly displaying this set of pathways:\n")
				print(matrix(data=plens[pleni], ncol=1, dimnames=list(names(plens)[pleni], "path length")))
				cat("\nSkipping display of these pathways:\n")
				print(matrix(data=plens[!pleni], ncol=1, dimnames=list(names(plens)[!pleni], "path length")))
				pathNames = names(plens)[pleni]
			}
		}

		
		return(pathNames)	
}

addPathwayImages<-function(results, 
													 study_name,
													 path_detail=path_detail, 
													 sigtest="hyperg_p_w_FDR", 
													 siglimit=0.05,
													 resSetName=NULL){
	#addPathwayImages()
	#creates pathway images for all aberrational data sets and makes them available as files for html outputs
	if(!is.null(resSetName)){
		tmpRes = results
		tmp = list()
		tmp[[resSetName]] = results[[resSetName]]
		results = tmp
	}
	#tmp = readline("Please assure cytoscape is on, and that\nthe CytoscapeRPC plugin has been activated.\nPress enter to continue")
	
	#check if imageSlots folder already exists
	relativePath_base = paste("/output/",study_name, sep="")
	graph_set = path_detail$graphite
	
	#for each index with aberration data, pull out all pathways and build path diagrams
	abi = grep(pattern="aberration",x=names(results))
	for(i in abi){	
		#declare/intitilize variables
		pfiles = c()
		pathNames = c()
		
		curSumName = names(results[i])
		cat("\n*************",curSumName,"*************\n")
		cursummary = results[[i]] #pull out an aberration summary
		#establish directory path
		relativePath = paste(relativePath_base,names(results)[i],"imageSlots/", sep="/")
		dir.create(paste(".", relativePath),recursive=T,showWarnings=F)
		#find the set of pathways enriched, as to the hypergeometric, 
		sigpaths = cursummary$pathsummary[,sigtest]<siglimit
		path_names = cursummary$pathsummary$path_id[sigpaths]
		#use these pathway names to build the networks using the paintPathway function
		for(p in path_names){
			cat("\nConverting pathway node identifiers...")
			gpath = convertIdentifiers(pathway=graph_set[[p]], type="symbol")
			if(length(nodes(gpath))>1){
				print(p)
				print(length(nodes(gpath)))
				diagram = paintPath(graphitePathway=gpath, ab_analysis=cursummary, paths_detail=path_detail)
				#sanatize pathname
				clean_p_name = gsub(pattern="[/;:*.~<>]", replacement="_", x=p)
				#make the image file name
				imageFileName = paste(clean_p_name,"_path_image2.png",sep="")
				#construct the absolute name that saveImage() needs
				diagram_fname = paste(getwd(),relativePath,imageFileName,sep="")
				saveImage(obj=diagram,image.type="png",file.name=diagram_fname, scale=2)
				#results[[i]][["imageSlots"]][[p]] = imageFileName
				pfiles = c(pfiles, imageFileName)
				pathNames = c(pathNames, p)
			}else{
				pathNames = c(pathNames, p)
				pfiles = c(pfiles, paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated."))
				# 				p = paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated.")
				# 				results[[i]][["imageSlots"]][[p]] = paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated.")
			}
		}#for each significant path
		imageSlots = matrix(data=pfiles,ncol=1,dimnames=list(pathNames=pathNames))
		results[[curSumName]][["imageSlots"]] = list()
		results[[curSumName]][["imageSlots"]][["imageSlots"]]= imageSlots
	}#for each aberration record
	
	if(!is.null(resSetName)){
		tmpRes[[resSetName]] = results[[1]]
		results = tmpRes
	}
	return(results)
}


paintOverlap<-function(results, study_name, paths_detail, specificPath=NULL){
	
	cat("\nOutputting network diagrams for pathway overlap analysis..\n")
	if(is.null(results$overlap_analysis)){
		cat("\n!!! Must run over lap analysis first (Choose \"Examine drug screen and aberrational pathway overlap\" from the main menu)\n")
		return(results)
	}
	graph_set = paths_detail$graphite
	
	if(is.null(specificPath)){
		### from overlap table, get overlapping pathways
		overlapTargeted = results$overlap_analysis[["Aberrationally enriched, containing drug targets"]]$path_id #this looks like a whole table, not just a list of pathway names
		# take path_ids out of the table
		overlapSensitive = results$overlap_analysis[["Aberration enriched, containing sensitive targets"]] #perhaps this is just a set of pathway names?
		# it comes from this: outlist[["Aberration enriched, containing sensitive targets"]] = as.matrix(conSensAndAb, ncol=1)
		if(!is.vector(overlapSensitive)){
			overlapSensitive = overlapSensitive[,1]
		}
		overlapSet = union(overlapSensitive, overlapTargeted)
	}else{
		overlapSet = specificPath
	}
	
	pfiles = c()
	pathNames = c()
	
	if(!length(overlapSet)) return(results)
	
	### establish directories for images
	relativePath_base = paste("/output/",study_name,"/overlap_analysis/", sep="")
	relativePath = paste(relativePath_base,"imageSlots/", sep="")
	relativewithdot = paste(".", relativePath,sep="")
	dir.create(relativewithdot,recursive=T,showWarnings=F)
	
	### take the overlapping pathways and send their records to paint path
	
	for(p in overlapSet){
		cat("Current path: ")
		print(p)
		gpath = convertIdentifiers(pathway=graph_set[[p]], type="symbol")
		if(length(nodes(gpath))>1){
			diagram = paintPath(graphitePathway=gpath, paths_detail=paths_detail,
													drug_analysis=results$overlap_analysis$functional_drug_screen_summary, 
													ab_analysis=results$overlap_analysis$combined_aberrations_summary)
			###save path image/deal with naming
			clean_p_name = gsub(pattern="[/;:*.~<>]", replacement="_", x=p) #sanatize pathname
			imageFileName = paste(clean_p_name,"_path_image2.png",sep="") #make the image file name
			diagram_fname = paste(getwd(),relativePath,imageFileName,sep="") #construct the absolute name that saveImage() needs
			saveImage(obj=diagram,image.type="png",file.name=diagram_fname, scale=2)
			pfiles = c(pfiles,imageFileName)
		}else{
			#p = paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated.")
			#dataSummary[["imageSlots"]][[p]] = paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated.")
			pfiles = c(pfiles,paste("The pathway \" ",p,"\"only has one node, thus the diagram of this pathway is obviated"))
		}
	}
	imageSlots = matrix(data=pfiles,ncol=1,dimnames=list(pathNames=overlapSet))
	#save the results	to the local file structure
	results$overlap_analysis[["imageSlots"]] = list()
	results$overlap_analysis[["imageSlots"]][["imageSlots"]]= imageSlots
	return(results)
}#paintOverlap

paintPath<-function(graphitePathway, ab_analysis, paths_detail, drug_analysis=NULL){
	
	
	### set up variables
	pathName = graphitePathway@title #used to extract data from path matrix
	cat("\nPath name:", pathName,"\n")
	cat("Nodes in pathway:",length(graphitePathway@nodes),"\n")
	pathway = graphitePathway #used to set up graph data structure
	paths = paths_detail$paths
	print("checking gene symbols are all official..")
	### correct symbols
	# 	pathway@nodes = corsym(symbol_set=pathway@nodes,
	# 												 verbose=F,
	# 												 symref=paths_detail$symtable)
	# 	print("gene symbols are all official")
	#first get the nodes of the pathway
	pnodes = colnames(paths)[paths[pathName,]]
	
	### set up aberration color vector
	#now, make parallel vectors for the pnodes vector, using the values in ab_analysis and drug_analysis
	#ab vector first
	abs_logic_vector = pnodes%in%rownames(ab_analysis$genesummary)
	#make color vector from logic vector
	abs_color_vector = rep("#e5e5e5",times=length(pnodes)) #set white as the default, non ab color
	#e5e5e5 is grey
	abs_color_vector[abs_logic_vector] = "#fbff00"#set yellow as the ab color
	#fbff00 is yellow
	names(abs_color_vector) = pnodes
	
	#make a color vector from the abs logic vector
	
	print("building graphNEL object")
	
	g <- buildGraphNEL(nodes=nodes(pathway), 
										 edges=edges(pathway), 
										 sym=FALSE)#, ...)
	nodes(g)<-corsym(symbol_set=g@nodes, symref=paths_detail$symtable,verbose=F)
	cat(".")
	g <- markMultiple(g)#marks edgeType attribute
	cat(".")
	g <- initEdgeAttribute(graph=g, 
												 attribute.name="edgeType", 
												 attribute.type="char", 
												 default.value="undefined")#Create the edge attribute slot that the Bioconductor graph class requires
	cat(".")
	g <- initEdgeAttribute(g, attribute.name="weight", 
												 attribute.type="numeric", 
												 default.value=1)
	cat(".")
	cy <- CytoscapeConnection()
	if (pathway@title %in% as.character(getWindowList(cy)))
		deleteWindow(cy, pathway@title)
	cat(".")
	w <- new.CytoscapeWindow(pathway@title, g)
	cat(".")
	print("Sending path to graphite.. ")
	
	try(displayGraph(w), silent=T)
	#print("About to re-display the graph.  .. ")
	#displayGraph(w)#puts all the nodes on the graph, but they are on top of eachother
	print("displayed the graph")
	node.colors <- abs_color_vector
	# 	data.values = noa(getGraph(w),"label")
	#make sure the node.colors are in the same order as the data.values
	# 	node.colors=node.colors[data.values]
	
	print("setting node color rule")
	
	if(pathName%in%ab_analysis$path_id){
		#if it's targeted, then pull out the targets
		targs = pnodes%in%rownames(ab_analysis$coverage_summary$genesummary)
	}
	
	w = setNodeColors_graphite(w=w,
										targeted=rownames(ab_analysis$coverage_summary$genesummary), 
										active=rownames(ab_analysis$genesummary), 
										default.color="#ffffff",
										node_color_choices=c("#a0a0a0","#ffffff","#fcff00"))
	
	setDefaultNodeSize(obj=w, new.size=40)
	w=setNodeSize(w=w,abAnalysis=ab_analysis, defaultSize=40,maxSize=100)
	
	print("setting edge target rule")
	tres = try(silent=T,expr={
		setEdgeTargetArrowRule(obj=w, 
													 edge.attribute.name="edgeType", 
													 attribute.values=c(edgeTypes, "multiple"), 
													 arrows=c(edgeArrows, "No Arrow"))
	})
	if(length(grep("error", tres, ignore.case=T))){
		print(tres)
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		print("NOTICE: edge attributes for path could not be set")
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	}
	
	targs = rep(F, times=length(pnodes))
	### check drug screen targeting
	# first, check if pathway is even targeted:
	dsPaths = drug_analysis$coverage_summary$pathsummary
	if(pathName%in%dsPaths$path_id){
		#if it's targeted, then pull out the targets
		targs = pnodes%in%rownames(drug_analysis$coverage_summary$genesummary)
	}
	#targNodes = targs[]
	
	targetedSymbols=rownames(drug_analysis$coverage_summary$genesummary)
	activeSymbols = rownames(drug_analysis$genesummary)
	#make sure active doesn't get written over by targeted
	targetedSymbols = setdiff(targetedSymbols, activeSymbols)
	#"#f90000" = RED ff3d3d=lighter red
	#"#0400f9" = BLUE 4d94ff=lighter blue
	#"#686868" = grey
	#"#a633ff" = light purple
	if(!is.null(drug_analysis)){
		w = setBorderColors(w=w,
												color_choices=c("#a0a0a0" ,"#4d94ff","#ff3d3d","#a633ff"),
												targeted=targetedSymbols, 
												active=activeSymbols)
	}
	
	setDefaultNodeBorderWidth(w, new.width=8)
	
	setDefaultNodeFontSize(obj=w, new.size=8)
	setDefaultBackgroundColor(w, '#ffffff')
	print("'bout to refresh layout..")
	layoutNetwork(w, "force-directed")
	print("bout to redraw layout. . ")
	redraw(w)
	print("layout redrawn.")
	return(w)
}#paintPath()


# setNodeColors(w=w,
# 							targeted=rownames(ab_analysis$coverage_summary$genesummary), 
# 							active=rownames(ab_analysis$genesummary), 
# 							default.color="ffffff",
# 							color_choices=c("#a0a0a0","ffffff","#fbff00"), 
# )
setNodeColors_graphite<-function(w, node_color_choices, targeted, active=NULL, default.color="#ffffff"){
	
	graph_nodes=noa(getGraph(w),"label")
	
	#set the default, known non-aberrational color
	colorVector = rep(default.color, times=length(graph_nodes))
	
	#set the aberrational color
	colorVector[graph_nodes%in%active] = node_color_choices[3]
	
	#if there are targeted genes, make the not targeted come up grey
	if(!is.null(targeted)){
		#readline("a targeted set was passed")
		colorVector[!graph_nodes%in%targeted] = node_color_choices[1]
	} 
	
	setNodeColorRule(w, 
									 node.attribute.name="label", 
									 control.points=noa(getGraph(w),"label"), 
									 colors=colorVector, 
									 mode="lookup", 
									 default.color=default.color)
	
	# 	setNodeBorderColorRule(obj=w,
	# 												 mode="lookup",
	# 												 colors=colorVector,
	# 												 control.points=noa(getGraph(w),"label"),
	# 												 node.attribute.name="label",
	# 												 default.color = default.color
	# 	)
	return(w)
}

setNodeSize<-function(w, abAnalysis, defaultSize=40, maxSize=100){
	print("Setting node sizes")
	# 	abAnalysis = results$sequence_capture_whole_exome_aberration_summary
	
	#get nodes available to examine
	graph_nodes=noa(getGraph(w),"label")
	
	#determine range of values
	gsum = abAnalysis$genesummary
	
	#subset the genes in the graph
	#get genes in graph and in aberration set
	abInPath  = intersect(graph_nodes, rownames(gsum))
	selgsums = gsum[abInPath,,drop=F]
	
	sizes = rep(x=defaultSize, times=length(graph_nodes))
	names(sizes) = graph_nodes
	#determine range for nodes in graph
	#min will be 1, max will be the max for the graph
	maxCount = range(selgsums)[2]
	minCount = range(selgsums)[1]
	sRange = maxSize-defaultSize
	if(maxCount!=minCount){
		geneFacts = (selgsums-minCount)/(maxCount - minCount)
		#this is how much of the max size a gene will be
		#if the gene is the largest, this will be 1
		#So, this should be multiplied by sRange to get the size 
		#addition; the resultant should then be added to the default
		#size to get the node size
		sizeSet =  (geneFacts*sRange) + defaultSize
		print(sizeSet)
		sizes[rownames(selgsums)] = sizeSet
	}	
	setNodeSizeRule(obj=w,
									node.sizes=sizes, 
									mode="lookup",
									control.points=noa(getGraph(w),"label"),
									node.attribute.name="label")
	
	return(w)
}

#'@title setNodeColors()
#'@param w the cytoscape connection
#'@param study A Study object. Used to access Path_Detail object, and determine which genes are in which pathways
#'@param node_color_choices a vector, length 3, containing the hex color codes for the targeted, the unactive and the active nodes
#'@param targeted the vector of nodes considered to be targeted
#'@param active the vector of nodes considered to be active
#'@param default.color the hex color used as the default color 
#'@export
setNodeColors<-function(w, study, node_color_choices, targeted, active=NULL, default.color=NULL){
	
	graph_nodes=noa(getGraph(w),"label")
	nodesInActive = intersect(graph_nodes, active)
  nodesInTargeted = intersect(graph_nodes, targeted)
	
	activeColorVector = rep(node_color_choices[3], times=length(nodesInActive))
	
	#set the targeted color
	targetedColorVector = rep(node_color_choices[1], times = length(nodesInTargeted))	
	fullColorVector = c(activeColorVector, targetedColorVector)
	fullNodeSet = c(nodesInActive, nodesInTargeted)
	
	#get the set of genes in the pathway
	pgenes = getGenesFromPaths(pids=w@title, STUDY=study)
	#find what genes are left
	rnodes = setdiff(pgenes, fullNodeSet)
	#if there are targed nodes, the rnodes will be nodes we know nothing about
	#if there are not targeted nodes, it is assumed we know the targeted nodes to be in their 'normal' state
	
	if(!length(targetedColorVector)){
		#if there are not targeted nodes provided, the coverage of the analysis is not reduced
		#and so all rnodes are assumed to be in their normal state
		fullColorVector=c(fullColorVector, rep(node_color_choices[2], times=length(rnodes)))
		fullNodeSet = c(fullNodeSet, rnodes)
	}#if(condition is false, the rnodes will persist in their default color)

	setNodeColorRule(w, default.color=default.color,
									 node.attribute.name="label", 
									 control.points=fullNodeSet, 
									 colors=fullColorVector, 
									 mode="lookup")
	
	return(w)
}

#"#f90000" = RED ff3d3d=lighter red
#"#0400f9" = BLUE 4d94ff=lighter blue
#"#686868" = grey
#"#a633ff" = light purple
# w = setBorderColors(w=w,
# 										color_choices=c("#a0a0a0" ,"#4d94ff","#ff3d3d","#a633ff"),
# 										targeted=targetedSymbols, 
# 										active=activeSymbols)
#'@title setBorderColors
#'@param w cytoscape connection object, 
#'@param color_choices vector of color 4 choices, the default, the targeted & not active nodes, the targeted and active nodes and the color to use if there is only coverage data
#'@param targeted vector of targeted genes (all strings)
#'@param active vector of active genes (all strings)
#'@param default.color the default
#'@export
setBorderColors<-function(w, 
													color_choices, 
													targeted=NULL, 
													active, 
													default.color="#a0a0a0"){
	
	graph_nodes=noa(getGraph(w),"label")
	
	#set the default color: grey (a0a0a0) ususally
	colorVector = rep(color_choices[1], times=length(graph_nodes))
	
	#if there is only coverage data (targeted not null, active null)
	if(!is.null(targeted)&is.null(active)){
		#set the active color to the last color choice
		colorVector[graph_nodes%in%targeted] = color_choices[4]
	}else if(!is.null(targeted)&!is.null(active)){
		#if active and targeted are there, set the active and the inactive colors
		colorVector[graph_nodes%in%targeted] = color_choices[2]
		colorVector[graph_nodes%in%active] = color_choices[3]
	}else{
		readline("setBorderColors() in an unusual state: the active genes are not null but the targeted is")
	}
	
	# 	#set the color for the nodes that are targeted
	# 	colorVector[graph_nodes%in%targeted] = color_choices[2]
	# 	
	# 	if(!is.null(targeted)){
	# 		readline("in set boarder colors, targeted is not null")
	# 		colorVector[graph_nodes%in%targeted] = color_choices[3]
	# 	}
	setNodeBorderColorRule(obj=w,
												 mode="lookup",
												 colors=colorVector,
												 control.points=noa(getGraph(w),"label"),
												 node.attribute.name="label",
												 default.color = default.color
	)
	return(w)
}

cytoscapePlot <- function(pathway, ...) {
	
	g <- buildGraphNEL(nodes(pathway), edges(pathway), FALSE)#, ...)
	g <- markMultiple(g)#marks edgeType attribute
	g <- initEdgeAttribute(graph=g, 
												 attribute.name="edgeType", 
												 attribute.type="char", 
												 default.value="undefined")#Create the edge attribute slot that the Bioconductor graph class requires
	g <- initEdgeAttribute(g, attribute.name="weight", 
												 attribute.type="numeric", 
												 default.value=1)
	g<-initNodeAttribute(g,
											 attribute.name="act", 
											 attribute.type="integer", 
											 default.value="off")
	nodeData(g,attr="act")<-c(1,1,1,0,0)
	
	#length(noa(getGraph(w),"label"))
	
	cy <- CytoscapeConnection()
	if (pathway@title %in% as.character(getWindowList(cy)))
		deleteWindow(cy, pathway@title)
	
	w <- new.CytoscapeWindow(pathway@title, g)
	displayGraph(w)#puts all the nodes on the graph, but they are on top of eachother
	node.colors <- c ("#0000AA", "#FFFF00")
	data.values = c(1,0)#noa(getGraph(w),"label")
	# 	setNodeColorRule(w, 
	# 									 node.attribute.name="act", 
	# 									 control.points=data.values, 
	# 									 colors=node.colors, 
	# 									 mode="interpolate", 
	# 									 default.color='#000000')
	setEdgeTargetArrowRule(w, 
												 edge.attribute.name="edgeType", 
												 attribute.values=c(edgeTypes, "multiple"), 
												 arrows=c(edgeArrows, "No Arrow"))
	layoutNetwork(w)
	redraw(w)
	return(w)
}#cytoscapePlot()

# Copyright 2011 Gabriele Sales <gabriele.sales@unipd.it>
#
#
# This file is part of graphite.
#
# graphite is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# graphite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with graphite. If not, see <http://www.gnu.org/licenses/>.


pathwayGraph <- function(pathway, edge.types=NULL) {
	buildGraphNEL(nodes(pathway), edges(pathway), TRUE, edge.types)
}

buildGraphNEL <- function(nodes, edges, sym, edge.types=NULL) {
	cat("A")
	if (!is.null(edge.types))
		edges <- selectEdges(edges, edge.types)
	cat("B")
	if (NROW(edges) == 0){
		cat("C1")
		g <- new("graphNEL", nodes, list(), "directed")
		cat("C11")
	} else {
		cat("C2")
		edges <- prepareEdges(as.matrix(edges), sym)
		print(dim(edges))
		cat("C201")
		# 		print(edges)
		# 		print(nodes)
		newedges = edgeList(nodes, edges)
		cat(" Edges extracted. ")
		g <- new("graphNEL", nodes=nodes, edgeL=newedges , edgemode="directed")
		cat("C202")
		edgeDataDefaults(g, "edgeType") <- "undefined"
		cat("C203")
		edgeData(g, edges[,1], edges[,2], "edgeType") <- edges[,3]
		cat("C22")
	}
	
	return(g)
}

selectEdges <- function(m, types) {
	unknownTypes <- setdiff(types, edgeTypes)
	if (length(unknownTypes))
		stop("the following edge types are invalid: ", paste(unknownTypes, collapse=", "))
	
	m[m[,4] %in% types,]
}

prepareEdges <- function(m, sym) {
	ns         <- canonicalEdgeNames(m)
	simplified <- matrix(unlist(tapply(1:NROW(m), ns, function(is) mergeEdges(m, is))),
											 ncol=4, byrow=T)
	
	if (sym)
		symmetricEdges(simplified)
	else
		simplified[, -3, drop=FALSE]
}

canonicalEdgeNames <- function(m) {
	apply(m, 1, function(e) {
		if (e[1] <= e[2])
			paste(e[1], e[2], sep="|")
		else
			paste(e[2], e[1], sep="|")
	})
}

mergeEdges <- function(m, is) {
	h <- m[is[1],]
	if (length(is) == 1)
		h
	else {
		if ("undirected" %in% m[is,3] || any(h[1]!=m[is,1]))
			dir <- "undirected"
		else
			dir <- "directed"
		
		c(h[1], h[2], dir, paste(unique(m[is,4]), collapse=";"))
	}
}

symmetricEdges <- function(m) {
	undirected <- m[m[,3]=="undirected" & m[,1]!=m[,2], c(2,1,4), drop=FALSE]
	
	if (NROW(undirected) > 0) {
		full <- m[, -3, drop=FALSE]
		stopifnot(is.null(dimnames(full)))
		rbind(full, undirected)
	} else
		return(m[, -3, drop=FALSE])
}

edgeList <- function(nodes, edges) {
	
	sapply(nodes,
				 function(n) list(edges=edges[edges[,1]==n, 2]),
				 simplify=FALSE,
				 USE.NAMES=TRUE)
	
}

markMultiple <- function(g) {
	d <- edgeData(g)
	if (length(d) == 0)
		return(g)
	
	ns <- names(d)
	for (i in 1:length(d)) {
		tp <- d[[i]]$edgeType
		if (length(grep(";", tp, fixed=T)) > 0) {
			nodes <- unlist(strsplit(ns[[i]], "|", fixed=T))
			edgeData(g, nodes[1], nodes[2], "edgeType") <- "multiple"
		}
	}
	
	return(g)
}


biopaxFileFromPathName<-function(){
	
	readline("function not yet implemented: biopaxFileFromPathName()")
	#check if we got it
	#if we dont have it get it
}

# Copyright 2011 Gabriele Sales <gabriele.sales@unipd.it>
#
#
# This file is part of graphite.
#
# graphite is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# graphite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with graphite. If not, see <http://www.gnu.org/licenses/>.


edgeTypes <- c("binding","catalysisIn(ACTIVATION)","catalysisOut(ACTIVATION)","catalysisOut(INHIBITION)","indirect","biochemicalReaction","complexAssembly","catalysisIn(INHIBITION)","binding/association","methylation","activation","expression","inhibition","phosphorylation","dephosphorylation","indirect effect","dissociation","ubiquination","repression","missing","state change","missing interaction","catalysisOut(INHIBITION-COMPETITIVE)","modulation(ACTIVATION)")

edgeArrows <- c("No Arrow", "Arrow","Arrow","T","No Arrow","Arrow","No Arrow","T","No Arrow","No Arrow","Arrow","Arrow","T","No Arrow","No Arrow","No Arrow","No Arrow","No Arrow","T","No Arrow","No Arrow","No Arrow","T","Arrow")

spiaAttributes <- c("activation", "compound", "binding/association", "expression", "inhibition", "activation_phosphorylation", "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation", "dissociation", "dephosphorylation", "activation_dephosphorylation", "state change", "activation_indirect effect", "inhibition_ubiquination", "ubiquination", "expression_indirect effect", "inhibition_indirect effect", "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation", "activation_binding/association", "indirect effect", "activation_compound", "activation_ubiquination")

spiaConv <- matrix(ncol=2, byrow=TRUE,
									 c("binding", "binding/association",
									 	"control(In(ACTIVATION))", "activation",
									 	"control(In(INHIBITION))", "inhibition",
									 	"control(Out(ACTIVATION))", "activation",
									 	"control(Out(INHIBITION))", "inhibition",
									 	"control(Out(INHIBITION-COMPETITIVE))", "inhibition",
									 	"control(Out(ACTIVATION_UNKMECH))", "activation",
									 	"control(Out(unknown))", "indirect effect",
									 	"control(indirect)", "indirect effect",
									 	"process", "activation",
									 	"process(BiochemicalReaction)", "activation",
									 	"process(activation)", "activation",
									 	"process(binding/association)", "binding/association",
									 	"process(dephosphorylation)", "dephosphorylation",
									 	"process(dissociation)", "dissociation",
									 	"process(expression)", "expression",
									 	"process(indirect effect)", "indirect effect",
									 	"process(indirect)", "indirect effect",
									 	"process(inhibition)", "inhibition",
									 	"process(missing interaction)", "indirect effect",
									 	"process(missing)", "indirect effect",
									 	"process(phosphorylation)", "phosphorylation",
									 	"process(repression)", "inhibition",
									 	"process(ubiquitination)", "ubiquination"))
colnames(spiaConv) <- c("type", "spiaType")
