

# test.addDrugsToGraph<-function(){
# 	
# 	splist = STUDY@results$overlap_analysis$"Aberration enriched, containing sensitive targets"
# 	dplist = STUDY@results$overlap_analysis$"Aberration enriched, not drug targeted"
# 	testPath1 = dplist$path_id[3]
# 	
# 	study=STUDY
# 	pname=testPath1
# 	
# 	pname = splist[1,1]
# 	
# }


addDrugsToGraph<-function(w, study, pname, drugID){
	
	###overall goal: for the path in question, show all the drugs targeting 
	###the genes in that path, and connections to the genes they target
	###
	
	####get drugs targeting the pathways 
	
	#path -> gene list -> drugSet
	geneList = getGenesFromPaths(pids=pname, STUDY=study)
	dtd = getAllDrugGeneAssociations(study=study)
	
	sif = dtd[dtd$geneID%in%geneList,]
	
	#drugSet -> targets in path -> sif
	colnames(sif)<-c("gene","drug")
	
	#add drug names 
	dmd = getDrugData()
	dmd = dmd[,c("DrugBank.ID","Name")]
	sif = merge(x=sif, y=dmd, by.x="drug", by.y="DrugBank.ID", all.x=T)
	sif$Name[is.na(sif$Name)] = sif$drug[is.na(sif$Name)]
	
	#add drug as nodes 	#add sif
	w = placeDrugs(sif=sif, w=w)

	#re-layout
	setLayout(tw=w, placeInOrganelles=F)
}

getAllDrugGeneAssociations<-function(study,
																		 fnames=c("./reference_data/drugDB/drugbank/all_target_ids_all.csv",
																		 				 "./reference_data/drugDB/drugbank/all_enzyme_ids_all.csv")){
	dtd2 = NULL
	for(fn in fnames){
		dtd1 = getDrugTargetData(hugo=study@studyMetaData@paths$symtable, 
														 fname=fn)
		dtd2 = rbind(dtd2,dtd1)
	}
	
	return(unique(dtd1))
}

placeDrugs<-function(sif, w, drugColor="#ffa500"){
	
	gtmp = w@graph 
	nlab = noa(graph=gtmp, node.attribute.name="label")
	
	gtmp = addNode(node=unique(sif$drug), object=gtmp)
	
	nodeData(gtmp, n=sif$drug, attr="label")<-sif$Name
	
	csn = c()
	ctn = c()
	for(i in 1:nrow(sif)){#for each drug-gene association
		cGene = sif$gene[i]
		cDrug = sif$drug[i]
		
		gi=which(nlab%in%cGene)
		targets = names(nlab)[gi]
		
		csn = c(csn,
						rep(cDrug, times=length(targets)))
		ctn = c(ctn, 
						targets)
	}
	
	for(e in 1:length(csn)){#for each edge
		addCyEdge(obj=w, 
							sourceNode=csn[e], 
							targetNode=ctn[e], 
							edgeType="drugTarget", 
							directed=T)
	}
	
	eNames = paste(csn, ctn, sep=" (drugTarget) ")
	
	udrug = sif[!duplicated(sif$drug),]
	
	setNodeColorDirect(obj=w, node.names=udrug$drug, new.color=drugColor)
	setNodeFontSizeDirect(obj=w, node.names=udrug$drug, new.sizes=24)
	setNodeLabelDirect(obj=w, node.names=udrug$drug, new.labels=udrug$Name)
	
	setEdgeColorDirect(obj=w, edge.names=eNames, new.value=drugColor)
	setEdgeLineStyleDirect(obj=w, edge.names=eNames, new.values="VERTICAL_SLASH")
	setEdgeSourceArrowShapeDirect(obj=w, edge.names=eNames, new.values="No Arrow")
	setEdgeTargetArrowColorDirect(obj=w, edge.names=eNames, new.colors=drugColor)
	setEdgeLineWidthDirect(obj=w, edge.names=eNames, new.value=5)
	

	
	redraw(w)

	return(w)
}

stylizeDrugs<-function(tw,dnames){
	
	nodeTable = inferPositions(w=tw, nodeTable=nodeTable)
	opos = getSpaceDividers(nodeTable=nodeTable,w=tw)
	
	lockNodeDimensions(obj=tw, new.state=FALSE)
	
	setNodeWidthDirect(obj=tw, node.names=opos$names, new.widths=as.numeric(opos$width))
	setNodeHeightDirect(obj=tw, node.names=opos$names, new.heights=opos$height)
	
	setNodeShapeDirect(obj=tw, node.names=opos$names, new.shapes="round_rect")
	
	tmp = noa(graph=tw@graph, node.attribute.name="nodeType")
	tmp[tmp==""] = "cellularCompartment"
	# noa(graph=tw@graph, node.attribute.name="nodeType")<-tmp
	setNodeAttributesDirect(obj=tw, 
													attribute.name="nodeType", 
													attribute.type='char', 
													node.names=names(tmp), 
													values=tmp)
	
	setNodeOpacityDirect(obj=tw, node.names=opos$names, new.values=0)
	setNodeBorderOpacityDirect(obj=tw, node.names=opos$names, new.values=255)
	setNodeBorderWidthDirect(obj=tw, node.names=opos$names, new.sizes=5)
	
	if('cytosol'%in%opos$names){
		#make the cytosol box not have a border
		setNodeBorderOpacityDirect(obj=tw, node.names='cytosol', new.values=0)
	}
	redraw(tw)
	
	setNodeLabelDirect(obj=tw, node.names=opos$names, new.labels=opos$names)
	setNodeLabelColorDirect(obj=tw, node.names=opos$names, new.colors=rep("#000000", times=length(opos$names)))
	setNodeFontSizeDirect(obj=tw, node.names=opos$names, new.sizes=rep(30, times=length(opos$names)))
	setNodeLabelOpacityDirect(obj=tw, node.names=opos$names, new.values=rep(100, times=length(opos$names)))
	
	# 	setNodeAttributesDirect(obj=tw, values=c(noa(graph=tw@graph, node.attribute.name="nodeType"), opos$names), 
	# 													attribute.name="label", 
	# 													attribute.type="char", 
	# 													node.names=names(noa(graph=tw@graph, node.attribute.name="label")))
	# 	setNodeBorderColorDirect(obj=tw, node.names=opos$names, new.color="#000000")
	setNodePosition(obj=tw, node.names=opos$names, x.coords=opos$centerX, y.coords=opos$centerY)
	redraw(tw)
	
}

	
getRecordsWithNode<-function(biopax){

	bpdf = biopax$df
	bpdf = unfactorize(bpdf)
	
	searchTerm = "(Complex2$)|(Complex2-)"
	
	classi = grepl(pattern=searchTerm, x=bpdf$class)
	idi = grepl(pattern=searchTerm, x=bpdf$id)
	propertyi = grepl(pattern=searchTerm, x=bpdf$property)
	propatti = grepl(pattern=searchTerm, x=bpdf$property_attr)
	propattvali = grepl(pattern=searchTerm, x=bpdf$property_attr_value)
	propvali = grepl(pattern=searchTerm, x=bpdf$property_value)
	res = bpdf[(classi|idi|propertyi|propatti|propattvali|propvali),]
	dim(res)
	
}

graphNELToCytoscape<-function(gnel, pathname="no path name provided as argument to graphNELToCytoscape()") 
{
	g = gnel
	g <- markMultiple(g)
	g <- initEdgeAttribute(g, "edgeType", "char", "undefined")
	g <- initEdgeAttribute(g, "weight", "numeric", 1)
	cy <- CytoscapeConnection()
	if (pathname %in% as.character(getWindowList(cy))) 
		deleteWindow(cy, pathname)
	w <- new.CytoscapeWindow(pathname, g)
	displayGraph(w)
	setEdgeTargetArrowRule(w, "edgeType", c(edgeTypes, "multiple"), 
												 c(edgeArrows, "No Arrow"))
	layoutNetwork(w)
	redraw(w)
	return(w)
}

test.setAberrationDataStyles<-function(){
	study=STUDY
	resSet = STUDY@results$functional_drug_screen_summary
	pname ="Abacavir metabolism"
}

genesNotInPath<-function(pname, study, w){
	#get the nodes not in the path
	#get the nodes in the path
	geneNames = getGenesFromPaths(pids=pname, STUDY=study)
	nodesInGraph = noa(graph=w@graph, node.attribute.name="label")
	notInPathProts = setdiff(x=nodesInGraph, y=geneNames)
	return(notInPathProts)
}

addToStyle<-function(style, geneNames, value){
	
	style$names = c(style$names, geneNames)
	if(length(value)==1){style$colors = c(style$colors, rep(value, times=length(geneNames)))
	}else if(length(value)==length(geneNames)){
		style$colors = c(style$colors, value)
	}else{
		print("*******************error: length of gene vector not equal to length of visual quality vector!!*******")
	}
	
	return(style)
}

setAberrationDataStyles<-function(pname, resSetNombre,
																	study, w, resSet, 
																	normalProteinColor="#4aa7ff",
																	notInPathColor="#F0F8FF",
																	abColor = "#fcff00", 
																	defaultColor="#686868"){ #default color will be equal to the "dark" color
	print("Setting styles for aberration data")
	style = list()
	style[["names"]] = c()#for readability
	style[["colors"]] = c()#for readability
	
# 	resSet = study@results[[resSetNombre]]
	nodesInGraph = noa(graph=w@graph, node.attribute.name="label")
	#### Assign the not-in-path color
	nipprots = genesNotInPath(pname=pname, study=study, w=w)
	names(nipprots)<-NULL
	#add to the color vector
	style = addToStyle(style=style, geneNames=nipprots, value=notInPathColor)

	#### Assign the aberration colors
	#get the aberrational genes
	abnodes = resSet$genesummary
	#pick out the aberrational genes in the graph
	abInGraph = rownames(abnodes)[rownames(abnodes)%in%nodesInGraph]
	style=addToStyle(style=style, geneNames=abInGraph, value=abColor)
	
	#### check if there is a coverage issue to handle
	if(!is.null(resSet$coverage_summary$genesummary)){
		readline("unit test the coverage coloration")
		#if there is a coverage issue, find the set difference between the covered
		#nodes and the nodes already in the nameVector (abInGraph?) these will get the
		#"normalProteinColor"
		coveredGenes = resSet$coverage_summary$genesummary
		coveredNormal = setdiff(coveredGenes, abInGraph)
		style=addToStyle(style=style, geneNames=coveredNormal, value=normalProteinColor)
	}else{#if there's no coverage issue, assign all remaining as normal
		#### Assign the color for the normal, non-ab proteins
		#get the names of the proteins in the path
		pgenes = intersect(getGenesFromPaths(pname, STUDY=study), nodesInGraph)
		normalGenes = setdiff(x=pgenes, y=style$names)
		style=addToStyle(style=style, geneNames=normalGenes, value=normalProteinColor)
	}
	
	#now send the colors to the graph
	setNodeColorRule(default.color=defaultColor, 
									 mode='lookup',
									 colors=style$colors,
									 control.points=style$names,
									 node.attribute.name="label",
									 obj=w)
}

setFunctionalDataStyle<-function(pname, study, w, resSet, resSetNombre,
																 default.color="#686868",
																 notCovered.color="#908b8b",
																 insensitive.color="#9cd3ff", 
																 sensitive.color="#ff6060", 
																 nonMonitoredBorder.width=2, 
																 monitoredBorder.width=15, 
																 default.width=1){
	
	style = list()
	style[["names"]] = c()#for readability
	style[["colors"]] = c()#for readability
	
	nodesInGraph = noa(graph=w@graph, node.attribute.name="label")

	nodesInPath = getGenesFromPaths(pids=pname, STUDY=study)
	nodesInPatNotInGraph = setdiff(nodesInPath, nodesInGraph)
	if(length(nodesInPatNotInGraph)){

		cat("\n\nWarning!!!!, not all of the expected genes were found in the current network\n\n")
		cat("The missing nodes is/are:", nodesInPatNotInGraph)
		nodesInPath = intersect(x=nodesInPath, y=nodesInGraph)
		
	}
	
	tmp = nodesInGraph
	names(tmp)=NULL
	nonPathNodes = setdiff(tmp, nodesInPath)

	coveredGenes = intersect(rownames(resSet$coverage_summary$genesummary), nodesInGraph)
													 
	notCoveredGenes = setdiff(nodesInPath, coveredGenes)
	
	sensitiveGenes = intersect(rownames(resSet$genesummary), nodesInGraph)
	
	insensitive = setdiff(coveredGenes, sensitiveGenes)
	
	
	
	########### set the colors for the active/sensitive
	style = addToStyle(style=style, geneNames=sensitiveGenes, value=sensitive.color)
	############ set the colors for the inacive/insensitive
	style = addToStyle(style=style, geneNames=insensitive, value=insensitive.color)
	############ set the not covered in analysis color
	style = addToStyle(style=style, geneNames=notCoveredGenes, value=notCovered.color)

	remainingNodes = setdiff(nodesInGraph, style$names)
	
	#need to convert node labels to node names
	
	style=addToStyle(style=style, geneNames=remainingNodes, value=default.color)
	
	colorVector = valuesForAllNodes(style=style, 
																	nodesInGraph=nodesInGraph, 
																	default=default.color)
	for(i in 1:length(names(nodesInGraph))){
		setNodeBorderColorDirect(obj=w,
														 node.names=names(nodesInGraph)[i], 
														 new.color=colorVector[i])
	}

	wstyle = list()
	wstyle[["names"]] = c()#for readability
	wstyle[["colors"]] = c()#for readability
	
	#get the nodes not in the gene set, and set their border width to 1
	wstyle = addToStyle(style=wstyle, geneNames=nonPathNodes, value=nonMonitoredBorder.width)
	#set the border widths of the nodes which are monitored
	wstyle = addToStyle(style=wstyle, geneNames=nodesInPath, value=monitoredBorder.width)
	
	wvalues = valuesForAllNodes(style=wstyle, nodesInGraph=nodesInGraph, default=1)

	setNodeBorderWidthDirect(obj=w, node.names=names(nodesInGraph), new.sizes=wvalues)

	setNodeBorderWidthRule(obj=w,
												 node.attribute.name="label",
												 attribute.values=wstyle$names, 
												 line.widths=wstyle$colors,
												 default.width=default.width)

	redraw(w)
	
}

valuesForAllNodes <- function (style, nodesInGraph, default) {
	styleDict = style$colors
	names(styleDict) = style$names
	colorSet = styleDict[nodesInGraph]
	colorSet[is.na(colorSet)] =  default
	names(colorSet) = NULL
	return(colorSet)
}


setBiologicalDataStyles<-function(pname, 
																	study, 
																	w, 
																	resSetName){
	
	print("Setting schema for display of experimental data...")
	#extract the correct part of the study
	studsub = study@results[[resSetName]]

	if(grepl(x=resSetName, pattern="aberration", ignore.case=T)){
		print("found aberration data")
		setAberrationDataStyles(pname=pname, 
														study=study, w=w, 
														resSet = studsub, 
														resSetNombre = resSetName)
			
	}else if(grepl(x=resSetName, pattern="function", ignore.case=T)){

		print("found functional data")
		setFunctionalDataStyle(pname=pname, monitoredBorder.width=20,
													 study=study, 
													 w=w, 
													 resSet = studsub) #, resSetName = resSetName)
	

	}else if(grepl(x=resSetName, pattern="overlap", ignore.case=T)){
		
		print("found overlap data")
		
		setAberrationDataStyles(pname=pname, 
														study=study, 
														w=w, 
														resSet = studsub$combined_aberrations_summary, 
														resSetNombre = "combined_aberrations_summary")
		
		setFunctionalDataStyle(pname=pname, monitoredBorder.width=20,
													 study=study, 
													 w=w, 
													 resSet = studsub$functional_drug_screen_summary, 
													 resSetNombre = "functional_drug_screen_summary")
		
	}

	print("experimental data color scheme established")
	
}

setDefaultCytoscapeStyle<-function(pname, study, w){
	
	print("setting default Cytoscape network display style...")
	default.size = 30
	
	# 	pname="Abacavir metabolism/test 2"
	#get node rows from the biopax object
	
	getVisualStyleNames(obj=w)
# 	layoutNetwork(obj=w, layout.name="force-directed")
	layoutNetwork(obj=w, layout.name="kamada-kawai")
#  	setVisualStyle(w, new.style.name="Solid")
	
	lockNodeDimensions(obj=w, new.state=FALSE)
	
	gsymbols = getGenesFromPaths(pids=pname, STUDY=study)

	print(pname)
	
	mygraph = getAllNodes(w)
	cgraph = getGraph(obj=w)
	nodeDict = noa(graph=cgraph, node.attribute.name="nodeType")
	nodelables = noa(graph=cgraph, node.attribute.name="label")

	labelsOfProteins = nodelables[nodeDict=="Protein"]
	
	setDefaultBackgroundColor(w, new.color="#FFFFFF")
	setDefaultNodeShape(w, new.shape="ellipse")
	setDefaultNodeSelectionColor(obj=w, new.color="#DFFFA5")
	
	#adjust the graph as to the node types:

	controlSubclasses= c("Catalysis", 
											 "Modulation", 
											 "TemplateReactionRegulation")
	conversionSubclasses = c("BiochemicalReaction", 
													 "ComplexAssembly", 
													 "Degradation",
													 "Transport",
													 "TransportWithBiochemicalReaction")
	interactionClasses = c("TemplateReaction", 
												 "MolecularInteraction", 
												 "GeneticInteraction", 
												 "Conversion", 
												 "Control")
	
	physicalEntities=c("Complex", 
										 "DNARegion",
										 "DNA", 
										 "Protein",
										 "RNA",
										 "RNARegion",
										 "SmallMolecule")

	allClasses = c(interactionClasses,conversionSubclasses,controlSubclasses,physicalEntities)
	
	shapeDict = rep("ellipse",times=length(allClasses))
	names(shapeDict)<-allClasses
	
	shapeDict[conversionSubclasses] = "triangle"
	shapeDict[interactionClasses] = "trapezoid"
	shapeDict[controlSubclasses] = "diamond"
	shapeDict[physicalEntities]  = c("hexagon", 
																	 "rect", 
																	 "rect", 
																	 "round_rect", 
																	 "rect", 
																	 "rect", 
																	 "ellipse")
	shapeDict["PhysicalEntity"]="hexagon"
	shapeDict["CellularCompartment"] = "round_rect"
	
	nn = names(shapeDict)
	nv = shapeDict
	names(nv)<-NULL
	setNodeShapeRule(obj=w,
									 attribute.values=nn,
									 node.shapes=nv,
									 node.attribute.name="nodeType")
	
	sizeDict = rep(default.size,times=length(allClasses))
	names(sizeDict)<-allClasses
	
	sizeDict[conversionSubclasses] = 25
	sizeDict[interactionClasses] = 15
	sizeDict[controlSubclasses] = 15
	sizeDict[physicalEntities]  = c(30, 
																	30, 
																	40, 
																	50, 
																	30, 
																	30, 
																	30)
	sizeDict["PhysicalEntity"]=40
	sizeDict[c(conversionSubclasses, interactionClasses, controlSubclasses)] = 20
	
	nn = names(nodeDict)
	nv = sizeDict[nodeDict]
	names(nv)<-NULL
	nv = as.numeric(nv)

	#this has got to be converted to node names because it's height and width direct
	
	setNodeWidthDirect(obj=w, node.names=nn, new.widths=nv)
	setNodeHeightDirect(obj=w, node.names=nn,new.heights=nv)

	###### set font sizes
	sizeTypeDict = rep(14, times=length(allClasses))
	names(sizeTypeDict) = allClasses
	sizeTypeDict[c(interactionClasses,conversionSubclasses,controlSubclasses)] = 10
# 	physicalEntities=c("Complex", 
# 										 "DNARegion",
# 										 "DNA", 
# 										 "Protein",
# 										 "RNA",
# 										 "RNARegion",
# 										 "SmallMolecule")
	sizeTypeDict[physicalEntities] = c(14,14,20,30,14,14,16)

	nodeLabelFontSizeDict = sizeTypeDict[nodeDict]
	names(nodeLabelFontSizeDict) = names(nodeDict)
	
	nodeLabelFontSizeDict[is.na(nodeLabelFontSizeDict)] = 14#incase I missed anything
	nn = names(nodeLabelFontSizeDict)
	names(nodeLabelFontSizeDict)<-NULL
	setNodeFontSizeDirect(obj=w, node.names=nn, new.sizes=nodeLabelFontSizeDict)
	
	###### set node opacities
	
	opacityTypeDict = rep(255, times=length(allClasses))
	names(opacityTypeDict) = allClasses
	opacityTypeDict[c(interactionClasses,conversionSubclasses,controlSubclasses)] = 100
	# 	physicalEntities=c("Complex", 
	# 										 "DNARegion",
	# 										 "DNA", 
	# 										 "Protein",
	# 										 "RNA",
	# 										 "RNARegion",
	# 										 "SmallMolecule")
	opacityTypeDict[physicalEntities] = c(150,
																				255,
																				255,
																				255,
																				255,
																				255,
																				200)
	
	nodeLabelOpacityDict = opacityTypeDict[nodeDict]
	names(nodeLabelOpacityDict) = names(nodeDict)
	
	nodeLabelOpacityDict[is.na(nodeLabelOpacityDict)] = 255#incase I missed anything
	nn = names(nodeLabelOpacityDict)
	names(nodeLabelOpacityDict)<-NULL

	setNodeLabelOpacityDirect(obj=w, node.names=nn, new.values=nodeLabelOpacityDict)
	setNodeOpacityDirect(obj=w, node.names=nn, new.values=nodeLabelOpacityDict)


	inOriginalModelDict = noa(graph=cgraph, node.attribute.name="inOriginalModelColumn")

	if(sum(!is.na(inOriginalModelDict))){
		inOriginalModelDict = inOriginalModelDict[grepl(pattern="^node in quesion", x=inOriginalModelDict, ignore.case=T)]
		if(length(inOriginalModelDict)){
			print(inOriginalModelDict)
			setNodeColorDirect(obj=w, node.names=names(inOriginalModelDict), new.color="#FF9933")
			setNodeSizeDirect(obj=w, node.names=names(inOriginalModelDict), new.sizes=70)
		}
	}

	setNodeTooltipRule(obj=w, node.attribute.name="Alternate.Names")

	###### Edges
	######
	######
	membershipLabels = c("component", 
											 "memberPhysicalEntity")
	catalysisLabels = c("controller", 
											"controlled")
	reactionLables = c("left", 
										 "right")
	
	edgeTypeLabels = c(membershipLabels, 
										 catalysisLabels, 
										 reactionLables)
	
	setDefaultEdgeColor(w, new.color="#C0C0C0")
	setDefaultEdgeLineWidth(obj=w, new.width=12)
	
	setEdgeLineWidthRule(obj=w, 
											 edge.attribute.name="edgeType",
											 default.width=12,
											 attribute.values=c(membershipLabels, catalysisLabels), 
											 line.widths=c(5,5,7,7))
	
	setEdgeLineStyleRule(obj=w, 
											 edge.attribute.name="edgeType", 
											 default.style="SOLID",
											 line.styles=c("DOT","DOT"),
											 attribute.values=catalysisLabels)
	
	setEdgeOpacityRule(obj=w, 
										 edge.attribute.name="edgeType", 
										 control.points=c("component","memberPhysicalEntity", "left","right", "controller","controlled"), 
										 opacities=c(255,255,150,150,150,150), 
										 mode="lookup")


	setEdgeColorRule(obj=w, mode='lookup', default.color="#C0C0C0",
									 edge.attribute.name="edgeType", 
									 control.points=c("component"), 
									 colors=c("#C0C0C0"))
	
	targetShapes = c("Circle",
									 "Circle", 
									 "No Arrow", 
									 "Delta",
									 "No Arrow",
									 "Delta")
	
	sourceShapes = c("Circle",
									 "Circle",
									 "Delta",
									 "No Arrow",
									 "Delta",									 
									 "No Arrow")
	
	setEdgeTargetArrowRule(obj=w, attribute.values=edgeTypeLabels,
												 arrows=targetShapes,
												 default="No Arrow",
												 edge.attribute.name="edgeType")

	setEdgeSourceArrowRule(obj=w, default="No Arrow",
												 arrows=sourceShapes,
												 attribute.values=edgeTypeLabels,
												 edge.attribute.name="edgeType")
	
	setEdgeTargetArrowColorRule(obj=w, 
															colors="#C0C0C0",
															default.color="#C0C0C0",
															edge.attribute.name="edgeType", 
															attribute.values="component")
	
	setEdgeSourceArrowColorRule(obj=w, 
															colors="#C0C0C0",
															default.color="#C0C0C0",
															edge.attribute.name="edgeType", 
															attribute.values="component")
	
	redraw(w)
}

#'@title nodeNamesFromLabels
#'@param w a cytoscape connection object
#'@param nodeLabels a vector of node labels
#'@return a vector of node names corresponding to the supplied node labels
nodeNamesFromLabels<-function(w, nodeLabels){
	dict = noa(graph=w@graph, node.attribute.name="label")
	retvals = names(dict)[dict%in%nodeLabels]
	return(retvals)
}

#'@title getPathwaysRecords
#'@description Opens or creates record of downlaoded biopax pathways
#'@param pwrecord.fileName the file name for the pathway records file. 
#'@return data frame with columns: dbID, path_name and download_date
getPathwaysRecords<-function(pwrecord.fileName=NULL){
	print("getting biopax pathway records..")
	if(!file.exists(pwrecord.fileName)){ 
		cat("\nPathway record file not found, copying default.\n")
		checkFileCopyDefault(fname=pwrecord.fileName)
		# 		checkFileCopyDefault(fname="./reference_data/paths/biopax/2161541.owl")
		
		# 		pwrecord=data.frame(matrix(nrow=0,ncol=3,dimnames=c(list(NULL,c("dbID", "path_name", "download_date")))), 
		# 												stringsAsFactors=F)
		
	}
	pwrecord = read.table(file=pwrecord.fileName, 
											strip.white=T,
											quote="",
											comment.char="",
											stringsAsFactors=F,
											header=T,sep="\t")
	pwrecord2=pwrecord
	return(pwrecord2)
	
}



#'@title Download biopax files for reactome pathways. 
#'@description Attempts to download biopax files for reactome pathways. First, pathway's individual database identifiers are obtained from bioMART, then restful calls are made to Reactome's restful interface. A record of successfully downloaded pathways is made in the file "./reference_data/paths/biopax/record_of_biopax_pathways.txt" so that pathways will not be repeatedly downloaded. 
#'@param study A \code{Study} object. 
#'@param pathNames The names of the pathways the program should attempt to download. 
#'@param verbose \code{logical} flag indicating if user should be prompted before attempted download of biopax pathways. 
#'@return Character vector containing the names of any pathways that could not be downloaded.
#'@export
#'@import RCurl
#'@import biomaRt
getReactomeBiopax<-function(study, pathNames, verbose=T){
	
	biopax.dir = "./reference_data/paths/biopax/"
	if(!file.exists(biopax.dir)) dir.create(path=biopax.dir, recursive=T, showWarnings=F)
	
	# find which pathways are needed
	pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	pwrecord = getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	neededPaths = pathNames[!pathNames%in%pwrecord$path_name]
	
	if(length(neededPaths)){
		if(verbose){
			if("d"!=readline("Some biopax pathway records cannot be found on this computer.\nTo attempt to download these pathways please enter d\nTo skip this step just press enter.")){
				cat("\nFiles describing the networks for these pathways are not available on this computer:\n")
				print(neededPaths)
				cat("\nTo manually install these pathway diagram files, please follow these steps: \n1)obtain their biopax/pathway diagram files\n2)Place the biopax files in the ./reference_data/paths/biopax/ directory\n3)Add records of the newly added biopax files to the file ./reference_data/paths/biopax/record_of_biopax_pathways.txt\n")
				return(neededPaths)
			}
		}
	}else{
		return(neededPaths)
	}
	#3 obtain the reactome dbIDs
	cat("\nDownloading Reactome database ids...\n")
	dbidDict = getBiomartDbIds()
	#2 find pathways unavailable by regular avenues
	notAvailableDbIds = neededPaths[!neededPaths%in%dbidDict$pid]#figure out if any are not available
	if(length(notAvailableDbIds)) printNotAvailable(notAvail=notAvailableDbIds)
	#4
	neededPaths = neededPaths[neededPaths%in%dbidDict$pid]#get those that are needed
	
	#	download any needed pathways
	# 	testQuery = "http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level2/109581"#the last numbers are the dbid
	#	testDBid = 15869
	current_restfulCallbase = "http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/"
	old_restfulCall = "http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/"
	old_restfulCall2 = "http://www.reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/"
	# 	apoptosisPwId = "109581"
	# neededPaths = "Metabolism of nucleotides"
	
	restfulCall=current_restfulCallbase
	for(p in neededPaths){
		cat("\nAttempting to obtain pathway \'",p,"\'\n")
		dbid = dbidDict[p,"dbid"]
		pwUrl = paste(restfulCall, dbid, sep="")
		cat("Trying URL :", pwUrl, "\n")
		
		biopaxFile = try(getURL(url=pwUrl), silent=T)
		if(is.error.message(biopaxFile)){
			cat("\nError, could not dl the path ", p, "at the URL", pwUrl,"... Skipping to the next..\n")
		}else{
			addBiopaxPath(pname=p, dbid=dbid, biopaxDat=biopaxFile)
		}
	}
	return(notAvailableDbIds)
}#getReactomeBiopax

#pname : the name of the pathway
#dbid : the database id of the biopax record (with or with out .owl)
#biopaxDat : the actual biopax data. If this is given as a file name, the file will be moved to the correct biopax dir)
#dldate : the date the biopax file was downloaded or obtained. . 
addBiopaxPath<-function(pname, dbid, biopaxDat, dldate = as.character(Sys.time())){
	
	dbid = gsub(pattern=".owl|.OWL", replacement="", x=dbid)
	biopax.dir = "./reference_data/paths/biopax/"
	
	if(!file.exists(biopax.dir)) dir.create(path=biopax.dir, recursive=T, showWarnings=F)
	# find which pathways are needed
	pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	addToBiopaxRecord(pname=pname, dbid=dbid, dldate=dldate)
	if(is.null(dim(biopaxDat))){#if just a file name, copy the file
		
		dbid = basename(path=dbid)
		sourcePath  = basename(biopaxDat)
		bpath = dirname(pwrecord.fileName)
		file.copy(from=biopaxDat, to=paste0(bpath, "/",sourcePath))
		
	}else{#if it's a whole biopax record, write the file
		
		write.table(file=paste(biopax.dir, dbid, ".owl",sep=""),
								x=biopaxDat, sep="", row.names=F, col.names=F, quote=F)
		
	}
	
}

addToBiopaxRecord<-function(pname, dbid, dldate, pwrecord.fileName = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"){
	pwrecord = getPathwaysRecords(pwrecord.fileName=pwrecord.fileName)
	if(sum(pname%in%pwrecord$path_name)){ #remove the biopax file if it's already there
		pwrecord = pwrecord[!(pwrecord$path_name==pname),]
	}
	pwrecord=rbind.data.frame(pwrecord,
								 c(dbid,
								 	pname,
								 	dldate))
	colnames(pwrecord)<-c("dbID", "path_name", "download_date")
	#save the pwrecord
	write.table(x=pwrecord, file=pwrecord.fileName, quote=F, sep="\t", row.names=F, col.names=T)
}


printNotAvailable<-function(notAvail){
	cat("\nBiopax files for these pathways could not be downloaded:\n")
	print(notAvail)
	cat("\nTo display network diagrams for these pathways, please go\n",
			"to the reactome website and find the data base ID for these \n",
			"pathways. The dbIds for these pathways can be entered into\n", 
			"this tab-delimited file:\n",
			"'./reference_data/paths/biopax/record_of_biopax_pathways.txt'\n", 
			"along with the path's name and the date.\n\n",
			"After this information is entered in, the pathway network diagram\n",
			"display can be re-run and the pathways will be displayed.")
}


getBiomartDbIds<-function(){
	
	reac=reactomeBiomart()
	reac = useDataset(dataset="pathway", mart=reac)
	bmres1 = getBM(mart=reac, 
								 uniqueRows=T,
								 filters=c("species_selection"),
								 values=c("Homo sapiens"),
								 attributes=c("_displayname",
								 						 "stableidentifier_identifier",
								 						 "pathway_db_id"))
	colnames(bmres1)<-c("pid", "sid", "dbid")
	rownames(bmres1) <- bmres1[,1]
	return(bmres1)
}
reactomeBiomart<-function(verbose=T){
	
	reac = useMart("REACTOME")
	
	if(verbose) print(listDatasets(mart=reac))
	
	return(reac)
}

biopaxFileNameFromPathName<-function(pathNames, pwrecord.file = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"){
	if(!file.exists(pwrecord.file)){ 
		pwrecord.file = system.file("extdata/record_of_biopax_pathways.txt",package = "packageDir")
	}
	
	pwrecord = read.table(file=pwrecord.file, 
												quote="",
												comment.char="",
												stringsAsFactors=F,
												header=T,sep="\t")
	
	#make a file-path association dictionary
	rownames(pwrecord)<-pwrecord$path_name
	
	fileRoot = "./reference_data/paths/biopax/"
	
	availablePaths = pathNames[pathNames%in%rownames(pwrecord)]
	
	dbids = pwrecord[availablePaths,"dbID"]
	filenames = paste(fileRoot, dbids, ".owl", sep="")
	names(filenames)<-availablePaths
	return(filenames)
}

# Metabolism of nucleotides    pathway

UserProvidedBiopax<-function(pathNames, pathMetaDataFile=NULL){
	
	localPathwayRecord = "./reference_data/paths/biopax/record_of_biopax_pathways.txt"
	
	neededPaths=character(0)
	
	if(is.null(pathMetaDataFile)){ #if there's no path meta data file, just add the single path
		for(pn in pathNames){
			
			blnk = readline(paste("Please select the biopax file for the pathway '",
														pn,
														"'\n(Press enter to continue to file selection prompt,\nor enter s to skip this path)"))
			if(blnk=="s"){
				neededPaths = c(neededPaths, pn)
			}else{
				fname = file.choose()
				addBiopaxPath(dbid=fname, 
											biopaxDat=fname,
											pname=pn, 
											dldate=file.info(fname)$mtime)
			}
		}
	}else{ #path meta data passed, copy all files accordingly 
		
		#open the pathway record
		pwrecord = getPathwaysRecords(pwrecord.fileName=pathMetaDataFile)
		#clean the record: remove quotes
		pwrecord$path_name = gsub(pattern="\"",replacement="",x=pwrecord$path_name)
		#remove any duplicated records
		pwrecord = pwrecord[order(as.Date(pwrecord$download_date)),]
		pwrecord = pwrecord[!duplicated(pwrecord$dbID, fromLast=T),]
		
		#check that the pw record has correct formatting.. .
		if( sum(!c("dbID","path_name")%in%colnames(pwrecord)) ) stop(paste0("Missing columns from biopax pathway record file\n",
																																										"Columns found:\n", 
																																										paste(colnames(pwrecord),collapse=" "),
																																										"\nColumns needed:\ndbID path_name\n",
																																										"file name\n",pathMetaDataFile))
		#check the file names of the .owl files
		pwrecord$dbID = gsub(pattern=".owl$|.OWL$", replacement="", x=pwrecord$dbID)
		
		#copy each of the pathway files to the correct directory
		sourcePath  = dirname(pathMetaDataFile)
		bpath = dirname(localPathwayRecord)
		
		if(!file.exists(bpath)) dir.create(path=bpath, recursive=T, showWarnings=F)
		
		for(partname in pwrecord$dbID){
			curfname = paste0(bpath,"/", partname, ".owl")
			cursourcename = paste0(sourcePath, "/", partname, ".owl")
			file.copy(from=cursourcename, to=curfname)
		}
		
		neededPaths = pathNames[!toupper(x=pathNames)%in%toupper(x=pwrecord$path_name)]
		#copy pathway record file and all .owl files to correct directory 
		write.table(x=pwrecord, file=localPathwayRecord, quote=F,sep="\t", row.names=F, col.names=T)
	}
	return(neededPaths)
}

pathwaysFromBiopax<-function(study, pathNames, resSetName){
	
	#check if they can be obtained from Reactome, and try to get them if they are
	if(fromReactome(study)){
		notAvailablePaths = getReactomeBiopax(study=study, pathNames=pathNames)
	}else{
		UserProvidedBiopax(pathNames=pathNames)
	}
	
	availablePaths = setdiff(pathNames, notAvailablePaths)
	#get the names of the biopax files aquired and saved by the getReactomeBiopax function; retreives this as a dictionary
	#with the names as the path names and the values the file names
	bpfnames = biopaxFileNameFromPathName(pathNames=availablePaths)
	
	#interface with cytoscape
	for(i in 1:length(bpfnames)){
		# 		bioPaxToCytoscape(fname=bpfnames[i], pathName=names(bpfnames)[i])
		w = try(sendNetToCytoscape(fname=bpfnames[i], 
															 pathname=names(bpfnames)[i],#paste(names(bpfnames)[i],": from R overlap analysis"),
															 study=study), silent=T)
		if(is.error(w)){
			cat("\nError while trying to display the",bpfnames[i],"pathway\n")
			logError(elist=list(time=Sys.Date(),
										name="trying sendNetToCytoscape", 
										errorText=as.character(w), 
										pathName = names(bpfnames)[i], 
										pathFileName=bpfnames[i]))
		}else{
			#color the nodes appropriately
			setDefaultCytoscapeStyle(pname=names(bpfnames)[i], study=study,
															 w=w)
			setLayout(tw=w, placeInOrganelles=F)
			setBiologicalDataStyles(pname=names(bpfnames)[i], 
															study=study, 
															w=w, 
															resSetName=resSetName)	
		}
	}#for each biopax file name
	print("Done outputting pathways from biopax files.")
	return(notAvailablePaths)
}#pathwaysFromBiopax



setLayout<-function(tw, placeInOrganelles=F){
	layoutNetwork(obj=tw, layout.name="kamada-kawai")
	if(placeInOrganelles) organelleLayout(tw=tw)
}


organelleLayout<-function(tw){
	
	nodeTable = getAllNodeAttributes(obj=tw)
	cat("\nInferring cellular location of nodes with location provided\n")
	nodeTable = inferPositions(w=tw, nodeTable=nodeTable)
	cat("Establishing organelle spaces\n")
	organLayout = getSpaceDividers(nodeTable=nodeTable,w=tw)
	if(nrow(organLayout)){
		cat("Placing organelles\n")
		placeOrgans(opos=organLayout, tw=tw)
		cat("Distributing network in appropriate organelles\n")
		placeNodesInOrgans(organLayout=organLayout, nodeTable=nodeTable, tw=tw)
	}
}

placeNodesInOrgans<-function(organLayout, nodeTable, tw){
	rownames(organLayout)  = organLayout$names
	allLocations = unique(nodeTable$cellularLocation)
	for(loc in allLocations){
# 		layoutNetwork(obj=tw, layout.name="circular")
#   	selectNodes(obj=tw, 
#   							preserve.current.selection=F,
# 								node.names=getNodeNames(attval=loc, nodeTable=nodeTable, attname="cellularLocation"))
#  		setLayoutProperties(obj=tw, 
#  												layout.name='force-directed', 
#  												properties.list=list(selected_only=TRUE, defaultNodeMass=10))
#  		layoutNetwork(obj=tw, layout.name="force-directed")
		nodePosTab  = getNodePositionTable(nodeTable=nodeTable, 
																			 attval=loc, 
																			 attname="cellularLocation", 
																			 w=tw)
		cur = organLayout[loc,]
		cur[2:length(cur)] = as.numeric(cur[2:length(cur)])
		destBounds = list(x=c(cur$x, cur$x+cur$width), y=c(cur$y, cur$y+cur$height))
		expandToFill(tw=tw, nodePosTab=nodePosTab, toBounds=destBounds)
		
	}

}

expandToFill<-function(tw, nodePosTab, toBounds=list(x=c(1000,2000),y=c(1000,2000))){
	
	xrange = range(nodePosTab$x)
	yrange = range(nodePosTab$y)
	
	newXpoints = scaleAndMove(oldBound = xrange, newBound = toBounds$x, points=nodePosTab$x)
	newYpoints = scaleAndMove(oldBound = yrange, newBound = toBounds$y, points=nodePosTab$y)
	
	setNodePosition(obj=tw, node.names=nodePosTab$id, x.coords=newXpoints, y.coords=newYpoints)
	
}

# getCanvasSize<-function(w){
# 	
# 	nodesPos = getNodePosition(obj=w, node.names=rownames(nodeTable))
# 	
# 	nodePosTab = matrix(data="", nrow=length(nodesPos), ncol=3, dimnames=list(names(nodesPos), c("id","x","y")))
# 	
# 	for(nn in names(nodesPos)){
# 		nodePosTab[nn,] = c(nn,nodesPos[[nn]]$x, nodesPos[[nn]]$y)
# 	}
# 	nodePosTab = as.data.frame(nodePosTab, stringsAsFactors=F)
# 	nodePosTab$x = as.numeric(nodePosTab$x)
# 	nodePosTab$y = as.numeric(nodePosTab$y)
# 	
# }

#get a table with the position of all nodes in the graph w
getNodePostionTable<-function(w){
	
	nodesPos = getNodePosition(obj=w, node.names=getAllNodes(obj=w))
	
	nodePosTab = matrix(data="", nrow=length(nodesPos), ncol=3, dimnames=list(names(nodesPos), c("id","x","y")))
	
	for(nn in names(nodesPos)){
		nodePosTab[nn,] = c(nn,nodesPos[[nn]]$x, nodesPos[[nn]]$y)
	}
	nodePosTab = as.data.frame(nodePosTab, stringsAsFactors=F)
	nodePosTab$x = as.numeric(nodePosTab$x)
	nodePosTab$y = as.numeric(nodePosTab$y)
	return(nodePosTab)
}


getSpaceDividers<-function(nodeTable, w){
	
	npt = getNodePostionTable(w)
	
	pSize = c(max(npt$x)-min(npt$x), max(npt$y)-min(npt$y))
	
	minSize = 150
	cytoHeight = 0
	spacing = 5
	locTab = table(nodeTable$cellularLocation)
	
	adjustForCytosol = "cytosol"%in%names(locTab)
	
	if(adjustForCytosol & length(locTab)==1){
		retval = spacing
		names(retval)="cytosol"
		return(data.frame())
	}
	vsplit = 0
	cytoRatio=0
	if(adjustForCytosol){
		vsplit = (1-(locTab["cytosol"]/sum(locTab)))*pSize[1]
		#first, pull out the cytoplasm, if it's there
		organs = locTab[names(locTab)!='cytosol']
	}else{
		organs = locTab
	}
	
	numberOfBoxes = length(organs)
	
	boxRatios = organs/sum(organs)
	boxRatios = boxRatios[order(names(boxRatios),decreasing=T)]
	
	freeSpace = pSize[2] - (minSize*length(boxRatios) + (spacing*length(boxRatios) + 1))

	organHeight = (boxRatios*freeSpace) + minSize
	
	yPosition = c()
	place = spacing
	
	for(on in names(organHeight)){
		yPosition = c(yPosition, place)
		place = place + organHeight[on] + spacing
	}
	names(yPosition) = names(organHeight)
	xPosition = rep(spacing, times=length(yPosition))
	organWidth = rep((vsplit-spacing),times=length(organHeight))
	
	if(adjustForCytosol){
		cspace = spacing
		csize = (pSize[2]-(spacing*2))
		names(cspace)<-'cytosol'
		names(csize)<-'cytosol'
		yPosition = c(yPosition, cspace)
		xPosition = c(xPosition, ( spacing + vsplit))
		organHeight = c(organHeight, csize)
		organWidth = c(organWidth, (pSize[1]-organWidth[1]-spacing))
	}
	#add center:
	retval = cbind.data.frame(names=names(organHeight),
														x=xPosition, 
														y=yPosition, 
														height=organHeight, 
														width=organWidth, 
														stringsAsFactors=F)
	retval = cbind.data.frame(retval,
														centerX = retval$x+(retval$width/2),
														centerY = retval$y+(retval$height/2))

	#add a cell wall:
	cwPos = c(pSize[1]/2, pSize[2]/2)
	
	retval = rbind.data.frame(retval, c("cell wall", -5, -5, pSize[2]+5, pSize[1]+5, cwPos))
	rownames(retval)[nrow(retval)]<-"cell wall"
	return(retval)
}


placeOrgans<-function(opos, tw){

	gtmp = tw@graph 
	nlab = noa(graph=gtmp, node.attribute.name="label")

	gtmp = addNode(node=opos$names, object=gtmp)
	
	nodeData(gtmp,n=opos$names, attr="label")<-opos$names
	for(nn in opos$names){
		addCyNode(obj=tw, nodeName=nn)
	}

	tw = setGraph(obj=tw, graph=gtmp)
	
	redraw(tw)
	lockNodeDimensions(obj=tw, new.state=FALSE)

	setNodeWidthDirect(obj=tw, node.names=opos$names, new.widths=as.numeric(opos$width))
	
	setNodeHeightDirect(obj=tw, node.names=opos$names, new.heights=opos$height)

	setNodeShapeDirect(obj=tw, node.names=opos$names, new.shapes="round_rect")

	tmp = noa(graph=tw@graph, node.attribute.name="nodeType")
	tmp[tmp==""] = "cellularCompartment"
	# noa(graph=tw@graph, node.attribute.name="nodeType")<-tmp
	setNodeAttributesDirect(obj=tw, attribute.name="nodeType", attribute.type='char', node.names=names(tmp), values=tmp)
	
	setNodeOpacityDirect(obj=tw, node.names=opos$names, new.values=0)
	setNodeBorderOpacityDirect(obj=tw, node.names=opos$names, new.values=255)
	setNodeBorderWidthDirect(obj=tw, node.names=opos$names, new.sizes=5)

	if('cytosol'%in%opos$names){
		#make the cytosol box not have a border
		setNodeBorderOpacityDirect(obj=tw, node.names='cytosol', new.values=0)
	}
	redraw(tw)

	setNodeLabelDirect(obj=tw, node.names=opos$names, new.labels=opos$names)
	setNodeLabelColorDirect(obj=tw, node.names=opos$names, new.colors=rep("#000000", times=length(opos$names)))
	setNodeFontSizeDirect(obj=tw, node.names=opos$names, new.sizes=rep(30, times=length(opos$names)))
	setNodeLabelOpacityDirect(obj=tw, node.names=opos$names, new.values=rep(100, times=length(opos$names)))
	
# 	setNodeAttributesDirect(obj=tw, values=c(noa(graph=tw@graph, node.attribute.name="nodeType"), opos$names), 
# 													attribute.name="label", 
# 													attribute.type="char", 
# 													node.names=names(noa(graph=tw@graph, node.attribute.name="label")))
# 	setNodeBorderColorDirect(obj=tw, node.names=opos$names, new.color="#000000")
	setNodePosition(obj=tw, node.names=opos$names, x.coords=opos$centerX, y.coords=opos$centerY)
	redraw(tw)

}


createGraph<-function(nodes){
	
	gr <- new("graphNEL", edgemode = "directed")
	for(n in nodes){
		print(n)
		gr<-graph::addNode(n, gr)
	}
	return(gr)
}

scaleAndMove<-function(oldBound, newBound, points){
	
	olen = oldBound[2]-oldBound[1]
	nlen = newBound[2]-newBound[1]
	
	targetCenter = mean(newBound)
	
	if(olen==0) olen=1
	#adjust zoom
	expansionFactor = (nlen/olen)*.8
	#first, subtract all the offsets from the points
	offset1 = min(points)
	points = points - min(points)
	#then expand them based on their positions * expansion factor
	points = points*expansionFactor

	#move the set of points to the correct place
	#get current center; the mean of the min and the max
	currentCenter = mean(c(min(points), max(points)))
	
	#adjust location
	offset = targetCenter - currentCenter
	
	points = points + offset
	
	return(points)
}

getNodeNames<-function(attval, attname, nodeTable){
	
	nodeSet = nodeTable$nodeID[nodeTable[[attname]]==attval]
	return(nodeSet)
}

#'@title Get position of nodes in cytoscape diagram. . 
#'@description Gets a table with the positions of all the nodes. 
#'@param nodeTable The table of nodes
#'@param attval Used to get positions of a subset of nodes. Ex: if the attname is cellularLocation, and the attval is cell_membrane, the nodes with the cellularLocation attribute set to cell_membrane would be returned. 
#'@param attname The name of the attribute to filter the nodes by
#'@param w The cytoscape window connection object. 
#'@return data table with three columns: node id, x position, y postion
getNodePositionTable<-function(nodeTable, attval, attname, w){
	
	nodeSet = nodeTable$nodeID[nodeTable[[attname]]==attval]
	nodesPos = getNodePosition(obj=w, node.names=nodeSet)
	
	nodePosTab = matrix(data="", nrow=length(nodesPos), ncol=3, dimnames=list(names(nodesPos), c("id","x","y")))
	
	for(nn in names(nodesPos)){
		nodePosTab[nn,] = c(nn,nodesPos[[nn]]$x, nodesPos[[nn]]$y)
	}
	nodePosTab = as.data.frame(nodePosTab, stringsAsFactors=F)
	nodePosTab$x = as.numeric(nodePosTab$x)
	nodePosTab$y = as.numeric(nodePosTab$y)
	return(nodePosTab)
}

moveNodeSet<-function(nodeTable, w, 
											att="cellularLocation", 
											attname="mitochondrial matrix", 
											quadrent = "I"){
	
	nodePosTab  = getNodePositionTable(nodeTable=nodeTable, attname=attname, attval=att, w=w)
	
	cent = getCenter(w)
	
	yoffset = nodePosTab$y/2
	xoffset = nodePosTab$x/2
	
	newCent = cent
	
	movement = list()
	movement$x = switch(quadrent, "I"=cent$x, "II"=0, "III"=0,"IV"=cent$x)
	movement$y = switch(quadrent, "I"=0, "II"=0, "III"=cent$y,"IV"=cent$y)
	
	ynew = yoffset + movement$y
	xnew = xoffset + movement$x
	
	setNodePosition(obj=w,node.names=nodePosTab$id, x.coords=xnew, y.coords=ynew)
}

findRows<-function(df, val){
	
	for(cn in colnames(df)){
		print(cn)
		rows = grep(pattern=val, df[,cn], ignore.case=T)
		if(length(rows)){
			
			print(df[rows,])
			
		}else{ print("   no matches found")}
		readline("press enter to continue to the next column")
	}
	
}

is.error<-function(testForError){
	grepl(pattern="error", x=class(testForError), ignore.case=T)
}


# logError(list(time=Sys.Date(),
# 							name="trying sendNetToCytoscape", 
# 							errorText=character(w), 
# 							pathName = names(bpfnames)[i], 
# 							pathFileName=bpfnames[i]))
logError<-function(elist,elogDir="./errorLog/"){
	
	dir.create(elogDir, showWarnings=F, recursive=T)
	efname = paste0(elogDir, 
									ifelse(test=is.null(elist$name), 
												 yes=as.character(Sys.Date()), 
												 no=elist$name), 
									".txt")
	cat("writing error to file", efname,"\n")
	for(i in 1:length(elist)){
		el=as.character(elist[[i]])
		print(el)
		write.table(x=paste(names(elist)[i], Sys.time()), 
								file=efname, 
								append=file.exists(efname), 
								sep="\t", 
								quote=F, 
								row.names=F, 
								col.names=F)
		write.table(x=el, 
								file=efname,
								append=T, 
								sep="\t", 
								quote=F, 
								row.names=F, 
								col.names=F)
	}
	
}



inferPositions<-function(w, nodeTable){
	
	#get the nodes with location not set
	nopos = nodeTable$nodeID[nodeTable$cellularLocation=="notProvidedInCellNetworkData"]
	
	#get location of connected nodes
	for(id in nopos){
		neigh = getFirstNeighbors(w, node.names=id)
		cellposV = nodeTable[neigh,"cellularLocation"]
		postab = table(cellposV)
		postab = postab[names(postab)!="notProvidedInCellNetworkData"]#make sure this doesn't end up being the majority vote
		consense = names(postab)[postab==max(postab)]
		nodeTable[id,"cellularLocation"] = consense[1]
	}
	return(nodeTable)
}


appendAltNames<-function(nodeTable, altnames){
	
	agres = aggregate(x=altnames$altName, by=list(altnames$nodeID), FUN=function(x){
		paste(x, collapse=" | ")
	})
	colnames(agres)<-c("nodeID", "Alternate.Names")
	out = merge(x=nodeTable, y=agres, by="nodeID", all.x=TRUE)
	out[is.na(out[,"Alternate.Names"]),"Alternate.Names"] = "No alternate names found in Biopax file"
	rownames(out)<-out$nodeID
	return(out)	
	
}

pullOutNodeTable<-function(df, nodeTypes=NULL){

	if(is.null(nodeTypes)){
		nodeTypes=c("BiochemicalReaction", "Catalysis", "Complex", "Protein", "SmallMolecule","PhysicalEntity")
		print(nodeTypes)
		readline("\nWarning, only using the above listed default node types\n")
	} 
	
	#get all the rows with the node types
	tdf = df[df$class%in%nodeTypes,]
	idType = unique(tdf[,c("id","class")])
	uids = idType$id
	
	##########################################
	#make an out data frame with these columns:
	##########################################
	#row names = nodeIDs
	#displayName
	#nodeType
	#cellular compartment
	#evidence: from these classes: PublicationXref
	#date: from the dates given by PublicationXrefs
	#comment: from the comments given for nodes
	#alternate.names: from references, names and display name
	####################################
	outcolnames = c("nodeID","nodeType", "displayName", "cellularLocation", "evidence", "date", "comment")
	#dataSource (property=="dataSource"/the value in the property_attr_value column will have a "#Provenance<#> id)
	outdf = data.frame(stringsAsFactors=F,
										 matrix(data="", 
														nrow=length(uids), 
														ncol=length(outcolnames), 
														dimnames=list(uids, outcolnames)))
	outdf$nodeID=rownames(outdf)
	outdf$nodeType = idType$class
	###################### set up the display names
	outdf = getDisplayNames(tdf=tdf, df=df, outdf=outdf)
	###################### set up the cellular locations
	outdf = setCellularLocations(tdf=tdf, outdf=outdf, df=df)
	################### set up the evidence
	nodeCitations = getPublicationRefs(df=df, tdf=tdf)
	outdf[nodeCitations$id,c("evidence","date")] = nodeCitations[,c("citation","date")]
	################### set up any comments
	nodeComments = getComments(tdf=tdf)
	outdf[nodeComments$id,"comment"] = nodeComments$comment
	
	return(outdf)	
}#pullOutNodeTable

setCellularLocations <- function (tdf, outdf, df) {
	cellularLocex = tdf[tdf$property=="cellularLocation",c("id", "property_attr_value")]
	outdf[cellularLocex$id,"cellularLocation"] = cellularLocex$property_attr_value
	#make a dictionary of the locations
	locDictTmp = df[df$class=="CellularLocationVocabulary"&df$property=="term",c("id","property_value")]
	locDict = locDictTmp$property_value
	#put the locations in the appropriate parts of the outdf
	names(locDict) = paste0("#",locDictTmp$id)
	
	outdf$cellularLocation = locDict[outdf$cellularLocation]
	outdf$cellularLocation[is.na(outdf$cellularLocation)] = "notProvidedInCellNetworkData"
	return(outdf)
}#setCellularLocations


getDisplayNames<-function(tdf, df, outdf){
	dispNameEx = tdf[tdf$property=="displayName",c("id", "property_value")]
	outdf[dispNameEx$id,"displayName"] = dispNameEx$property_value
	#adjust the display names of the Catalyses
	catTmp = tdf[tdf$class=="Catalysis"&tdf$property=="controlType",c("id","property_value")]
	if(nrow(catTmp)){
		outdf[catTmp$id,"displayName"] = catTmp$property_value
	}
	#if any display names were missing, fill them in with the nodeID
	outdf[outdf$displayName=="","displayName"] = outdf$nodeID[outdf$displayName==""]
	return(outdf)
}



#'@title getPublicationRefs()
#'@param df the main biopax data frame
#'@param tdf a subsegment of the biopax data frame from which citations and citation dates are to be found. Must have columns "property", "property_value" and "id".
#'@return data frame with columns "id"       "citation" "date"
getPublicationRefs<-function(df, tdf){
	
	# 	dfpex = df[df$property=="xref"&df$class=="BiochemicalReaction",]
	#get the xref rows
	tdfxref = tdf[tdf$property=="xref",]
	#get the pub rows
	pxr = df[df$class=="PublicationXref",]
	
	tdf$property_attr_value = gsub(pattern="^#",replacement="",x=tdf$property_attr_value)
	#find the overlap after stripping any hashes off of the property_attr_value
	tdfxref$property_attr_value = gsub(pattern="^#",replacement="",x=tdfxref$property_attr_value)
	#the rows with publication information specific to the nodes in the pathway
	pubRows = pxr[pxr$id%in%tdfxref$property_attr_value,]
	#now group 
	ppxr = aggregate(x=pubRows$property_value, by=list(pubRows$id), 
									 FUN=function(x){paste(x, collapse=" | ", sep=" | ")})
	#connect the publications to the nodes
	nodesPlusCitation = merge(x=tdf, y=ppxr, by.x="property_attr_value", by.y="Group.1")
	
	agCitation = aggregate(x=nodesPlusCitation$x, by=list(nodesPlusCitation$id), FUN=function(x){paste(x, collapse=" || ", sep=" || ")})
	colnames(agCitation)<-c("id", "citation")
	
	#get the dates related to nodes
	pxrDates = pubRows[pubRows$property=="year",c("id", "property_value")]
	names(pxrDates) = c("id", "year")
	#merge them with the nodes
	datesMergedToNodes = merge(x=tdf, y=pxrDates, by.x="property_attr_value", by.y="id")
	datesMergedToNodes = datesMergedToNodes[,c("id","year")]
	
	agDates = aggregate(x=datesMergedToNodes$year, by=list(datesMergedToNodes$id), FUN=function(x){paste(x, collapse=" ; ", sep=" ; ")})
	names(agDates)<-c("id","date")
	
	citAndDate = merge(x=agCitation, y=agDates, by="id", all.x=T)
	
	return(citAndDate)
}#getPublicationRefs


#'@title getComments()
#'@description Gets and pastes together comments for all elements in the data.frame provided by rBiopaxParser's readBiopax function.
#'@param tdf data.frame, the section of the biopax data frame to pull comments about nodes from. Must have columns "property", "property_value" and "id"
#'@return data frame with two columns "id" and "comment" where multiple comments are sepparated by pipes (ie, "|")
getComments<-function(tdf){
	#pull the comments
	commentRows = tdf[tdf$property=="comment",]
	nodeToComments = aggregate(x=commentRows$property_value, by=list(commentRows$id), FUN=function(x){paste(x, collapse=" | ", sep=" | ")})
	colnames(nodeToComments)<-c("id","comment")
	return(nodeToComments)
}

test.sendNetToCytoscape<-function(){
	
	#bioPaxToCytoscape
	pathName="Abacavir metabolism test"
	fname = "/Users/samhiggins2001_worldperks/tprog/main_131219/biopax_test/abacavir_metabolism_Reactome_2161522_biopax3.owl"
	
	sendNetToCytoscape(fname=fname, pathname=pathName)
	
	if(is.null(fname)) fname = file.choose()
	biopax = readBiopax(file=fname, verbose=T)
	print(biopax)
	
	
	pathName2="Nucleotide metabolism test internal"
	fname2 = "/Users/samhiggins2001_worldperks/tprog/main_131219/reference_data/paths/biopax/15869.owl"
	sendNetToCytoscape(fname=fname2, pathname=pathName2)
	
# 	directPath = "/Users/samhiggins2001_worldperks/tprog/main_131219/reference_data/paths/biopax/499943.owl"
	print("Trying pathway from other source...")
	
	directPath= "./biopax_test/NetworkFromPathwayCommons.owl"
	sendNetToCytoscape(fname=directPath, pathname="NetworkFromPathwayCommons", study=STUDY)
	
	directPath2="/Users/samhiggins2001_worldperks/tprog/main_131219/biopax_test/pathsBetween_camk2d_IFNG.OWL"
	sendNetToCytoscape(fname=directPath2, pathname="NetworkFromPathwayCommons", study=STUDY)

}


sendNetToCytoscape<-function(fname, pathname, study){
	cat("\nOutputting pathway", pathname, "to Cytoscape,\n")
	cat("from file:", fname, "\n")
	cat("Loading Biopax file...\n")
	biopax = readBiopax(file=fname, verbose=T)

	#adding 
	gnelmodel = biopaxToGraphNEL(biopax=biopax, STUDY=study)

	while(T){
		cw = try(graphNELToCytoscape(gnel=gnelmodel, pathname=pathname), silent=T)
		if(is.error(cw)){
			print("Error in sending pathway to Cytoscape..")
			print(cw)
			if(grepl(pattern="connect", x=as.character(cw), ignore.case=T)){#if it's a connection error
				print(as.character(cw))
				uintmp = readline("It appears there was a problem with the program's connection to Cytoscape.\nPlease check that Cytoscape is running and the cytoscapeRPC plugin is activated.")
			}else{#if it's not a connection error, pass the error on up
				break
			}
		}else{
			break
		}
	}
	return(cw)
}#sendNetToCytoscape


findAllInteractionRows<-function(df, nodeTypes){

	ex1 = df[df$class%in%nodeTypes,]
	ex1ids = unique(ex1$id)
	ex1ids = paste0("#", ex1ids)
	
	irows = df[(df$property_attr_value%in%ex1ids)&(df$property!="comment")&df$class%in%nodeTypes,]
	
	return(irows)
}


biopaxToGraphNEL<-function(biopax, STUDY){
	
	pname = names(biopax$file)
	df = unfactorize(biopax$df)
	
	physicalEntities=c("Complex", 
										"DNARegion",
										"DNA", 
										"Protein",
										"RNA",
										"RNARegion",
										"SmallMolecule",
										"PhysicalEntity")
	
	controlSubclasses= c("Catalysis", 
											 "Modulation", 
											 "TemplateReactionRegulation")
	
	conversionSubclasses = c("BiochemicalReaction", 
													 "ComplexAssembly", 
													 "Degradation",
													 "Transport",
													 "TransportWithBiochemicalReaction")
	
	interactionClasses = c("TemplateReaction", 
												 "MolecularInteraction", 
												 "GeneticInteraction", 
												 "Conversion", 
												 "Control")
	
	nodeTypes  = c(physicalEntities, controlSubclasses, conversionSubclasses, interactionClasses)
# 							"PhysicalEntity", 
 #TemplateReaction, TemplateReactionRegulation GeneticInteraction Conversion (subclass of degredation)
	
	#get the names of the nodes
	nodeTable0 = pullOutNodeTable(df=df,nodeTypes=nodeTypes)
	
	nodeTable = fixGeneNames(nodeTable=nodeTable0, df=df, STUDY=STUDY, pname=pname)

	############## make the graphNEL object
	#set the nodes
	gnelmodel = graphNEL(nodes=nodeTable$nodeID, 
											 edgemode="directed")
	#set the edges
	edgeDat = setEdges(df=df, nodeTypes=nodeTypes, pre=gnelmodel)
	gnelmodel = addEdge(from=edgeDat[,"subject"], 
											to=edgeDat[,"object"], 
											graph=gnelmodel, 
											weights=1)

	gnelmodel = initEdgeAttribute(graph=gnelmodel, 
																attribute.name="edgeType", 
																attribute.type="char", 
																default.value="associated")

	gnelmodel = initEdgeAttribute(graph=gnelmodel, attribute.name="weight", 
																attribute.type="numeric", default.value=1)

	############## establish edge attributes
	edgeData(self=gnelmodel, 
					 from=edgeDat[,"subject"],
					 to=edgeDat[,"object"],
					 attr="edgeType")<-edgeDat[,"direction"]

	############## establish node attributes
	colnames(nodeTable)[colnames(nodeTable)=="displayName"] ="label"
	for(cn in colnames(nodeTable)){
		gnelmodel = initNodeAttribute(graph=gnelmodel, attribute.name=cn, attribute.type="char", default.value="")
		nodeData(self=gnelmodel, n=nodeTable$nodeID, attr=cn)<-nodeTable[[cn]]
	}
	
	return(gnelmodel)
}#biopaxToGraphNEL


attemptFixSmallMoleculeNames<-function(xref, nodeTable){
	
	#get all the smallMolecule nodeNames
	smnn = nodeTable$nodeID[nodeTable$nodeType=="SmallMolecule"]
	
	#get all the xref rows corresponding to small molecule node names
	xrefSub = xref[xref$nodeID%in%smnn,]
	#find all the rows which are endorsed by [ChEBI:
	eRows = xrefSub[grepl(pattern="[[]ChEBI:", x=xrefSub$altName),]
	if(nrow(eRows)){
		eRows$altName = gsub(pattern="[ ]*[[]ChEBI:[0-9]+[]]", replacement="", x=eRows$altName )
		nodeTable[eRows$nodeID,"displayName"] = eRows$altName
	}
	return(nodeTable)
}
# 
# attemptFixSmallMoleculeNames<-function(xref, nodeTable){
# 	
# 	#get all the smallMolecule nodeNames
# 	smnn = nodeTable$nodeID[nodeTable$nodeType=="SmallMolecule"]
# 	
# 	#get all the xref rows corresponding to small molecule node names
# 	xrefSub = xref[xref$nodeID%in%smnn,]
# 	#find all the rows which are endorsed by [ChEBI:
# 	eRows = xrefSub[grepl(pattern="[[]ChEBI:", x=xrefSub$altName),]
# 	
# 	eRows$altName = gsub(pattern="[ ]*[[]ChEBI:[0-9]+[]]", replacement="", x=eRows$altName )
# 	
# 	nodeTable[eRows$nodeID,"displayName"] = eRows$altName
# 	return(nodeTable)
# }

fixGeneNames <- function (nodeTable, df, STUDY, pname) {
	
	pathGenes = getGenesFromPaths(pids=pname, STUDY=STUDY)
	
	allAltNames = as.data.frame(getAllAltNames(nodeIDs=nodeTable$
																						 	nodeID, df=df), stringsAsFactors=F)
	
	nodeTable = appendAltNames(nodeTable=nodeTable, altnames=allAltNames)
	nodeTable = attemptFixSmallMoleculeNames(xref=allAltNames, nodeTable=nodeTable)
	
	#see if any of these cannot be found in the biopax displayNames
	pathIdsMissingFromBiopax  = setdiff(pathGenes, nodeTable$displayName)
	
	if(length(pathIdsMissingFromBiopax)==0) return(nodeTable) #if they're all correct, return the original node table
	
	#pull the possible matches out of the biopax df
	dontCorrect = c("Degradation", 
									"Catalysis", 
									"BiochemicalReaction", 
									"SmallMolecule", 
									"TemplateReactionRegulation", 
									"Modulation", 
									"Transport", 
									"TransportWithBiochemicalReaction")
	#look at only those rows which are not obviously not genes, and which are not found in the path
	toBeCorrected = !nodeTable$nodeType%in%dontCorrect
	notInPathi = !nodeTable$displayName%in%pathGenes
	
	#figure out which nodes need attention
	toCorrectInNodeTablei=toBeCorrected&notInPathi#the index of displayNames not in the pathway, with no spaces and not in the pathway
	nonHugoNodeIds = nodeTable$nodeID[toCorrectInNodeTablei]
# 	nonHugoSymbols = nodeTable$displayName[toCorrectInNodeTablei]
# 	idSymDict = nonHugoSymbols
# 	names(idSymDict) = nonHugoNodeIds #idSymDict: dictionary linking the node ids to their errant, non-hugo symbols

	idnameXref = allAltNames[allAltNames$nodeID%in%nonHugoNodeIds,]

	idnameXref[,2]=toupper(x=idnameXref[,2])

	#check if any of the pathIdsMissingFromBiopax are in the idnameXref
	nodeTable = checkUseAlternates(missing=pathIdsMissingFromBiopax, xref=idnameXref, nodeTable=nodeTable)
	
	foundi = pathIdsMissingFromBiopax%in%idnameXref[,2]

	pathIdsMissingFromBiopax  = setdiff(pathGenes, nodeTable$displayName)
	if(length(pathIdsMissingFromBiopax)==0){
		nodeTable=addModelMembershipColumn(nodeTable=nodeTable, pathGenes=pathGenes, stillmissing=pathIdsMissingFromBiopax, toCorrectInNodeTablei=toCorrectInNodeTablei)
		return(nodeTable)#if that fixes it, exit
	} 
	#if here, there are still nodes from the path that cant be found in the biopax file
	#next, try to correct the hugo symbols
	# 	update the nonHugoNodeIds, nonHugoSymbols and idSymDict to reflect those still not in the path
	notInPathi = !nodeTable$displayName%in%pathGenes
	notHugoi = !nodeTable$displayName%in%STUDY@studyMetaData@paths$symtable$Approved.Symbol
	toCorrectInNodeTablei=toBeCorrected&notInPathi&notHugoi
	nonHugoNodeIds = nodeTable$nodeID[toCorrectInNodeTablei]
	#		attempt to correct the alternate symbols to HUGO
	idnameXref[,2] = corsym(symbol_set=idnameXref[,2], symref=STUDY, verbose=F)
	idnameXref = unique(idnameXref)#make the rows unique

	nodeTable=checkUseAlternates(missing=pathIdsMissingFromBiopax, xref=idnameXref, nodeTable=nodeTable)
	
	pathIdsMissingFromBiopax  = setdiff(pathGenes, nodeTable$displayName)
	if(length(pathIdsMissingFromBiopax)==0){
		nodeTable=addModelMembershipColumn(nodeTable=nodeTable, pathGenes=pathGenes, stillmissing=pathIdsMissingFromBiopax, toCorrectInNodeTablei=toCorrectInNodeTablei)
		return(nodeTable)#if that fixes it, exit
	}
	print("!!!!!!!!!!!!!!!!!!!!!!!!!note: not all gene symbols from the original path model could be found in the biopax model provided!!!!!!!!!!!!!!!")
	#if the function gets to this point there are still path ids not found, so deal with that appropriately 
	#figure out which symbols are still missing from the model
	stillmissing = setdiff(pathGenes, nodeTable$displayName)
	cat("These are the gene symbols not found in the biopax file:\n",stillmissing,"\n")
	#add in-model column
	#append column to the nodeTable
	nodeTable=addModelMembershipColumn(nodeTable=nodeTable, pathGenes=pathGenes, stillmissing=pathIdsMissingFromBiopax, toCorrectInNodeTablei=toCorrectInNodeTablei)
	
	return(nodeTable)
}#fixGeneNames

#adds column to nodeTable indicating if the graph node is found in the 
addModelMembershipColumn <- function (nodeTable, pathGenes, stillmissing, toCorrectInNodeTablei) {

	inOriginalModelColumn = rep("entity symbol not in pathway gene set", times=nrow(nodeTable))
	names(inOriginalModelColumn)<-nodeTable$nodeID
	inModeli = nodeTable$displayName%in%pathGenes
	inOriginalModelColumn[inModeli] = "entity symbol in pathway gene set"
	if(length(stillmissing)){#if some of the gene set is missing, fill in the in gene set column appropriately 
		problemIndexes = toCorrectInNodeTablei&!nodeTable$displayName%in%pathGenes
		inOriginalModelColumn[problemIndexes] = "node in quesion: one or more symbols from the pathway gene set could not be found in the biopax model"
	}
	nodeTable = cbind.data.frame(nodeTable, inOriginalModelColumn, stringsAsFactors=F)
	return(nodeTable)
}

checkUseAlternates<-function(missing, xref, nodeTable){
	#see if any of the gene identifiers are now in the path
	foundi = missing%in%xref[,2]
	if(sum(foundi)){#if the gene names in xref match any of the pathGenes not found in the biopax
		#then: 
		#		fix the biopax names
		fixedGeneNames = missing[foundi]
		#				pull the gene names out of the Xref
		xrefsForFixedGeneNames = xref[xref[,2]%in%missing,]
		#				add the appropriate reparied symbols to nodeTable
		nodeTable[xrefsForFixedGeneNames[,1],"displayName"] = xrefsForFixedGeneNames[,2]
		#		update the needed path genes
	}
	return(nodeTable)
}


getAllAltNames<-function(nodeIDs, df){
	
	#find possible corrections:
	#attempt to correct the names in the nonHugoNodeIds
	#find alternate names directly associated with the node
	idnameXref = getAlternateNames(nodeIds = nodeIDs, df=df)
	#get protein reference table; find alternate names via the protein cross reference
	reftab = as.matrix(getProteinReferenceTable(nodeIDs=nodeIDs, df=df))
	#combine the two sources of alternate names
	colnames(reftab)<-colnames(idnameXref)
	idnameXref = rbind(idnameXref, reftab)
	idnameXref = unique(idnameXref)
	
	return(idnameXref)
}

getProteinReferenceTable<-function(nodeIDs, df){
	#get node->reference dictionary
	rtab = df[df$property=="entityReference"&df$id%in%nodeIDs,]
	rtab$property_attr_value =	gsub(pattern="^#",replacement="", x=rtab$property_attr_value)
	rtab = rtab[,c("id","property_attr_value")]
	colnames(rtab)<-c("protein_id","ref_id")
	
	protrefs = df[df$id%in%rtab$ref_id&df$property=="name",]
	protrefs = merge(x=rtab, y=protrefs, by.x="ref_id", by.y="id")

	return(protrefs[,c("protein_id","property_value")])
}

# allowUserToFixBiopaxNames<-function(nodeTable, df, STUDY, pname){
# 	#for each non-matched name from the path, give the user
# 	#the genes remainig in the path
# 	#the set of alternate names
# 	#an option to input a new name
# 	pRem = setdiff(pathGenes, nodeTable$displayName)
# 	if(!length(pRem)) return(nodeTable)
# 	cat("\nThese genes in the current pathway do not have exact matches in the biopax file:\n")
# 	print(prem)
# 	symbolCorrectionsTable = matrix(data="",nrow=0,ncol=2)
# 	colnames(symbolCorrectionsTable)<-c("oldSymbol", "newSymbol")
# 	
# 	for(i in 1:length(pRem)){#for each index in the genes remaining in the pathway
# 		cat("Does the gene symbol,", pRem[i], ", from the current pathway\n",
# 				pname,
# 				"\nmatch any of the above described nodes in the biopax file?")
# 		print()
# 		
# 		uin = promptNumeric(prompt="If you find a match, enter the number corresponding to the match here.\nIf a match cannot be found enter a blank line.")
# 		if(uin==""){
# 			#if the user enters a blank, dont make a correction
# 		}else{
# 			#if the user enters a number make the corresponding correction
# 			badSym = nonHugoSymbols[uin]
# 			nonHugoSymbols[uin] = pRem[i]
# 			#then add to symbol corrections table
# 			symbolCorrectionsTable = rbind(symbolCorrectionsTable, c(badSym, nonHugoSymbols[uin]))
# 		}
# 	}
# 	if(nrow(symbolCorrectionsTable)){
# 		cat("\nThis is the list of symbols that were corrected\n")
# 		print(symbolCorrectionsTable)
# 		if(readline("Would you like to save this set of gene symbol corrections? \n(enter y to save, anything else not to save)")=="y"){
# 			#save the new corrections to the main symbol correction table
# 			addCorrections(new_corrections=symbolCorrectionsTable,correctionsfile="./testCorrectionsFile.txt")
# 		}
# 		
# 	}
# }


#'@title getAlternateNames
#'@description looks through the data frame assocaited with biopax file for alternate names associated with a node/entity
#'@param nodeIds The ids for the nodes whose alternate names are being searched for. 
#'@param df The data frame containing the biopax data, as given by the rBiopaxParser package function readBiopax
#'@return data.frame with columns nodeID and altName.
getAlternateNames<-function(nodeIds, df){
	ex1 = df[df$id%in%nodeIds & (df$property=="name"|df$property=="displayName"),]
	nlist = aggregate(x=ex1$property_value, by=list(ex1$id), FUN=paste)
	mout = matrix(data="", ncol=2,nrow=0, dimnames=list(NULL, c("nodeID","altName")))
	for(i in 1:nrow(nlist)){
		cur = nlist[i,]
		curmat = cbind(rep(cur[,1], time=length(cur[[1,2]])), cur[[1,2]])
		mout = rbind(mout,curmat)
	}
	return(mout)
}

addExperimentalDataToGNEL<-function(gnelmodel){
	gnelmodel = initNodeAttribute(graph=gnelmodel, 
																attribute.name="functionallyAffected",
																attribute.type="integer", 
																default.value=0)
	gnelmodel = initNodeAttribute(graph=gnelmodel, 
																attribute.name="aberrationallyAffected",
																attribute.type="integer", 
																default.value=0)
	return(gnelmodel)
}

graphNELToCytoscape<-function(gnel, pathname){
	
	cw <- new.CytoscapeWindow(pathname, graph=gnel, overwriteWindow=T)
	
	displayGraph(cw)
	
	
	setNodeLabelDirect(obj=cw, 
										 node.names=names(noa(graph=gnel, node.attribute.name="label")), 
										 new.labels=noa(graph=gnel, node.attribute.name="label"))
	redraw(cw)
	
	return(cw)
}


findInDf<-function(df, term="PhysicalEntity1", retval=F){
	
	indexes = c()
	for(cn in colnames(df)){
		print(cn)
		indexes = c(indexes, grep(pattern=term, x=df[[cn]], ignore.case=T))
	}
	
	print(df[indexes,])
	
	if(retval) return(df[indexes,])
}

setEdges<-function(df, pre, nodeTypes){
	
	edgeInteractionRows = findAllInteractionRows(df=df, nodeTypes=nodeTypes)
	
	edgeEx = edgeInteractionRows[,c("class", "id", "property", "property_attr_value")]
	edgeEx[,4] = gsub(pattern="^#", replacement="", x=edgeEx[,4])
	colnames(edgeEx)<-c("class","subject","direction","object")

	return(edgeEx)
}

checkForDups<-function(edg){
	
	rev = edg[,c(2,1,3)]
	
	d1  = duplicated(x=edg)
	d2  = duplicated(x=edg, fromLast=T)
	ad = d1|d2
	
	edg[ad,]
}


setNodeColors<-function(w, nodeLabels, color, defaultColor, study){
	#assure the nodes are in the pathway?
	
	
	setNodeColorRule(obj=w, node.attribute.name="label", 
									 default.color=defaultColor,
									 mode="lookup",
									 control.points=nodeLabels, 
									 colors=rep(color, times=length(nodeLabels)) )
	
}




