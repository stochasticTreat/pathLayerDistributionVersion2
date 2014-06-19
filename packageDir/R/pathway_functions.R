#pathway_functions.R
####################



#'@title Retreive the default set of cellular pathways. 
#'@description Retreives the default set of cellular pathways. These are the gene sets for all Reactome pathways as downloaded on 4/6/14 . 
#'@param path_file string giving the relative or absolute file path to an accessible set of cellular pathways. Set of pathways must first have been recorded in the reference_data/paths/pathMetaData.txt file. 
#'@return A \code{Path_Detail} object. 
#'@export
#'@examples
#'cellularPathwayGeneSets = getDefaultPaths()
getDefaultPaths<-function(path_file=system.file("extdata/Reactome.2014.04.06.12.52.27.txt", package = "packageDir")){

	path_detail_tmp<-getPaths(path_file=path_file, verbose=F)
	return(path_detail_tmp)
	
}

#pathsFromGraphite function
#takes path name
#returns nodes as HUGO symbols
getNodes<-function(path){
	#returns nodes in a pathway with their identifiers converted to hugo symbols
	return(nodes(convertIdentifiers(pathway=path, type="symbol")))
}

# importAllGraphite()
#'@title Import all pathway repositories from the graphite package to this package. 
#'@description Gives options to import pathway repositories from Reactome, Spike, NCI, KEGG and Biocarta. 
#'@param repositories A \code{list} object, each slot containing a pathway repository from graphite that is to be imported
#'@param symtab A HUGO symbol look up table, as provided by the package HGNChelper or from genenames.org. 
#'@param verbose A flag indicating if this function should be run in an interactive mode. This allows individual pathway repositories to be selected, insted of automatically selecting all pathways.
#'@import graphite
#'@import tools
importAllGraphite<-function(repositories = list(reactome=reactome, spike=spike, nci=nci, kegg=kegg, biocarta=biocarta), 
														symtab=NULL, 
														verbose=F){
	if(is.null(symtab)) symtab = getHugoSymbols()
	if(verbose){
		print(matrix(data=names(repositories), ncol=1, dimnames=list(1:length(repositories), "Repository Name")))
		while(T){
			line = readline("Please enter the number(s) of the pathway set(s) you would like to import from graphite\n (separated by spaces): ")
			line = try(as.numeric(strsplit(x=line, split=" ")[[1]]), silent=T)
			if(!sum(is.na(line))&!sum(!line%in%1:length(repositories))) break
			cat("\nSorry, that input was not understood, please try again.\n")
		}
		repositories = repositories[line]
	}
	for(pws in repositories){
		
		preped_paths = importFromGraphite(db=pws)
		#correct gene symbols
		colnames(preped_paths$paths)<-corsym(symbol_set=colnames(preped_paths$paths), 
																				 symref=symtab, 
																				 verbose=T)
		
		#setPathMetaData
		preped_paths = manualPathMetaData(preped_paths=preped_paths, symtab=symtab)
		
		#establish path set
		recordPathSet(ps=preped_paths)
	}
}

importFromGraphite<-function(db){

	outset = list()
	date = c()
	dbname = c()
	for(i in 1:length(db)){
		pname = names(db)[i]
		cat(pname, length(nodes(db[[pname]])), "\n")
		outset[[pname]] = getNodes(path=db[[pname]])
		
		date = union(date, db[[pname]]@timestamp)
		dbname	= union(dbname, db[[pname]]@database)
	}
	bip = list_to_table(pth=outset)
	
# 	neededSlots = c("name", "date", "source")
	preped_paths = list()
	preped_paths[["date"]] = paste(db[[1]]@timestamp, sep=" ", collapse=" ")
	preped_paths[["source"]] = paste(db[[1]]@database, "via the bioconductor graphite package", sep=" ", collapse=" ")
	preped_paths[["name"]] = db[[1]]@database
	
	preped_paths[["paths"]] = as.matrix(bip)
	preped_paths[["original_file_or_source"]] = paste(db[[1]]@database, "from Bioconductor package 'graphite', version", packageVersion("graphite") )
	preped_paths[["original_file_creation_date"]] = paste("graphite utilized on",as.character(Sys.time()))
	return(preped_paths)
}

#parseReactomeUniprot()
# loads paths from file in this format: <path name><\t><gene><space><gene>....
# everything after and including "firstGeneColum" is assumed to be a gene"
#the first column is assumed to be the path name or other unique identifier for the pathway
#
#takes: firstGeneColum: the first tab stop with gene names in it, usually 2 or 3
#				keepLowerCaseGenes:		keep genes who are all lower case
parseReactomeUniprot<-function(firstVectorColum, fname = NULL ){
	
	pset = read.table(fname, sep="\t",comment.char="", stringsAsFactors=F, quote="")
	outlist = list()
	
	for(i in 1:nrow(pset)){
		cur = pset[i,]
		#pull out the "[19 processes]:" part
		paths = gsub(pattern="[[0-9]+ processes]: ",replacement="",x=cur[firstVectorColum])
		gene = cur[1,1]
		pathsV = strsplit(paths,split="; ")[[1]]
		outlist[[gene]] = pathsV
	}
	return(outlist)
}


#switchIds
#takes set of ids (vector)
#returns set of ids with as many switched as possible

switchIds<-function(idv,symtab=NULL){
	
	#see which can be switched
	#idv%in%symtab$
	dict = c("testSymbol")
	names(dict) = c("UniProt:P12956")
	for(i in 1:length(idv)){
		swap = dict[idv[i]]
		if(!is.na(swap)){
			idv[i] = swap
		}
	}
	return(idv)
}

switchPathSymbols<-function(po){
	#first switch the node names
	po@nodes = switchIds(idv=po@nodes)
	#second switch the names in the edges
	po@edges$src = switchIds(idv=po@edges$src)
	po@edges$dest = switchIds(idv=po@edges$dest)
	return(po)
}

#'@title Convert list of vectors to bipartate graph format. 
#'@description Used to restructure list of vectors to bipartate graph format. 
#'@param pth A list object, with each slot containing a vector. 
#'@param uni The full, unique set of elements in all the vectors in the list passed in the argument, pth. If not provided this will be automatically recalculated.
#'@return \code{matrix} object with logical values representing a bipartate graph; rows are set as the slot names of the input list. Columns are set as all the unique elements in the vectors the list contains. 
#'@export
list_to_table<-function(pth, uni=NULL){
	#switch from list pathway format to logic table pathway format
	#takes: pathways in list format
	#returns: pathways in logic table format
	if(is.null(uni)){uni = getUniverse(pth)}
	
	#initilize the bipartate graph table with all Fs (non membership)
	pmat = matrix(nrow=length(pth), ncol = length(uni), data = rep(x=F, times=length(pth)*length(uni)))  
	rownames(pmat) <- names(pth)
	colnames(pmat) <- uni
	
	for(i in 1:length(pth)){
		row = pth[[i]]
		#remove any missing data (ie, the NAs)
		row = row[which(!is.na(row))]
		#set appropriate membership of columns (gene names) in current path
		pmat[i,row] = T 
	}
	return(pmat)
}

#like matrix to list but does not produce empty lists
#takes a bipartite graph/logical matrix with row and column names the pathways and elements, respectively
#outputs a list, with each member a vector path members
matrix_to_list<-function(m){
	
	out = list()
	
	for(i in 1:nrow(m)){
		name = rownames(m)[i]
		members = as.character(colnames(m)[m[i,]])
		if(length(members)){
			out[[name]] = members
		}
	}
	return(out)
}


getUniverse<-function(gsv)
{
	####getUniverse
	#takes pathway list (list$pathname = [node1,node2,node3])
	#returns "universe", the list of unique nodes in all pathways
	num = length(names(gsv))
	uni = c()
	tot = 0
	for(i in 1:num){
		tot = tot + length(gsv[[i]])
		uni = union(uni, gsv[[i]])
	}
	#print(paste("sum of numbers of nodes in all pathways:", tot))
	#print(paste("total number of unique nodes:", length(uni)))
	return(uni)
}

test.getPaths<-function(){
	symtab = getHugoSymbols()
	path_file=NULL
	referenceFileName = system.file("extdata/pathMetaData.txt", package = "packageDir")
	testPths = getPaths(symtab=symtab, referenceFileName=referenceFileName)
	print(dim(testPths$paths))
	print(testPths$paths[1:2,1:2])
}

#getPaths
#main function to obtain pathways
#requires path meta data be available in file "./reference_data/paths/pathMetaData.txt"
# takes: path_file: the file name of the pathways to be loaded; if provided, interactive pathway selection will not be used
#				 referenceFileName: the name of the path meta data file
#				 symtab: the symbol look up table for hugo symbols 
# returns path_detail object

#'@title Get set of cellular pathways
#'@description Get \code{Path_Detail} object, housing gene sets for cellular pathways. 
#'@param path_file String giving path to file from which pathways should be loaded. Alternatively, a \code{Study} object can be provided here and pathways will automatically be retreived from it. 
#'@param referenceFileName Optional. String giving path to file where pathway meta data (date loaded, path repository name, etc..) can be found. 
#'@param symtab A gene identifier lookup table. If not provided, one will be automatically obtained (may take longer)
#'@param verbose Logical flag indicating if interactive user prompts and additional information should be displayed. 
#'@return \code{Path_Detail} object. 
#'@details 
#'Pathways are stored in the local working directory in the folder \code{./reference_data/paths/}. 
#'Pathways must be imported into this package so that meta-data and import records can be maintained. 
#'Path meta data is stored in the file \code{./reference_data/paths/pathMetaData.txt} . Any time a path
#'repository is imported, a record of its import is added to the path meta data file. When loading a path
#'repository for use in a study, this program checks the path meta data file, for available path repositories
#'and for path repository file locations. 
#'@export
#'@examples
#'
#'\dontrun{
#'#select and/or load pathways interactively
#'path_detail = getPaths()
#'}
#'
#'#Get pathways from a study object:
#'study=getTestStudyObject()
#'path_detail=getPaths(study)
#'
#'#select pathways from their file name:
#'path_detail = getPaths(path_file="./reference_data/paths/Reactome.2014.04.06.12.52.27.txt", 
#'											verbose=FALSE)
getPaths<-function(path_file=NULL, 
									 referenceFileName = "./reference_data/paths/pathMetaData.txt", 
									 symtab=NULL, 
									 verbose=T){
	
	if(class(path_file)=="Study") return(path_file@studyMetaData@paths)
	#initilize
	if(is.null(symtab)) symtab = getHugoSymbols(verbose=verbose)
	
	if(!length(path_file)) path_file=NULL
	
	if(!is.null(path_file)){
		if(!file.exists(path_file)){
			cat("\n\nThe path file\n",path_file,"\ncould not be found.\n\n")
			# 		blnk = readline("Press enter to continue and select a different pathway gene set file")
			path_file = NULL
		}
	}
	#provide user selections and path file selection
	while(T){#allows repeated loading of pathways 
		pRecord = choosePaths(ref=referenceFileName, path_file=path_file)
		if(!is.null(pRecord)) break
		#if excecution gets here, no pathway was selected and the option to import new pathway sets will be provided
		importPathways(symtab)
	}
	#load the actual paths
	paths = loadPaths(pRecord=pRecord, symtab=symtab)
	
	return(paths)
}


getPathMetaData<-function(ref="./reference_data/paths/pathMetaData.txt"){
	
	if(file.exists(ref)){
		fullTab = read.table(file=ref, header=T, sep="\t", stringsAsFactors=F, comment.char="")
	}else{
		
		fullTab = read.table(file=system.file("extdata/pathMetaData.txt",package = "packageDir"), 
												 header=T, sep="\t", stringsAsFactors=F, comment.char="")
		
		if(!file.exists(dirname(path=ref))) dir.create(path=dirname(path=ref), recursive=T, showWarnings=F)
		file.copy(from=system.file("extdata/pathMetaData.txt",package = "packageDir"), to=ref)
		file.copy(from=system.file("extdata/Reactome.2014.04.06.12.52.27.txt",package = "packageDir"), 
							to=paste0(dirname(ref), "/Reactome.2014.04.06.12.52.27.txt"))
		
	# write.table(x=fullTab, file=ref, quote=F, sep="\t", col.names=T, row.names=F)
		
		return(fullTab)
	}
}

#choosePaths()
#allows interactive selection of pathways
#takes: 	ref: path meta data file (default:"./reference_data/paths/pathMetaData.txt")
#					path_file: the file name of a pathway set (must be in GSEA format)
#returns: list: 1 row from the path meta data file given by the ref argument
choosePaths<-function(ref="./reference_data/paths/pathMetaData.txt", 
											path_file=NULL){
	
	if(!length(path_file)) path_file = NULL
	outRecord = NULL
	fullTab = getPathMetaData(ref=ref)
	if(is.null(path_file)){
		refTab=fullTab[,c(1:4)]
		colnames(refTab)<-gsub(pattern="[._]", replacement=" ",x=colnames(refTab))
		colnames(refTab)<-c("Repository source", "procurement date", "Number of paths", "Number of genes")
		refTab = cbind(refTab, 1:nrow(refTab))
		colnames(refTab)[ncol(refTab)]<-"Selection Number"
		
		cat("\nCurrently available pathway repositories:\n")
		print(refTab)
		while(T){
			line = readline(paste("Please enter the selection number for the pathway repository you would like to use.\n",
														"Or, if you would like to import a different pathway repository, (for example: a custom repository)\n",
														"please enter i: ", sep=""))
			if(line=="i"){
				return(NULL)
			}else if(!is.na(as.numeric(line))){
				if(as.numeric(line)%in%1:nrow(refTab)){
					outRecord = fullTab[as.numeric(line),]
					break
				}
			}
			print("Sorry, that input was not recognized, please try again.")
		}
	}else{
		#check if the path file can be found.. 
		outreci = which(fullTab$file==path_file)[1]
		if(is.na(outreci)) outreci = which(basename(fullTab$Original.file.name)==basename(path_file))[1]
		if(is.na(outreci)) outreci = which(basename(fullTab$file)==basename(path_file))[1]
		if(is.na(outreci)){
			warning(paste("Pathway file could not be found.... \nPath file:",path_file,"Files in pathway records:\n",paste(fullTab$file,sep="\n",collapse="\n")))
			return(choosePaths())
		} 
		outRecord = fullTab[outreci,]
	}
	return(outRecord)
}#choosePaths

getPathObject<-function(symtab, pRecord, pData, repset=NULL){
	#getPathObject
	#creates path list object with all possible slots filled
	# symtab: the hugo lookup table
	# pData: the bipartate graph of pathways
	# pRecord: the row of data from the path meta data file
	# repset: depricated: the set of repositories
# 	pd=list()
	pd = Path_Detail$new()
	pData = as.matrix(pData)
	colnames(pData)  = gsub(pattern="\\.",replacement="-",x=colnames(pData))#replace periods with dashes
	if(pRecord$Gene_Id_Type=="HUGO"){
		colnames(pData) = corsym(symbol_set=colnames(pData),verbose=F,symref=symtab)
	}else{
		cat("\nNon-HUGO symbols used in paths thus skipping gene symbol correction\n")
	}

	pd = setPathMetaData(pd=pd, 
											 symbol_type=pRecord$Gene_Id_Type,
											 p.file=pRecord$file, 
											 symtab=symtab, 
											 p.info = paste(pRecord$Set_name, pRecord$date_procured),
											 p.date = pRecord$date_procured,
											 p.source=pRecord$Set_name, 
											 p.paths=pData)
	pd$paths = pData
	return(pd)
}

setPathMetaData<-function(symtab, 
													p.paths, 
													p.file, 
													p.date, 
													p.source, 
													p.name=NULL, 
													p.info=NULL, 
													pd=list(), 
													symbol_type="HUGO"){
#setPathMetaData
	#establishes path meta data
	#takes: symtab: the symbol lookup table
	#				p.paths: the set of pathways in bipartate graph format
	#				p.file: the file containing the pathways
	#				p.date: the date the pathways were produced
	#				p.source: the source of the pathways
	#				p.name: optional, the name of the set of pathways 
	#								(if not provided, a formatted version of the pathname will be used)
	#				p.info: identifying data about the pathways
	#				pd: a nascent path object; providing this will preserve any custom list slots that have already been filled. 
	#returns: path_details object with meta data filled in
	p.info = ifelse(is.null(p.info), 
									yes=p.file, 
									no=p.info)
	p.name = ifelse(is.null(p.name), 
									yes=gsub(pattern=" [0-9]*.[0-9]*.[0-9]* [0-9]*.[0-9]*.[0-9]*.txt", replacement="", x=basename(path=p.file)), 
									no=p.name)
	
	p.paths = as.matrix(p.paths)
	
	pd[["name"]] = p.name
	pd[["file"]] = p.file
	pd[["info"]] = p.info
	pd[["date"]] = p.date
	pd[["source"]] = p.source
	pd[["gene_overlap_counts"]] = rep(T,nrow(p.paths))%*%p.paths
	pd[["full_path_length"]] = p.paths%*%rep(T,ncol(p.paths))
	pd[["symtable"]] = symtab
	pd[["paths"]] = p.paths
	pd[["symbol_type"]] = symbol_type
	pd[["graphite"]] = list()
	
 	return(pd)
}#setPathMetaData


#loadPaths
#used by getPaths()
#loads a set of paths from set of files at ./reference_data/paths/working_paths/
#will try to load path in GSEA format from the file denoted by file_name
#takes: paths pRecord: a list with one slot named "file" and other slots to be used by getPathObject
#returns: pathset
loadPaths<-function(pRecord, symtab=NULL){
	
	if(is.null(symtab)) symtab = getHugoSymbols()

	tmptab = list_to_table(pth=loadPathsAsSets(firstGeneColum=3, fname=pRecord$file))
	print("getting path object")
	pathset = getPathObject(symtab=symtab, pRecord=pRecord, pData=tmptab)
	print("got path object..")
	return(pathset)
}

#imports set of pathways in bipartate format
#takes: optional: fname: the file name of the bipartate path set
#													if fname is not provided, user will be prompted to select a set of pathways
#returns: bipartate graph of the path set
bipartateImport<-function(fname=NULL){
	preped_paths = list()
	cat("\nPlease select the tab-delimited file containing the pathway set you would like to use.\n")
	cat("Pathways in file should be in bipartate graph format, with the columns the\n")
	cat("colums the path members, the rows the path names and the values the values all\n")
	cat("T or F indicating path membership.\n")
	while(T){
		if(is.null(fname)) fname = file.choose()
		cat("\nFile selected: ", fname, "\n")
		preped_paths[["paths"]] = try(expr=read.table(file=fname, header=T, row.names=1, sep="\t"), silent=T)
		if(is.data.frame(preped_paths[["paths"]])) break
		if(!is.null(fname)) print("Alert! the file provided could not be read as a bipartate graph.")
		if(!is.null(fname)) return(NULL)
		blnk=readline("Sorry, the file chosen could not be read as a bipartate graph.\nPlease press enter and select another file.")
	}
	cat("\nFile created on", as.character(file.info(fname)$mtime),"\n")
	preped_paths[["paths"]] = read.table(file=fname, header=T, row.names=1, sep="\t")
	preped_paths[["paths"]] = as.matrix(preped_paths[["paths"]])
	
	preped_paths[["original_file_or_source"]] = fname
	preped_paths[["original_file_creation_date"]] = as.character(file.info(fname)$mtime)
	
	return(preped_paths)
}

#imports set of pathways in bipartate format
#takes: optional: fname: the file name of the GSEA format path set
#													if fname is not provided, user will be prompted to select a set of pathways
#returns: bipartate graph of the path set
GSEAImport<-function(fname=NULL){
	preped_paths = list()
	cat("\nPlease select the tab-delimited file containing the pathway set you would like to use.\n")
	cat("Pathways in file should be in GSEA format, with each row containing a path name followed by\n")
	cat("path members, all sepparated by tabs.\n")
	
	while(T){
		if(is.null(fname)) fname = file.choose()
		cat("\nFile selected: ", fname, "\n")
		preped_paths[["paths"]] =  try(expr=list_to_table(pth=loadPathsAsSets(firstGeneColum=3, fname=fname)),silent=T)
		if(is.data.frame(preped_paths[["paths"]])|is.matrix(preped_paths[["paths"]])) break
		if(!is.null(fname)) print("Alert! the file provided could not be read as a GSEA-format path set.")
		if(!is.null(fname)) return(NULL)
		blnk=readline("\nSorry, the file chosen could not be read as a GSEA-format path set.\n(Ex: <path_name><path_source><gene1><gene2>...)\nPlease press enter and select another file.")
	}
	
	cat("\nFile created on", as.character(file.info(fname)$mtime),"\n")
	preped_paths[["paths"]] = as.matrix(preped_paths[["paths"]])
	preped_paths[["original_file_or_source"]] = fname
	preped_paths[["original_file_creation_date"]] = as.character(file.info(fname)$mtime)
	
	return(preped_paths)
}#GSEAImport

#importPathways
#allows import of pathways either interactively or by choice and file name
#takes: symtab : table of symbol lookups
#				choice: currently b for bipartate graph or g for GSEA format
#				fname: the file name of the pathway repository
#								must match the format indicated by the choice argument
#returns: nothing
importPathways<-function(symtab=NULL, choice=NULL, fname=NULL){
	if(is.null(symtab)) symtab = getHugoSymbols()
	
	#initialize
	preped_paths = NULL#will hold the forming path object

	#allow selection
	while(T){
		if(!is.null(choice)) if(choice%in%c("b", "g")) break
		choice = readline("To provide a set of pathways in bipartate graph format, enter b\nTo provide a set of pathways in GSEA format, enter g: ")
	}
	
	#&open path file
	# conversion to bipartate format
	if(choice == "b"){
		preped_paths = bipartateImport(fname=fname)
	}else if(choice=="g"){
		preped_paths = GSEAImport(fname=fname)
	}
	
	usesHUGO = readline("Does the selected pathway set employ HGNC/HUGO symbols? (y/n)")
	
	if(usesHUGO=="y"){
		symbol_type="HUGO"
		#correct gene symbols
		colnames(preped_paths$paths)<-corsym(symbol_set=colnames(preped_paths$paths), 
																				 symref=symtab, 
																				 verbose=T)
	}else{
		symbol_type=readline("Please enter the type of gene identifiers employed (ex: UniProt) ")
		symtab=NULL
	}

	#setPathMetaData
	preped_paths = manualPathMetaData(preped_paths=preped_paths, symtab=symtab, symbol_type=symbol_type)

	#establish path set
	recordPathSet(ps=preped_paths)
	
}#importPathways

#manualPathMetaData
#allows user to interactively establish path meta data
#used by importPathways() function
#takes: preped_paths: nascent path_detail object with file name and paths slots filled in. 
#				symtab: the symbol lookup table
#returns: preped_paths: the path_detail object with meta data filled in
manualPathMetaData<-function(preped_paths, symtab=NULL, symbol_type="HUGO"){
	if(is.null(symtab)) symtab = getHugoSymbols()
	pathsFolder = "./reference_data/paths/"
	
# 	allSlots = c("paths","name", "file", "info", "date", "source", "gene_overlap_counts", "full_path_length", "symtable", "original_file_or_source", "original_file_creation_date")
	neededSlots = c("name", "date", "source")
	missingSlots = neededSlots[!neededSlots%in%names(preped_paths)]
	addedSlots = list()
	
	while(T&sum(!neededSlots%in%names(preped_paths))){
		defaultName = basename(path=preped_paths$original_file_or_source)
		defaultName = gsub(pattern=paste(".",file_ext(defaultName),"$",sep=""),replacement="",x=defaultName)
		defaultName = gsub(pattern="[._]", replacement=" ", x=defaultName)
		if("source"%in%missingSlots) addedSlots[["source"]] = readline(paste("Please enter the source of the pathway set.\nTo use the default,\"",defaultName,"\", just press enter: "))
		if("name"%in%missingSlots) addedSlots[["name"]] = readline(paste("Please enter the name of the set of pathways just added: \nTo use the default,\"",defaultName,"\", just press enter: "))
		if("date"%in%missingSlots) addedSlots[["date"]] = readline("Please enter the date the new pathway set was generated or obtained\n(from the providers of the pathway repository): ")
		for(s in names(addedSlots)) if(addedSlots[[s]]=="") addedSlots[[s]] = defaultName
		cat("\nPlease verify these entries are all correct: \n")
		for(r in names(addedSlots)) cat("For \"", r, "\" you entered: ", addedSlots[[r]], "\n",sep="")
		if(readline("If the above entries are all correct, press enter. \nTo change one or more, enter \"c\"")=="") break
	}

	preped_paths = c(preped_paths, addedSlots)
	fileName = paste(pathsFolder, preped_paths$name," ", gsub(pattern="[-:]", replacement=".", x=as.character(Sys.time())),".txt", sep="")
	cat("\nusing file name:", fileName,"to save the current set of pathways.\n")
	preped_paths = setPathMetaData(symtab=symtab, symbol_type=symbol_type,
																 pd=preped_paths,
																 p.info=paste(preped_paths[["name"]], preped_paths[["date"]]),
																 p.name=preped_paths[["name"]],
																 p.source=preped_paths[["source"]],
																 p.paths=preped_paths$paths, 
																 p.file=fileName, 
																 p.date=preped_paths[["date"]])
	return(preped_paths)
}#manualPathMetaData


#establishPathSet
# loads set of pathways into program memory
# a)gets meta data to put in path summary table
# b)corrects gene names
# c)makes GSEA format list of paths
# d)saves GSEA format paths file
# e)adds to path records file

test.recordPathSet<-function(){
	
	pathset=path_detail$paths
	
}

# 
# #obtain path meta data
# ps = setPathMetaData(symtab=hugo, 
# 										 p.date=p.date, 
# 										 p.source=p.source,
# 										 p.file=fileName, 
# 										 p.paths=pathset)
recordPathSet<-function(ps){
	pathset = ps$paths
	#make gsea data structure
	gsea = toGSEAformat(bip=pathset, psetName=ps$name)
	#save gsea data structure
	write.table(x=gsea, file=ps$file, quote=F, sep="\t", row.names=F,col.names=F)
	#add to pth records file
	prec = addToPathRecords(pfull=ps)
}#recordPathSet


#addToPathRecords
#used by recordPathSet() function to add path repository meta data to 
#takes: pfull: the completed path_detail object
#				precFname: the name of the path meta data/path records file
#returns: prec
addToPathRecords<-function(pfull, precFname = "./reference_data/paths/pathMetaData.txt"){
	#	data: path set name
	#				date
	#				number of paths
	#				number of genes
	# 			number of paths with zero nodes
	#				number of unapproved hugo symbols
	#				pathFileName
	#open file
	prec = NULL
	if(file.exists(precFname)) prec = read.table(file=precFname, header=T, sep="\t", quote="", stringsAsFactors=F)
	#assemble input data
	newRow = c(pfull$name, 
						 pfull$date, 
						 nrow(pfull$paths), 
						 ncol(pfull$paths), 
						 sum(pfull$full_path_length==0), 
						 sum(!colnames(pfull$paths)%in%pfull$symtable$Approved.Symbol), 
						 pfull$file, 
						 pfull$original_file_or_source, 
						 pfull$original_file_creation_date, 
						 pfull$symbol_type)
	prec = rbind.data.frame(prec, newRow)
	colnames(prec)<-c("Set_name", 
										"date_procured", 
										"Number_of_pathways", 
										"Number of genes in pathways", 
										"Number of paths with zero genes", 
										"Number of path nodes with out HUGO symbols", 
										"file", 
										"Original file name", 
										"Path_Import_Date", 
										"Gene_Id_Type")
	#check if more than one record use the same source path file
	prec=checkMultipleRecordsOneRepository(prec=prec)
	#write file
	write.table(x=prec, file=precFname, quote=F, sep="\t", col.names=T, row.names=F)
	return(prec)
}#addToPathRecords

#checkMultipleRecordsOneRepository()
#assures multiple records weren't produced from single repository files
#takes: the set of pathway records to be written to disk
#returns: the pathway records with any unwanted records removed
checkMultipleRecordsOneRepository<-function(prec){
	
	prec2=unique(prec)
	
	bnames = basename(prec2[["Original file name"]])
	duprows = c(which(duplicated(x=bnames)), which(duplicated(x=bnames, fromLast=T)))
	
	if(!length(duprows)) return(prec2) #return if there's nothing to correct

	#see if everything except the file name is the same
	dupTest = prec[,!names(prec)%in%"file"]
	if(nrow(dupTest)>nrow(unique(dupTest))){
		remFile = prec$file[duplicated(prec$file,fromLast=T)]
		prec2 = prec2[!duplicated(prec$file,fromLast=T),]
		file.remove(remFile)
		
	}

	reprows = prec2[duprows,]
	rownames(reprows)<-1:nrow(reprows)
	cat("\nIt appears the same set of pathways was imported multiple times.\n",
			"(repeated \"original file name\"s were found)\n",
			"These are the records that appear to be repeated:\n")
	print(reprows)
	while(T){
		line = readline("Please enter the line numbers of the path sets you would like to keep\n(sepparated by spaces): ")
		line = strsplit(x=line, split=" ")[[1]]
		line = try(as.numeric(line), silent=T)
		if(!sum(is.na(line))&!sum(!line%in%1:nrow(reprows))) break
		print("Error, input not recognized, please try again.")
	}
 
	deleteRows = duprows[!1:length(duprows)%in%line]
	#pull out the gsea file to be deleated
	todel = prec2[deleteRows,]
	#remove the GSEA file that was created
	file.remove(todel$file)
	# now limit prec2 to the selections
	prec3 = prec2[!1:nrow(prec2)%in%deleteRows,]
	
	return(prec3)
}#checkMultipleRecordsOneRepository()

test.toGSEAformat<-function(){
	
	pset = basename(STUDY@studyMetaData@paths$info)
	gsea = toGSEAformat(bip=STUDY@studyMetaData@paths$paths, psetName=pset)
}

toGSEAformat<-function(bip, psetName){
	out = matrix(data="", nrow = nrow(bip), ncol=3)
	for(r in 1:nrow(bip)){
		
		pname = rownames(bip)[r]
		genes = paste( colnames(bip)[bip[r,]] , collapse="\t")
		
		out[r,] = c(pname, psetName, genes)
		
	}
	return(out)
}






#loadPathsAsSets()
# loads paths from file in this format: <path name><\t><gene><space><gene>....
# everything after and including "firstGeneColum" is assumed to be a gene"
#the first column is assumed to be the path name or other unique identifier for the pathway
#
#takes: firstGeneColum: the first tab stop with gene names in it, usually 2 or 3
#				keepLowerCaseGenes:		keep genes who are all lower case
loadPathsAsSets<-function(firstGeneColum, 
													keepLowerCaseGenes=F, 
													fname="./reference_data/paths/ReactomePathways_dl_June_4_13.gmt.txt", 
													verbose=T){
	cat("\nLoading pathway sets..\n")
	#figure out what happens when this gets passed a bipartate graph.. 
	pset = read.table(fname, sep="\n",comment.char="", stringsAsFactors=F, quote="")
	outlist = list()
	excluded = list()
	
	for(i in 1:nrow(pset)){
		currow = strsplit(pset[i,],split="\t")[[1]]
		pname = currow[1]
		if(length(currow)>2){
			genes = currow[firstGeneColum:length(currow)]
			if(!keepLowerCaseGenes){
				#check gene names to see if any are all lower case
				badInds = !grepl(pattern="[A-Z]+",x=genes)
				if(sum(badInds)){
					print(genes)
					if(verbose) cat(i,"Lowercase gene names found. From",pname,"excluding",sum(badInds),"gene(s):", genes[which(badInds)], "indexes:",which(badInds),"\n")
					excluded[[pname]] = union(excluded, genes[badInds])
				} 
				genes = genes[!badInds]
			}
			outlist[[pname]] = genes
		}else{
			
			cat("\nThe pathway \n\"",pname,"\"\nwas found not to have any genes under the currently employed type of gene symbols.\n", sep="")
		}
		
	}
	if(length(excluded)){
		cat("\nExcluded gene symbols:\n")
		print(excluded)
	}

	return(outlist)
}#loadPathsAsSets


####getTargetMatrix
#get a bipartate graph indicating which genes are targets in which pathways
#paths with no targets are removed
#takes tgenes: vector of targeted genes
#      paths: set of pathways in logic matrix/bipartite graph format
#returns: bipartate graph: rows = paths columns = genes
getTargetMatrix<-function(tgenes, paths){
	cat("\ntargetMatrix(): ")
	nmat = paths #make copy of the paths
	
	tcols = colnames(paths)%in%tgenes#find how many of the path genes are in the set of target genes
	cat("genes from current set that are found the pathways: ", sum(tcols),"\n")
	if(sum(tcols)==0){
		print("There are no target genes in the set of paths")
		return(NULL)
	}else if(sum(tcols)==1){#if there is only one gene in the target genes that is in the set of paths
		print("There is only one target gene in the set of paths")
		#if there is only one gene, make a one column matrix out of it
		colreduced = nmat[,tcols]
		rowreduced = as.matrix(colreduced[colreduced],ncol=1)
		colnames(rowreduced) = colnames(paths)[tcols]
		return(rowreduced)
	}
	#  tcol = rep(T, nrow(nmat))%*%nmat>0
	nmat = nmat[,tcols,drop = FALSE]
	trow = nmat%*%rep(T, ncol(nmat))>0
	nmat = nmat[trow, ,drop = FALSE]
	
	return(nmat)
}

#appends column indicating if the gene is in the path repository
GenesInPaths<-function(ordgenelist, path_detail){
	print("inside GenesInPaths()")
	print(path_detail$info)
	pmemb = rep("No", times=nrow(ordgenelist))
	if(is.numeric(ordgenelist[,1])){
		pmemb[rownames(ordgenelist)%in%colnames(path_detail$paths)]="Yes"
		ordgenelist = cbind.data.frame(ordgenelist, pmemb)
		colnames(ordgenelist)<-c("Nubmer of patients","In the current set of pathways?")
		
		return(ordgenelist)
	}
	pmemb[ordgenelist[,1]%in%colnames(path_detail$paths)]="Yes"
	ordgenelist = cbind.data.frame(ordgenelist, pmemb)
	colnames(ordgenelist)<-c("Gene","Number of patients","In the current set of pathways?")
	return(ordgenelist)
}

getGenesFromPaths<-function(pids, STUDY){
	
	if(class(STUDY)=="Study") pths = STUDY@studyMetaData@paths$paths
	if(class(STUDY)=="Path_Detail") pths = STUDY$paths
	
	allGenes = c()
	for(p in pids){
		allGenes = union( colnames(pths)[pths[p,]], allGenes)
	}
	return(allGenes)
}

#getPathsForGenes()
#gives quick report on set of pathways an input set of genes belongs to
#takes: glist: a vector of genes symbols
#path_detail: a path_detail object
getPathsForGenes<-function(glist=NULL, path_detail=path_detail, verbose=F){
	
	if(is.null(glist)){
		glist = strsplit(readline("Enter set of genes sepparated by spaces: "),split=" ")[[1]]
	}
	tm = getTargetMatrix(tgenes=glist,paths=path_detail$paths)
	#add row and column totals
	
	gsums = rep(T,nrow(tm))%*%tm
	psums = tm%*%rep(T, ncol(tm))###note, this contains the grand total
	#reorder it all
	gorder = order(gsums,decreasing=T)
	porder = order(psums,decreasing=T)
	tm = tm[porder,gorder,drop=F]
	gsums = gsums[gorder]
	psums = psums[porder]
	
	tm = rbind(tm,gsums)
	tmppsums = c(psums,sum(psums))
	tm = cbind(tm,tmppsums)
	#correct marginal's names
	rownames(tm)[nrow(tm)] = "Marginal total: paths gene belongs to"
	colnames(tm)[ncol(tm)] = "Marginal total: input genes found in path"
	if(verbose){
		cat("\nGenes from input list not found in current pathway repository:\n")
		cat(setdiff(x=glist, colnames(path_detail$paths)))
		cat("\nSummary of genes from input list found in current pathway repository:\n")
		print(tm)

	}
	return(tm)
}#getPathsForGenes()

#getFrequencyMatrix
#takes:   genecounts: 2 column table <gene_names> <counts across cohort of patients inwhich that gene is active>
#         paths:      pathway target matrix in bipartate logic matrix form
#returns: frequencymatrix: pseudo bipartate graph: rows: path names, 
#																								columns: genes, 
#																								 values: number of times gene is found mutated across cohort
#getFrequencyMatrix(genecounts=genesum, paths=somtargetmatrix)
getFrequencyMatrix<-function(genecounts, paths){
	rownames(genecounts)<-genecounts[,1]	#make it like a dictionary so that gene extraction can be done
	inpath = genecounts[colnames(paths),]	#extract only the genes which are actually in the pathways
	out = NULL
	for(i in 1:nrow(paths)){#for each path
		nrow = inpath[,2]*paths[i,]#multiply every value in the i'th path by the count for that gene. . . NO!
		out = rbind(out, nrow)
	}
	rownames(out)<-rownames(paths)
	return(out)
}

#getMutCountMatrix: makes pseudo bipartate graph containing number of times gene is found mutated across cohort
#										as the values for each member in each pathway
#takes:   genecounts: 2 column table <gene_names> <counts across cohort of patients inwhich that gene is active>
#         paths:      pathway target matrix in bipartate logic matrix form
#returns: frequencymatrix: pseudo bipartate graph: rows: path names, 
#																								columns: genes, 
#																								 values: number of times gene is found mutated across cohort
#getFrequencyMatrix(genecounts=genesum, paths=somtargetmatrix)
getMutCountMatrix<-function(genecounts, paths){
	cat(" In getMutCountMatrix()")
	rownames(genecounts)<-genecounts[,1]	#make it like a dictionary so that gene extraction can be done
	inpath = genecounts[colnames(paths),]	#extract only the genes which are actually in the pathways
	out = matrix(data=0,ncol=ncol(paths),nrow=nrow(paths),dimnames=list(rn=rownames(paths),cn=colnames(paths)))
	
	for(gene in colnames(paths)){#for each path
		ncol = inpath[gene,2]*paths[,gene]#multiply every value in the i'th path by the count for that gene. . . NO!
		out[,gene] = ncol
	}
	rownames(out)<-rownames(paths)
	return(out)
}##getMutCountMatrix:



test.prepPathListForSave<-function(){
	matrixToSave = prepPathListForSave(gset=reactome)
	print(matrixToSave[1:2])
	
	write.table(x=matrixToSave, file="./reference_data/paths/GraphiteReactomeUniprotFormat.txt", 
							quote=F, 
							row.names=F, 
							col.names=F)
	
	pset = read.table("./reference_data/paths/GraphiteReactomeUniprotFormat.txt",
										header=F,
										sep="\n",
										comment.char="", 
										stringsAsFactors=F, 
										quote="")
	pset = matrix(data=pset[,1], ncol=1)
	checkEquals(target=matrixToSave, current=pset)
	GSEAImport(fname="./reference_data/paths/GraphiteReactomeUniprotFormat.txt")
}

#'@title prepPathListForSave()
#'@description preps paths in gene-set format (list with names=path names and values = vector of gene names) for saving as a tab delimited file
#'@param gset The list object of paths
#'@param path_source A string desribing the source of the pathways, to be pasted as the second element in each line of the output in the manner as "Reactome pathway" is seen here: 
#'       Abacavir metabolism	Reactome pathway	UniProt:P49902	UniProt:Q16774
#'@return one column matix, with each row containing the tab-delimited set: path name, path source and path members
#'@export
#'@import graphite
#'@examples 
#'library(graphite)
#'plistsAsDataFrame = prepPathListForSave(gset=reactome, path_source="Reactome pathway")
prepPathListForSave<-function(gset, path_source="Reactome pathway"){
	names(gset)
	pasted =sapply(X=gset, FUN=function(x){paste(nodes(x), sep="", collapse="\t")})
	dout = c()
	for(i in 1:length(pasted)){
		
		dout = c(dout,paste(names(pasted)[i], path_source, pasted[i], sep="\t"))
		
	}
	return(matrix(data=dout, ncol=1))
}#prepPathListForSave
