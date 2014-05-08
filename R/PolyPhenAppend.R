#deal with OHSU PolyPhen results
# source('./PolyPhenProcessing_functions.R')


preProcessPdat<-function(fname){
	cat('\nPreprocessing PolyPhen output... ')
	#removes/corrects unwanted formatting in initial PolyPhen output and prepares file for regular input
	prepdat = read.table(file=fname, header=F,sep="\n",comment.char="", stringsAsFactors=F)
	#first add an extra column header to the first row
	if(!length(grep(pattern="ref.col",prepdat[1,]))){
		prepdat[1,] = paste(prepdat[1,], "pp.ref.dat","map.col", sep="\t")
	}
	postind = grepl(pattern="^## ",x=prepdat[,1])
	prepdat=prepdat[!postind,,drop=F]
	newFname = paste(fname, ".reformatted.txt",sep="")
	write.table(x=prepdat, file=newFname, sep='\n',row.names=F, col.names=F, quote=F)
	cat('Preprocessing complete\n')
	return(newFname)
}

# AddPolyPhenToVariantTypes
# unique(ohsuseq$PolyPhen)
loadPolyPhenResults<-function(fname){
	cat(" Loading PolyPhen results from file\n",fname,"\n......")
	newFname = preProcessPdat(fname)
	pdat = read.table(file=newFname, 
										sep="\t",
										header=T,
										comment.char="", 
										stringsAsFactors=F)
	#first make a map col: 
# 	mc = paste(indexCol,"|",
# 									 datf1$Hugo_Symbol, "|",
# 									 datf1$Start_Position, "|",
# 									 datf1$Reference_Allele,
# 									 "/",
# 									 datf1$Tumor_Seq_Allele1,
# 									 sep="")
	#parse out the map.col
	mc = pdat$map.col
	
	#get new map column
	mcs = data.frame(matrix(byrow=T,data=unlist(strsplit(mc, split="\\|")), ncol=4), stringsAsFactors=F)
	colnames(mcs)<-c("index", "Symbol", "Start.Pos", "pid")
	pdat = cbind.data.frame(pdat, mcs, stringsAsFactors=F)
	#add keys
	rkeys = rep("", times=nrow(pdat))
	for(i in 1:nrow(pdat)){
		rkeys[i] = paste(pdat[i,c("pid","Start.Pos","Symbol")], collapse="")
	}
	pdat = cbind.data.frame(pdat,rkeys, stringsAsFactors=F)
	#fix the white space in prediction
	pdat$prediction = gsub(pattern="^ *", replacement="", x=pdat$prediction)
	cat("results loaded..\n ")
	return(pdat)
}

# getMapColumn(pdat)<-function(){
# 	
# 	mc = cbind.data.frame(pdat$)
# 	
# }


getNumericScoreColumn<-function(xcol){
	cat("recoding scores as numeric...")
	#add column indicating numeric score for each polyphen score
	ppmn = c("benign","probably damaging","possibly damaging","unknown", "\\N")
	ppsc = c(1,3,2,0,-1)
	names(ppsc)<-ppmn
	nscol = rep(0, nrow(xcol))
	for(n in ppmn){
		nscol[grep(pattern=n, x=xcol$prediction)] = ppsc[n]
	}
	cat("scores recoded\n")
	return(nscol)
}

PPMultiMapRedux<-function(pdat){
	cat("\nTaking max polyphen score for each gene...\n")
	#takes the formatted polyphen output and reduces the rows to only include
	# the max score for each index
	#return: data frame, same structure as input; same data, but rows removed; no multi mapping to indexes
	#extract the needed columns
	
	xcol = pdat[,c("index","pid","Start.Pos","Symbol","prediction","rkeys")]
	
	numscore = getNumericScoreColumn(xcol = xcol)
	xcol = cbind(xcol, numscore)
	#find the unique set of indexes
	ui = unique(xcol$rkeys)
	gindexes = c()
	for(i in ui){# for each index, find the max polyphen score
		seti = which(xcol$rkeys == i)
		sel = which(xcol$numscore[seti] == max(xcol$numscore[seti]))[1]
		gindexes = c(gindexes, seti[sel])
	}
	redux = xcol[gindexes,]	
	cat("\nMax polyphen score aquired for each gene...\n")
	return(redux)
}

# scoresout = splitScoresOut(seqdat=ohsuseq)

appendPPScoresOHSU<-function(oseqdat){
	li = oseqdat$PolyPhen!=""
	npol = paste("PolyPhen_", oseqdat$PolyPhen[li], sep="")
	nvarcol = paste(oseqdat$Consequence[li], npol, sep=";")
	oseqdat$Consequence[li] = nvarcol
	return(oseqdat)
}

#deal with TCGA polyphen results

#open the polyphen results and use loadPolyPhenResults to break out the reference data

getKeyedPredictionList<-function(fname){
	cat("getting PolyPhen predictions.. \n")
	pdat = loadPolyPhenResults(fname=fname)
	fpdat = PPMultiMapRedux(pdat = pdat)
	print(colnames(fpdat))
	keyedPredictions = fpdat$prediction
	rkeys = rep("", times=nrow(fpdat))
	for(i in 1:nrow(fpdat)){
		rkeys[i] = paste(fpdat[i,c("pid","Start.Pos","Symbol")], collapse="")
	}
	names(keyedPredictions)<-rkeys
	cat('... predictions aquired\n')
	return(keyedPredictions)
}

#make keys for maf data
MafDataKeys<-function(mafData){
	cat("getting keys from .maf data... ")
	res = rep("", times=nrow(mafData))
	for(i in 1:nrow(mafData)){
		res[i] = paste(mafData[i,c("pid","Start_Position", "Hugo_Symbol")], collapse="")
	}
	cat("keys aquired\n")
	return(res)
}

# kpd[c("","TCGA.AB.2934118536561ABCG4","TCGA.AB.2934118536561ABCG4")]

AlterVariantClassification<-function(kpd, mdk, mafData){
	cat("Altering variant classifications in .maf data...")
	preds = kpd[mdk]
	preds_mod = paste("_PolyPhen_", preds, sep="")
	preds_mod[is.na(preds)] = ""
	names(preds_mod) = NULL
	mafData = cbind(mafData, preds_mod)
	mafData$Variant_Classification = paste(mafData$Variant_Classification,preds_mod, sep="")
	cat("classifications altered\n")
	return(mafData)
}

addPolyPhenResults<-function(mafData, tracker, s){
	
	s=setting(s=s, prompt="Would you like to include PolyPhen analysis results in this analysis? (y/n) ")
	uin=s$.text
	if(uin=="y"){
		s=setting(s=s,prompt="\nPlease select PolyPhen output file\n")
		fname = s$.text		
		tracker[["PolyPhen results file used"]] = fname
		kpd = getKeyedPredictionList(fname=fname)
		mdk = MafDataKeys(mafData=mafData)
		mafData = AlterVariantClassification(kpd=kpd, mdk=mdk, mafData=mafData)

	}
	return(list(mafData=mafData, tracker=tracker, s=s))
}


