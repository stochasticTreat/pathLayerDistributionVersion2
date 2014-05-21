#deal with OHSU PolyPhen results
# source('./PolyPhenProcessing_functions.R')


#loads poly phen results
#-runs preprocessing function() to make the table loadable
#-extracts the pid/start.pos/symbol from the map.col
#-cleans up white space from the predictions
#returns the cleaned, preped table of polyphen predictions
loadPolyPhenResults<-function(fname){
	cat(" Loading PolyPhen results from file\n",fname,"\n......")
	newFname = preProcessPdat(fname=fname) #clean up the polyphen output so that the table can be correctly read in as a table with the correct nubmer of columns
	pdat = read.table(file=newFname, 
										sep="\t",
										header=T,
										comment.char="", 
										stringsAsFactors=F)

	mc = pdat$map.col
	
	#get new map column
	mcs = data.frame(matrix(byrow=T,data=unlist(strsplit(mc, split="\\|")), ncol=4), stringsAsFactors=F)
	colnames(mcs)<-c("index", "Symbol", "Start.Pos", "pid")
	pdat = cbind.data.frame(pdat, mcs, stringsAsFactors=F)
	#add keys
	rkeys = rep("", times=nrow(pdat))
	
	pdatex = pdat[,c("pid","Start.Pos","Symbol")]
	rkeys = apply(X=pdatex, MARGIN=1, FUN=function(x){paste(x, collapse="")})
	
	pdat = cbind.data.frame(pdat,rkeys, stringsAsFactors=F)
	#fix the white space in prediction
	pdat$prediction = gsub(pattern="^ *", replacement="", x=pdat$prediction)
	cat("results loaded..\n ")
	return(pdat)
}#loadPolyPhenResults

PPMultiMapRedux<-function(pdat){
	cat("\nFinding max polyphen score for each gene...\n")
	#takes the formatted polyphen output and reduces the rows to only include
	# the max score for each index
	#return: data frame, same structure as input; same data, but rows removed; no multi mapping to indexes
	#extract the needed columns
	
	xcol = pdat[,c("index","pid","Start.Pos","Symbol","prediction","rkeys")]
	
	numscore = getNumericScoreColumn(xcol = xcol, pdat=pdat)
	xcol = cbind(xcol, numscore)
	#find the unique set of indexes
	
	agres = aggregate(x=numscore, by=list(ukeys=xcol$rkeys), FUN=max)
	colnames(agres)[2]<-"numscore"
	
	scoreDict = c("benign","probably damaging","possibly damaging","unknown", "\\N")
	names(scoreDict)<-as.character(c(1,3,2,0,-1))
	
	predictions = scoreDict[as.character(agres$numscore)]
	pout = cbind.data.frame(rkeys=agres$ukeys, predictions = predictions, stringsAsFactors=F)
	
# 	
# 	ui = unique(xcol$rkeys)
# 	gindexes = c()
# 	for(i in ui){# for each index, find the max polyphen score
# 		seti = which(xcol$rkeys == i)
# 		sel = which(xcol$numscore[seti] == max(xcol$numscore[seti]))[1]
# 		gindexes = c(gindexes, seti[sel])
# 	}
# 	redux = xcol[gindexes,]
# 	

	cat("\nMax polyphen score aquired for each gene...\n")
	return(pout)
}


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
	cat("\ngetting PolyPhen predictions.. \n")
	pdat = loadPolyPhenResults(fname=fname)
	fpdat = PPMultiMapRedux(pdat = pdat)
	cat("\nMultimapping delt with\n")
	print(colnames(fpdat))
	# 	keyedPredictions = fpdat$prediction
	# 	rkeys = rep("", times=nrow(fpdat))
	# 	for(i in 1:nrow(fpdat)){
	# 		rkeys[i] = paste(fpdat[i,c("pid","Start.Pos","Symbol")], collapse="")
	# 	}
	# 	names(keyedPredictions)<-rkeys
	# 	
	# 	
	kpd2 = fpdat$predictions
	names(kpd2)<-fpdat$rkeys
	cat('... predictions aquired\n')
	return(kpd2)
}

#make keys for maf data
MafDataKeys<-function(mafData){
	cat("getting keys from .maf data... ")

	mres = mapply(FUN=function(a, b, c){paste(a,b,c,collapse="",sep="")}, mafData$pid, mafData$Start_Position, mafData$Hugo_Symbol)

	# 	res = rep("", times=nrow(mafData))
	# 	for(i in 1:nrow(mafData)){
	# 		res[i] = paste(mafData[i,c("pid","Start_Position", "Hugo_Symbol")], collapse="")
	# 	}
	
	cat("keys aquired\n")
	return(mres)
}

# kpd[c("","TCGA.AB.2934118536561ABCG4","TCGA.AB.2934118536561ABCG4")]

AlterVariantClassification<-function(kpd, mdk, mafData){
	cat("Altering variant classifications in .maf data...")
	preds = kpd[mdk]
	preds_mod = paste("_PolyPhen_", preds, sep="")
	preds_mod[is.na(preds)] = ""
	names(preds_mod) = NULL
	mafData$Variant_Classification = paste(mafData$Variant_Classification, preds_mod, sep="")
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


