#integration test
cleanSlots<-function(tslots, cnamesToSkip = "(full_path_length)|(mgsa)"){
	
	
	if(is.list(tslots)&!is.data.frame(tslots)){
		for(cn in names(tslots)){
			print(cn)
			cur = tslots[[cn]]
			cur = cur[,!grepl(pattern=cnamesToSkip,
												x=colnames(cur), 
												ignore.case=T, 
												perl=T)]
			# 		if(class(cur)=="matrix") cur = as.data.frame(cur, stringsAsFactors=F)
			tslots[[cn]] = cur
		}
	}else{
		tslots = tslots[,!grepl(pattern=cnamesToSkip,
											x=colnames(tslots), 
											ignore.case=T, 
											perl=T)]
		colnames(tslots)<-gsub(pattern="[().]", replacement="_",colnames(tslots))
	}

	return(tslots)
}

test.abacavirMetabolismOverlapSets<-function(){
	
	require(RUnit)
	#first, run the study from the abacavir metabolism settings
	#load the data
	#------>the following two are both done automatically by the function
	#checkFileCopyDefault
	#--->move drug screen file from test data to working dir
	#--->move somatic mutations processing from test data to working dir. 
	
	# 	if(!file.exists("./output/")) dir.create(path=dirname(path="./output/"), recursive=T, showWarnings=F)
	# 	file.copy(from=system.file("extdata/pathMetaData.txt",package = "packageDir"), to=ref)
	
	# 	save(abacavirIntegrationTest, file="./output/abacavirIntegrationTestStudyObject.rda")
	# 	abacavirSettings = STUDY@studyMetaData@settings
	
	#		save(abacavirSettings, file="../packageDir/inst/extdata/abacavirSettings.rda")
	
	varname = load(system.file("extdata/abacavirSettings.rda",package="packageDir"), verbose=T)
	integrationSettings = get(varname[1])
	#initialize the study
	s1 = getStudyObject(path_detail=getDefaultPaths(), 
											settings=integrationSettings,
											study.name="integrationTestAvacavirMetabolism")
	s1 = loadBasicArms(STUDY=s1)
	#run the analysis from loaded settings
	s2 = autoRunFromSettings(study=s1)
	
	#test study slots
	#test somatic mutation slots
	varname2 = load(file=system.file("testData/intTest_subsomaticPathSum.rda",package = "packageDir"), verbose=T)
	res1 = get(varname2[1])
	cur1 = s2@results$somatic_mutation_aberration_summary$pathsummary
	cur1 = cur1[rownames(res1),]#make sure the rows are in the same order. 
	#remove the full path length column and normalize the column names
	cur1 = cleanSlots(cur1)
	res1 = cleanSlots(res1)
	checkEquals(target=res1, current=cur1)
	
	#get the current result
	cur2 = s2@results$functional_drug_screen_summary$pathsummary
	#get the gold standard
	var3 = load(file=system.file("testData/intTest_subfunPathSum.rda",package = "packageDir"), verbose=T)
	res2 = get(var3[1])
	#clean stuff up
	cur2 = cleanSlots(cur2)#remove the full path length column
	cur2 = cur2[rownames(res2),]#make sure the rows are in the same order. 
	#remove the full path length column and normalize the column names
	cur2 = cleanSlots(cur2)
	res2 = cleanSlots(res2)
	checkEquals(target=res2, current=cur2)
	
	#check the overlap analysis
	#get the current overlap analysis slots
	cur3 = s2@results$overlap_analysis
	nombres = c("Paths containing drug-sensitive genes",
							"Overlap between of genes in aberration enriched, not drug targeted paths",
							"Enriched for aberration and enriched for sensitive drug targets",
							"Drug targeted, not aberrationally enriched",
							"Aberration enriched, containing sensitive targets",
							"Aberration enriched, not drug targeted",
							"Aberrationally enriched, containing drug targets")
	cur3 = cur3[nombres]

	#get the saved overlap analysis slots
	var4 = load(file=system.file("testData/intTest_overlapAnalysisSlots.rda", package="packageDir"),verbose=T)
	res3 = get(var4[1])
	
	cur3 = cleanSlots(tslots=cur3, cnamesToSkip ="(full_path_length)|(mgsa)|(proportion_of_cohort)")
	res3 = cleanSlots(tslots=res3, cnamesToSkip ="(full_path_length)|(mgsa)|(proportion_of_cohort)")
	
	for(sl in names(res3) ){
		print(sl)
		cur = cur3[[sl]]
		ctarg = res3[[sl]]
		if(is.vector(cur)){
			if(length(cur)){
				cur = cur[order(cur)]
				ctarg = ctarg[order(ctarg)]
			}
		}else{
			if( nrow(cur)&ncol(cur) ){
				cur = cur[order(cur[,1]),]
				ctarg = ctarg[order(ctarg[,1]),]
				colnames(ctarg)<-colnames(cur)
				rownames(cur)<-NULL
				rownames(ctarg)<-NULL
			}
		}
		checkEquals(target=ctarg, current=cur)
		#	checkEquals.DataFrame(target=ctarg, current=cur)
	}
	
}