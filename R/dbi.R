#dbi

dbiWrite<-function(tabla, name="test_table", append=F, dbname="interactome", port="5432", clean = T ){
  if(name =="test_table")
  {
    print("warning, writing table data to test_table")
  }
  
  if(clean) colnames(tabla) = gsub(pattern="[.]", replacement="_", x=colnames(tabla))
  library("DBI")
  #library("RMySQL")
  #drv <- dbDriver("MySQL")
  #con <- dbConnect(drv, dbname="HNSCCdb", port=3306, user="root")
  
  library("RPostgreSQL")
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname=dbname, password="soyyo", port=port, user="postgres")
  
  #write the table
  dbWriteTable(conn=con, name=name, row.names = 0, value=tabla, overwrite=!append, append=append)
  print(dbExistsTable(con, name=name))
  dbDisconnect(con)
}

# addForeignKey<-function(tableName, dbName, keyColumns=c()){
# 	
# }
	

dbiRead<-function(query,dbname="metagenomics"){
	library("RPostgreSQL")
	drv <- dbDriver("PostgreSQL")
	con <- dbConnect(drv, 
									 dbname=dbname, 
									 password="soyyo", 
									 port=5432, 
									 user="postgres")
	rs <- dbSendQuery(con,query)
	resset = fetch(rs,n=-1)
	dbDisconnect(con)
	return(resset)
} #database interface, returns table

dbCall<-function(sql, dbname="interactome", port="5432"){
	library("DBI")
	#library("RMySQL")
	#drv <- dbDriver("MySQL")
	#con <- dbConnect(drv, dbname="HNSCCdb", port=3306, user="root")
	library("RPostgreSQL")
	drv <- dbDriver("PostgreSQL")
	con <- dbConnect(drv, dbname=dbname, password="soyyo", port=port, user="postgres")
	dbSendQuery(conn=con, statement=sql)
	dbDisconnect(con)
}



test.FileToDb<-function(){
	fname="./drugDB/drugbank/all_target_ids_all.csv"
	tableName="targets"
	sep=","
	header=T
	quote="\""
	dbname = "drugs"
	port=5432
	app=F
	delim=","
	fileToDb(tableName=tableName, app=app, dbname=dbname, fname=fname, sep=sep, header=header, quote=quote, port=port, delim=delim)
}

fileToDb<-function(tableName, dbname,  fname, app=F, sep=NULL, header=NULL, quote="", port=5432, delim=","){

# 	tab = read.table(file=fname, sep=",", header=T, quote="\"")
	#first open the file
	tab = read.table(file=fname, sep=sep, header=header, quote=quote)
	#now make a shell out of it: 
	tshell = tab[c(),]
	dbiWrite(tabla=tab, name=tableName, dbname=dbname, port=port, append=F)
	
	readline("There is an issue where double quotes in the palce of a blank cell will screw up the copy to db.")
	#now copy from the file to the table
	query = paste("copy ", tableName, " from '", fname, "' DELIMITER '", delim, "' NULL '\"\"' CSV", sep="")
	#COPY genes FROM '/Users/samhiggins2001_worldperks/tprog/main_131219/drugDB/drugbank/all_target_ids_all.csv' DELIMITER ',' CSV;
	dbCall(sql=query, dbname=dbname, port=port)
	
}

fname="./drugDB/drugbank/all_target_ids_all.csv"
tab = read.table(file=fname, sep=sep, header=header, quote=quote, stringsAsFactor=F)


#findHGNC
#processAllTargetIds
#first, select only the human genes: 

htab = tab[tab$Species == "Homo sapiens",]
hugo = STUDY@studyMetaData@paths$symtable
# > dim(tab)
# [1] 4141   15
# > dim(htab)
# [1] 2106   15

#next, attempt-correction of gene Name column 
tmpName = corsym(symbol_set=htab$Gene.Name,verbose=F, symref=hugo)
htab$Gene.Name = tmpName
################################# check

notHugoIndex = which(!tmpName%in%hugo$Approved.Symbol)
notHugo= htab[notHugoIndex,]

hextract = hugo[hugo$HGNC.ID%in%notHugo$HGNC.ID,]
rownames(hextract)<-hextract$HGNC.ID

#now obtain the hugo symbols by their row names
htab$Gene.Name[notHugoIndex] = hextract[notHugo$HGNC.ID,]$Approved.Symbol

sum(htab$Gene.Name=="")

htab$Gene.Name[is.na(htab$Gene.Name)] = ""

#gene names are now fixed

#extract the GSEA format lines now: 
targetTable = htab[,c("Gene.Name", "Drug.IDs")]

makeStacked<-function(dfin){
	
	ltmp = list()
	#first make a list out of it
	for(i in 1:nrow(dfin)){
		ltmp[[dfin[i,1]]] = strsplit(x=dfin[i,2], split="; ")[[1]]
	}
	tabd = list_to_table(pth=ltmp)
	
	#now, make the output data frame and fill it
	dfout = data.frame(matrix(ncol=2,nrow=sum(tabd), dimnames=list(1:sum(tabd), c("geneID","drugID"))))
	
	curi = 1
	for(i in 1:nrow(tabd)){
		if(sum(tabd[i,])){
			lindex = curi + sum(tabd[i,]) - 1
			if(lindex == curi){
				dfout[curi,1] = rownames(tabd)[i]
				dfout[curi,2] = colnames(tabd)[tabd[i,]]
			}else{
				dfout[curi:lindex,1] = rownames(tabd)[i]
				dfout[curi:lindex,2] = colnames(tabd)[tabd[i,]]
			}
			curi=lindex + 1
		}
	}
	return(dfout)
}


targetTable = makeStacked(dfin=targetTable)

dbiWrite(tabla=targetTable, name="targets", append=F, dbname="drugs")

drugFname = "./drugDB/drugbank/drug_links.csv"
drugTab = read.table(drugFname, sep=",", header=T)

dbiWrite(tabla=drugTab, name="drug", append=F, dbname="drugs")
#what am I doing?
#just got the targeting table writen, 
#now check that the gene and drug tables are being made

#now load the

#load data into pisces database
#source('./dbi.R')


# dbiWrite(tabla=results$somatic_mutation_aberration_summary$unfiltered_data,
# 				 name="somatic_mutation_data", 
# 				 dbname="pisces")




