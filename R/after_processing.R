
library(XML)
library(rentrez)
library(RCurl)
library(GEOquery)
library(Biobase)
library(feather)
library("rols")
library(rjson)

######## datatable_to_use ########
geo=list.files("data",pattern = "GSE")
geo_id <- c()
for(i in 1:length(geo)){
	geo_id[i] <- rentrez::entrez_search(db="gds",term=geo)$ids[[1]]
}

GeoToSra <- function(ids) {
	fileURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=", paste(ids, collapse=","),sep="")
	geo_summary=xmlRoot(xmlTreeParse(getURL(fileURL)))

	GeoToSra_df <- data.frame(gse_id=character(), study=character(), sra_ftp=character(),taxon=character(),title=character(),
					summary=character(),stringsAsFactors=FALSE)
	for(i in 1: length(geo_summary)) {
		GeoToSra_df[i,] <- data.frame(gse_id=xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item[@Name='Accession']")[1][[1]]),
		study=xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item//Item//Item[@Name='TargetObject']")[1][[1]]),
		sra_ftp = xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item//Item//Item[@Name='TargetFTPLink']")[1][[1]]),
		taxon=xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item[@Name='taxon']")[1][[1]]), 
		title= xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item[@Name='title']")[1][[1]]),
		summary= xmlValue(getNodeSet(geo_summary[[i]], "//DocSum//Item[@Name='summary']")[1][[1]]),
		stringsAsFactors=FALSE)
	}
	return(GeoToSra_df)
}	
x=GeoToSra(geo_id)

datatable <- matrix(, nrow = length(geo), ncol = 6)
for ( i in seq_along(geo)) {
	print(i)
	load(paste0("data/",geo[i],"/eset.RData"))
	metadata <- pData(eset)
		
	datatable[i,] <- cbind(
	geo[i],
	length(unique(unlist(lapply(rownames(metadata), function(x) unlist(strsplit(x,"_"))[2])))),
	x[x$gse_id %in% geo[i], 4],
	x[x$gse_id %in% geo[i], 5], 
	x[x$gse_id %in% geo[i], 6],
	as.integer(nrow(metadata)))
}
datatable_to_use <- as.data.frame(datatable, stringsAsFactors=FALSE)
colnames(datatable_to_use) <- c("GEO accession", "Number of samples", "Species", "Title", "Study summary", "sra_runs")
write_feather(datatable_to_use, paste0("data/datatable_to_use.csv"))

######## GREIN datatable ########
datatable_to_use <- as.data.frame(read_feather(paste0("data/datatable_to_use.csv")),stringsAsFactors=F)
datafrm <- data.frame(datatable_to_use[,-6] ,check.names = FALSE)
datafr2 <- datafrm[order(datafrm$Species),]
datafr2 <- data.frame('GEO accession'= paste0('<button id="blah" type="button" class="btn btn-primary action-button" value=', datafr2[,1], '>', datafr2[,1], '</button>'),
				datafr2[,-1], stringsAsFactors=F, check.names=F)

write_feather(as.data.frame(datafr2), paste0("data/GRIEN_datatable.csv"))

######## all_samples_feather ########
dir <- list.files("data",pattern = "GSE")
if(file.exists("data/all_samples_feather.csv")){
	all_samples <- as.data.frame(read_feather(paste0("data/all_samples_feather.csv")),stringsAsFactors=F)
} else {
	all_samples <- NULL
}
new_samples <- dir[which(!dir %in% unique(all_samples[,"gse_acc"]))]
datalist = list()
for (i in 1:length(new_samples)) {
	print(i)
	gse_acc <- new_samples[i]
	print(gse_acc)
	#setwd(paste0(destdir,"/",gse_acc))
	load(paste0("data/",gse_acc,"/eset.RData"))
	meta <- pData(eset)
	if("geo_accession" %in% names(meta)){
		meta2 <- meta[,c("geo_accession","Run","Sample")]
	} else {
		meta$geo_accession <- NA
		meta2 <- meta[,c("geo_accession","Run","Sample")]
	}
	meta2$gse_acc=rep(gse_acc,nrow(meta2))
	rownames(meta2) <- NULL
	datalist[[i]] <- meta2 
}
all_samples = rbind(all_samples, do.call(rbind,datalist))
write_feather(all_samples, paste0("data/all_samples_feather.csv"))

######## samples_ontology_full_feather ########
all_samples <- as.data.frame(read_feather(paste0("data/all_samples_feather.csv")),stringsAsFactors=F, check.names=F)
metaSRA=fromJSON(file="data/metasra.v1-4.json")
x=metaSRA
common_metasra <- x[which(names(x) %in% all_samples[,"Sample"])]
if(file.exists("data/samples_ontology_full_feather.csv")){
	datatable <- as.data.frame(read_feather(paste0("data/samples_ontology_full_feather.csv")),stringsAsFactors=F)
	common_metasra_new <- common_metasra[which(!names(common_metasra) %in% datatable[,"Sample"])]
	datatable_new <- matrix(, nrow = length(common_metasra_new), ncol = 3)
	ont_terms <- olsNamespace(Ontologies())
	for(i in 1:length(common_metasra_new)){
		print(i)	
		mot <- common_metasra_new[[i]]$"mapped ontology terms"
		if(length(mot)>0){
			tm <- as.character(sapply(mot, function(x) unlist(strsplit(x,":",""))[1]))
			avail <- data.frame(mot=mot, len=sapply(tm, function(x) length(grep(x, ont_terms, ignore.case=T))))
			avail_term <- as.character(avail[which(avail[,2]>0),1])
			common_metasra_new[[i]]$"mapped ontology terms" <-  sapply(avail_term, function(x) termLabel(term(Ontology(unlist(strsplit(x, ":",""))[1]),x)))
			datatable_new[i,] <- cbind(names(common_metasra_new[i]), paste0(paste0(names(common_metasra_new[[i]][[1]]),":", common_metasra_new[[i]][[1]]), collapse=" ; "),
								 paste0(common_metasra_new[[i]][[3]], collapse=",")
								)
		} else {
			datatable_new[i,] <- cbind(names(common_metasra_new[i]), NA, paste0(common_metasra_new[[i]][[3]], collapse=","))	
		}
	}
	df <- datatable_new
	df <- na.omit(df)
	df_new <- data.frame(df, stringsAsFactors = F)
	colnames(df_new) <- c("Sample", "mapped_ontology_terms","Sample_type")
	write_feather(df_new, paste0("/opt/raid10/genomics/naim/Thesis/geo/samples_ontology_full_feather.csv"))

} else {
	datatable <- NULL
	common_metasra_new <- common_metasra
	datatable_new <- matrix(, nrow = length(common_metasra_new), ncol = 3)
	ont_terms <- olsNamespace(Ontologies())
	for(i in 1:length(common_metasra_new)){
		print(i)	
		mot <- common_metasra_new[[i]]$"mapped ontology terms"
		if(length(mot)>0){
			tm <- as.character(sapply(mot, function(x) unlist(strsplit(x,":",""))[1]))
			avail <- data.frame(mot=mot, len=sapply(tm, function(x) length(grep(x, ont_terms, ignore.case=T))))
			avail_term <- as.character(avail[which(avail[,2]>0),1])
			common_metasra_new[[i]]$"mapped ontology terms" <-  sapply(avail_term, function(x) termLabel(term(Ontology(unlist(strsplit(x, ":",""))[1]),x)))
			datatable_new[i,] <- cbind(names(common_metasra_new[i]), paste0(paste0(names(common_metasra_new[[i]][[1]]),":", common_metasra_new[[i]][[1]]), collapse=" ; "),
								 paste0(common_metasra_new[[i]][[3]], collapse=",")
								)
		} else {
			datatable_new[i,] <- cbind(names(common_metasra_new[i]), NA, paste0(common_metasra_new[[i]][[3]], collapse=","))	
		}
	}
	df <- rbind(as.matrix(datatable), datatable_new)
	df <- na.omit(df)
	df_new <- data.frame(df, stringsAsFactors = F)
	write_feather(df_new, paste0("/opt/raid10/genomics/naim/Thesis/geo/samples_ontology_full_feather.csv"))
}
