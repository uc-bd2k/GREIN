
### Run this to process a GEO RNA-seq data
#Install the packages if not already installed
install.packages(c("devtools", "XML", "parallel", "utils", "rentrez", "RCurl", "rjson", "feather")
source("https://bioconductor.org/biocLite.R")
biocLite(c("GEOquery", "Biobase", "rols", "tximport", "EnsDb.Hsapiens.v86", "EnsDb.Rnorvegicus.v79", "EnsDb.Mmusculus.v79",
"AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))

library(XML)
library(rentrez)
library(RCurl)
library(GEOquery)
library(Biobase)
library(devtools)
library(parallel)
library(tximport)	
library(EnsDb.Hsapiens.v86)
library(EnsDb.Rnorvegicus.v79)
library(EnsDb.Mmusculus.v79)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(rols)
library(rjson)
library(feather)

install_github("uc-bd2k/GREP2")
library(GREP2)

logdir <- "data/user_geo_request/"
destdir <- "data/user_geo_request/"

cat(paste("STEP 1: Processing starts... ","\n",sep=""),file=paste0(logdir,"/",geo_series_acc,"/log.txt"))	
## Change the parameters accordingly
process_geo_rnaseq (geo_series_acc=geo_series_acc,destdir="data/user_geo_request",
	ascp=TRUE,prefetch_workspace="path_to_prefetch_workspace",
	ascp_path="path_to_aspera",get_sra_file=FALSE,trim_fastq=FALSE,
	trimmomatic_path=NULL,index_dir="path_to_indexDir",
	species=species,countsFromAbundance="lengthScaledTPM",n_thread=2)
	
load(paste0("data/",geo_series_acc,"/counts_data_list.RData", verbose=T)
load(paste0("data/",geo_series_acc,"/metadata.RData", verbose=T)

gene_ensembl= function(species) {
	if (species == "Homo sapiens") {
		return(org.Hs.eg.db)
	}
	else if (species == "Rattus norvegicus") {
		return(org.Rn.eg.db)
	}
	else if (species == "Mus musculus") {
		return(org.Mm.eg.db)
	}
	else {
		return(NULL)
	}
}
countsTable <- counts_data_list$gene_counts
metadata <- metadata
annot <- AnnotationDbi::select(gene_ensembl(species),keys=rownames(countsTable),columns=c("SYMBOL","SYMBOL", "GENENAME"),keytype="ENSEMBL")
annot <- annot[!duplicated(annot[,1]),]
rownames(annot) <- annot[,1]
fdata <- annot[rownames(countsTable),]

phenodata <- new("AnnotatedDataFrame",data=metadata)
featuredata <- new("AnnotatedDataFrame",data=fdata)
eset <- ExpressionSet(assayData=data.matrix(countsTable),phenoData=phenodata, featureData=featuredata)
save(eset, file=paste("data/",geo_series_acc,"/eset.RData", sep=""), compress=F)

closeAllConnections()
