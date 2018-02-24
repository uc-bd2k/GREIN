
create_eset <- function(gse_acc,species, destdir) {

	# Generate eset.
	#
	# Args:
	#   gse_acc: GEO accession number.
	#	species: Either Homo sapiens, Mus musculus, or Rattus norvegicus.
	#   destdir: Directory to save the results.
	#
	# Returns:
	#   eset.

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
	setwd(paste0(destdir,"/",gse_acc))
	countsTable <- read.csv(file="countsTable.csv", header=TRUE, row.names=1)
	metadata <- read.csv(file="metadata.csv", header=TRUE, row.names=1)
	annot <- AnnotationDbi::select(gene_ensembl(species),keys=rownames(countsTable),columns=c("SYMBOL","SYMBOL", "GENENAME"),keytype="ENSEMBL")
	annot <- annot[!duplicated(annot[,1]),]
	rownames(annot) <- annot[,1]
	fdata <- annot[rownames(countsTable),]
	
	phenodata <- new("AnnotatedDataFrame",data=metadata)
	featuredata <- new("AnnotatedDataFrame",data=fdata)
	eset <- ExpressionSet(assayData=data.matrix(countsTable),phenoData=phenodata, featureData=featuredata)
	cat("saving eset\n")
	save(eset, file=paste(destdir,"/", gse_acc, "/","eset.RData", sep=""))
}