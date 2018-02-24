
tximport2counts <- function(gse_acc,species,destdir){

	## Transcripts to gene counts
	setwd(paste0(destdir,"/",gse_acc))
	metadata <- read.csv(file="metadata.csv", header=TRUE, row.names=1)

	edb= function(species) {
		if (species == "Homo sapiens") {
			get(load("/data/Tx_ensemble_hs.RData"))
		}
		else if (species == "Rattus norvegicus") {
			get(load("/data/Tx_ensemble_rn.RData"))
		}
		else if (species == "Mus musculus") {
			get(load("/data/Tx_ensemble_mm.RData"))
		}
		else {
			return(NULL)
		}
	}
	assign('Tx.ensemble', edb(species))
	Tx.ensemble <- get('Tx.ensemble')	
	tx2gene<- Tx.ensemble[,c(1,2)]
	all_files <- file.exists(paste(destdir, "/",gse_acc, "/salmon/", rownames(metadata),"_transcripts_quant/",rownames(metadata), "_quant_new.sf", sep=""))
	files <- file.path(paste(destdir, "/",gse_acc, "/salmon/", rownames(metadata)[all_files],"_transcripts_quant/",rownames(metadata)[all_files], "_quant_new.sf", sep=""))

	if(length(files) == dim(metadata)[1]){
	
		txi <- tximport(files, type = "salmon", tx2gene = tx2gene, importer = read.delim, dropInfReps=TRUE, countsFromAbundance = "lengthScaledTPM")
		save(txi, file=paste(destdir,"/", gse_acc, "/salmon/","txi.RData", sep=""))
		countsTable <- txi$counts
		colnames(countsTable) <- rownames(metadata)
		cat("generating counts table\n")
		write.csv(countsTable, file = paste0(destdir,"/", gse_acc,"/countsTable.csv"))
		
		transcripts <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE, dropInfReps=TRUE, importer = read.delim, countsFromAbundance = "lengthScaledTPM")
		save(transcripts, file=paste(destdir,"/", gse_acc, "/salmon/","transcripts.RData", sep=""))
		transcripts_counts <- transcripts$counts
		colnames(transcripts_counts) <- rownames(metadata)
		write.csv(transcripts_counts, file = paste0(destdir,"/", gse_acc,"/transcripts_counts.csv"))
	} else {
		stop('salmon file missing')
	}
}