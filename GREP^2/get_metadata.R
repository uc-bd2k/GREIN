
get_metadata <- function(gse_acc, sra_study_acc, destdir ) {

	# Download experiment and sample metadata.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	sra_study_acc: SRA study accession number.
	#   destdir: Directory to save the results.
	#
	# Returns:
	#   Combined metadata.
	
	setwd(paste0(destdir,"/",gse_acc))
	system(paste("rm -rf *matrix.txt.gz *_metadata.csv"))
	
	repeat {
		geo_list_prim <- try(lapply(getGEO(gse_acc,GSEMatrix=TRUE,destdir = getwd() ,getGPL=FALSE), function(x) pData(phenoData(x))), silent=TRUE)
		if(!is(geo_list_prim, 'try-error')) break
		Sys.sleep(3)
	}
	geo_list_sec <- lapply(geo_list_prim, function(x) x[,Reduce(intersect, lapply(geo_list_prim, function(x) colnames(x)))])
	geo_df <- do.call(rbind, geo_list_sec)
	closeAllConnections()
	metadata_geo <- data.frame(lapply(geo_df, as.character), stringsAsFactors=FALSE)
	
	# sample specific metadata
	system(paste("wget -O ",gse_acc,
	"_metadata.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=",sra_study_acc,"'", sep=""))
	metadata_sra <- data.frame(lapply( read.csv(file=paste(gse_acc,"_metadata.csv", sep=""), header=TRUE)
	, as.character), stringsAsFactors=FALSE)
		
	if(nrow(metadata_geo[which(metadata_geo$library_strategy=="RNA-Seq"),])>1){
		if(length(which(metadata_geo$geo_accession %in% metadata_sra$SampleName))>1){
			metadata <- merge(metadata_geo, metadata_sra, by.x = "geo_accession", by.y="SampleName")
			rownames(metadata) <- paste0(metadata$Run,"_",metadata$geo_accession)
			metadata <- metadata[which(metadata$LibraryStrategy=="RNA-Seq"),]
			metadata <- metadata[!is.na(metadata$Run),]
			metadata <- metadata[ , !apply( metadata , 2 , function(x) all(is.na(x)) ) ]
			
			cat("save metadata\n")
			setwd(paste0(destdir,"/",gse_acc))
			write.csv(metadata, file = "metadata.csv")
			system(paste("rm -rf *matrix.txt.gz *_metadata.csv"))
		} else {
			x=as.character(unlist(metadata_geo[,names(Filter(function(u) any(grepl('BioSample: https://www.ncbi.nlm.nih.gov/',u)), metadata_geo))]))
			metadata_geo$BioSample_2 <- sub("BioSample: https://www.ncbi.nlm.nih.gov/biosample/","",x[grep("SAMN",x)])
			metadata_sra$BioSample_sra <- metadata_sra$BioSample
			metadata_geo$BioSample_geo <- metadata_geo$BioSample_2
			metadata <- merge(metadata_geo, metadata_sra, by.x = "BioSample_2", by.y="BioSample")[,-1]
		}
	} else {
		setwd(destdir)
		system(paste0("rm -rf ",paste0(destdir,"/",gse_acc)))
	}
}
