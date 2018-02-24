
downloadGEO_rnaseq <- function(gse_acc, destdir, prefetch_workspace_sra, source_code_dir, n_thread, dir_size) {
	
	# Download data Pipeline.
	#
	# Args:
	#   gse_acc: GSE accession number
	#	destdir: Directory where the data are saved
	#	prefetch_workspace_sra: Initial location of the downloaded sra files created by the sra-toolkit itself
	#	source_code_dir: Code directory
	#	n_thread: Number of threads
	#   dir_size: Size of the data directory
	#
	# Returns:
	#   Downloaded sra files and metadata

		if(as.numeric(unlist(strsplit(gsub("T.*","",system(paste0("du -s ", destdir), intern=T)), "\t"))[1])/1000000000 < dir_size){
		if (!dir.exists(paste(destdir,"/",gse_acc,sep=""))){
			if(!file.exists(paste(destdir,"/",gse_acc,"/", "download_completed",sep=""))) {
		
				library(GEOquery)
				library(Biobase)
				library(parallel)
				
				load("data/geo_QuerynBrowser_merged_jan28_2018.RData")
				sra_study_acc <- geo_QuerynBrowser_merged_jan28_2018[which(geo_QuerynBrowser_merged_jan28_2018$gse_id==gse_acc), "study"]
				species <- geo_QuerynBrowser_merged_jan28_2018[which(geo_QuerynBrowser_merged_jan28_2018$gse_id==gse_acc), "taxon"]
				
				if(!dir.exists(paste0(destdir,"/",gse_acc)) & !dir.exists(paste0(destdir,"/",gse_acc))){
					system(paste0("mkdir ",paste0(destdir,"/",gse_acc)))
				}
				system(paste0("chmod -R 777 ",paste0(destdir,"/",gse_acc)))
				
				setwd(paste0(destdir,"/",gse_acc))
				###### Download metadata ######
				source(paste0(source_code_dir,"/get_metadata.R"))
				if (!file.exists("metadata.csv")) {
					metadata <- get_metadata(gse_acc, sra_study_acc,destdir)
				} else {
					metadata <- read.csv(file="metadata.csv", header=TRUE, row.names=1)
				}
				if(length(list.files(path = ".", full.names = TRUE, recursive = FALSE))==0){
					system(paste0("rm -rf ",paste0(destdir,"/",gse_acc)))
				}

				if(nrow(metadata)>1) {
					###### Download sra files for each run ######
					source(paste0(source_code_dir,"/get_sra.R"))
					get_sra(gse_acc,sra_run_acc=metadata$Run,destdir,prefetch_workspace_sra,n_thread, dir_size)
					
					empty <- setdiff(metadata$Run, sub("\\/.*", "", list.files(pattern = "*.sra$", recursive = TRUE)))
					if (length(empty)==0) {
						system(paste0("touch download_completed"))				
					} else {
						get_sra(gse_acc,sra_run_acc=empty,destdir,prefetch_workspace_sra,n_thread, dir_size)
						system(paste0("touch download_completed"))
					}
					
				} else {
					cat(paste0(gse_acc," has only 1 sample","\n"))
					setwd(destdir)
					system(paste0("rm -rf ",paste0(destdir,"/",gse_acc)))
				}					
			} else {
				cat(paste0("download has already been completed for ",gse_acc,"\n"))
			}
		} else {
			cat(paste0("download has already been completed for ",gse_acc,"\n"))
		}
	} else {
		stop('Disk is full')	
	}
}