
get_sra <- function(gse_acc,sra_run_acc=metadata$Run,destdir,prefetch_workspace_sra, n_thread, dir_size) {
	
	# Download sra files.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	sra_run_acc: SRA run accession number.
	#   destdir: Directory to save the results.
	#   prefetch_workspace_sra: Default directory of the sra-toolkit to save the sra files. This can be changed by using ./vdb-config -i.
	#	n_thread: Number of processors to use.
	#	dir_size: Size of the directory
	#
	# Returns:
	#   .sra files for each run.

	mclapply(1: length(sra_run_acc),function(i) {
		if(as.numeric(unlist(strsplit(gsub("T.*","",system(paste0("du -s ", destdir), intern=T)), "\t"))[1])/1000000000 < dir_size){
			system(paste0("mkdir ",paste0(destdir,"/",gse_acc,"/",sra_run_acc[i])))
			setwd(paste(destdir,"/",gse_acc,"/",sra_run_acc[i],sep=""))
			if(!file.exists(paste0(destdir,"/",gse_acc,"/",sra_run_acc[i],"/",sra_run_acc[i],".sra"))) {
			#the following command is going to save the run files in 'sra' folder within prefetch workspace and ref seq within 'refseq' folder.
				system(paste0("rm ",prefetch_workspace_sra, "/", "*.lock"))
				system(paste0("prefetch -X 500G --ascp-path '/home/.aspera/connect/bin/ascp|/home/.aspera/connect/etc/asperaweb_id_dsa.openssh' --ascp-options '-k 1 -QT -l 400m' ",	sra_run_acc[i]))
				system(paste0("rm ",prefetch_workspace_sra, "/", "*.lock"))
				system(paste0("mv ",prefetch_workspace_sra, "/",sra_run_acc[i],".sra* ", getwd()))
			} else {
				cat(paste0("download has already been completed for ",sra_run_acc[i],"\n"))
			}
		} else {
			stop('Disk is full')
		}

	}, mc.cores=n_thread)
		
	# If any of the runs is not downloaded for any technical reason, then do the following:
	setwd(paste0(destdir,"/",gse_acc))
	empty <- setdiff(sra_run_acc, sub("\\/.*", "", list.files(pattern = "*.sra$", recursive = TRUE)))
	
	if (length(empty)!=0) {
		repeat {
			mclapply(1: length(empty),function(i) {
				system(paste0("rm ",prefetch_workspace_sra, "/", "*.lock"))
				setwd(paste(destdir,"/",gse_acc,"/",empty[i],sep=""))
				system(paste0("prefetch -X 500G --ascp-path '/home/.aspera/connect/bin/ascp|/home/.aspera/connect/etc/asperaweb_id_dsa.openssh' --ascp-options '-k 1 -QT -l 400m' ", empty[i]))
				system(paste0("rm ",prefetch_workspace_sra, "/", "*.lock"))
				system(paste0("mv ",prefetch_workspace_sra,"/",empty[i],".sra* ", getwd()))
			}, mc.cores=n_thread)
			setwd(paste0(destdir,"/",gse_acc))
			empty <- setdiff(sra_run_acc, sub("\\/.*", "", list.files(pattern = "*.sra$", recursive = TRUE)))
			if (length(empty)==0) break
			cat(paste("All SRA files are downloaded successfully. ",Sys.time(),"\n",sep=""))
		}
	} else {
		cat(paste("All SRA files are downloaded successfully. ",Sys.time(),"\n",sep=""))
	}
}
