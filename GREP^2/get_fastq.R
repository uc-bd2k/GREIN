
get_fastq <- function(gse_acc,sra_run_acc=metadata$Run, library_layout=metadata$LibraryLayout, destdir, n_thread, dir_size) {
	
	# Generate fastq files from the sra files.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	sra_run_acc: SRA run accession number.
	#	library_layout: If the run is single end or paired.
	#   destdir: Directory to save the results.
	#	n_thread: Number of processors to use.
	#
	# Returns:
	#   .fastq files for each run.
	
	
	mclapply(1: length(sra_run_acc),function(i) {
		if(as.numeric(unlist(strsplit(gsub("T.*","",system(paste0("du -s ", destdir), intern=T)), "\t"))[1])/1000000000 < dir_size){
			if (library_layout[i]=="SINGLE") {
				setwd(paste(destdir,"/",gse_acc,"/",sra_run_acc[i],sep=""))
				if(length(list.files(pattern=".fastq"))==1) {
					cat("processing next sample")
				} else {
					system (paste("fastq-dump --outdir ",paste(destdir,"/",gse_acc,"/", sra_run_acc[i], sep=""),
					" --skip-technical  --readids --read-filter pass --dumpbase --split-spot --clip ",
					paste("./",sra_run_acc[i],".sra", sep=""), sep=""))
				}
			} else {
				setwd(paste(destdir,"/",gse_acc,"/",sra_run_acc[i],sep=""))
				if(length(list.files(pattern=".fastq"))==2) {
					cat("processing next sample")
				} else {
					do.call(file.remove, list(list.files(pattern = "*.fastq")))
					system (paste("fastq-dump --outdir ",paste(destdir,"/",gse_acc,"/", sra_run_acc[i], sep=""),
					" --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip ",
					paste("./",sra_run_acc[i],".sra", sep=""), sep=""))
				}
			}
		} else {
			stop('Disk is full')
		}
	}, mc.cores=n_thread)
	cat(paste("All fastq files are generated successfully. ",Sys.time(),"\n",sep=""))
}