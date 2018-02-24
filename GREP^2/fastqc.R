
fastqc <- function(gse_acc,sra_run_acc=metadata$Run, destdir, n_thread ) {
	
	# Generate run wise FastQC report and a combined FastQC report using MultiQC.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	sra_run_acc: SRA run accession number.
	#   destdir: Directory to save the results.
	#	n_thread: Number of processors to use.
	#
	# Returns:
	#   A combined html report of all the fastqc runs' report along with the single reports. 

	# Location of FastQC software for pre alignment QC
	
	setwd(paste0(destdir,"/",gse_acc))
	fastqc.path = "/home/FastQC/fastqc"
	fastq.files = list.files(paste(destdir,"/",gse_acc,"/",sra_run_acc,sep=""), pattern=".fastq$", full=TRUE)
	fastqc_f <- list.files(paste(destdir,"/",gse_acc,"/fastqc",sep=""), pattern=".zip$", full=F)
	
	if(!dir.exists("fastqc")){
		system(paste("mkdir ",paste(destdir,"/",gse_acc,sep=""),"/fastqc",sep=""))
	}
	setwd(paste(destdir,"/",gse_acc,"/fastqc/",sep=""))
	#if(file.exists("multiqc_report.html")) {
	if(length(fastq.files)==length(fastqc_f)){
		cat("file exists")
	} else {	
		mclapply(fastq.files,function(x) {
			if(!file.exists(paste0(destdir,"/",gse_acc,"/fastqc/",unlist(strsplit(basename(x), ".fastq")),"_fastqc.zip"))){
				#run FastQC on each of the files
				system(paste(fastqc.path," ",x,sep=""))
				#move the files created by FastQC to the output directroy
				system(paste("mv ",sub(".fastq","_fastqc",x),"* ",out.dir=paste(destdir,"/",gse_acc,"/fastqc",sep=""),sep=""))
			} else {
				cat('file exists')
			}
		}, mc.cores=n_thread)  	
	}
	setwd(paste(destdir,"/",gse_acc,sep=""))
	system(paste0("touch fastqc_completed"))	
}