
trim_fastq <- function(gse_acc,sra_run_acc=metadata$Run,instrument_model= metadata$instrument_model,library_layout=metadata$LibraryLayout, destdir,n_thread){

	# Trim the fastq files to improve it's quality.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	sra_run_acc: SRA run accession number.
	#	instrument_model: Instrument used to generate raw reads. (e.g., Illumina HiSeq, MiSeq etc.)
	#	library_layout: If the run is single end or paired.
	#   destdir: Directory to save the results.
	#	n_thread: Number of processors to use.
	#
	# Returns:
	#   Trimmed fastq files for each run.
	
	adapters<- function(instrument_model) {
		if (grepl("HiSeq|MiSeq",instrument_model)) {
			return("TruSeq3")
		}
		else if (grepl("GA|Genome Analyzer",instrument_model)) {
			return("TruSeq2")
		}
		else if (grepl("NextSeq",instrument_model)) {
			return("NexteraPE")
		}
		else {
			stop()
		}
	}

	path.adaptors="/home/Trimmomatic-0.36/adapters/"
	mclapply(1: length(sra_run_acc),function(i) {
		if (library_layout[i]=="SINGLE") {
			setwd(paste(destdir,"/",gse_acc,"/",sra_run_acc[i],sep=""))
			fq = paste(sra_run_acc[i],"*.fastq",sep="")
			fq_se =  paste(sra_run_acc[i],"_se.fastq",sep="")
			system(paste("java -jar /home/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ",fq, " ", fq_se, " ",	"ILLUMINACLIP:",
			path.adaptors,adapters(instrument_model[i]),"-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",sep=""))
		} else {
			setwd(paste(destdir,"/",gse_acc,"/",sra_run_acc[i],sep=""))
			fq1 = paste(sra_run_acc[i],"*_1.fastq",sep="")
			fq2 = paste(sra_run_acc[i],"*_2.fastq",sep="")
			fq1_paired = paste(sra_run_acc[i],"_1_paired.fastq",sep="")
			fq2_paired = paste(sra_run_acc[i],"_2_paired.fastq",sep="")
			fq1_unpaired = paste(sra_run_acc[i],"_1_unpaired.fastq",sep="")
			fq2_unpaired = paste(sra_run_acc[i],"_2_unpaired.fastq",sep="")
			
			system(paste("java -jar /home/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ",
			fq1," ",fq2," ",fq1_paired," ",fq1_unpaired," ",fq2_paired," ",fq2_unpaired," ","ILLUMINACLIP:",
			path.adaptors,adapters(instrument_model[i]),"-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",sep=""))
		}
		
	}, mc.cores=n_thread)
}