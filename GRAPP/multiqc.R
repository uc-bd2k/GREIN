
multiqc <- function(gse_acc, destdir) {
	
	# Generate combined QC report for Salmon and FastQC.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#   destdir: Directory to save the results.
	#
	# Returns:
	#   A combined html report. 

	if(!file.exists(paste0(destdir,"/",gse_acc,"/fastqc/multiqc_report_salmon.html"))) {
		if(file.exists(paste0(destdir,"/",gse_acc,"/fastqc/multiqc_report.html"))){
			setwd(paste0(destdir,"/",gse_acc,"/fastqc"))
			system("rm multiqc_report.html")
		}
		setwd(paste0(destdir,"/",gse_acc))
		cat(paste("Creating MultiQC report.\n",sep=""))
		system(paste("multiqc ./fastqc ./salmon ","-o ",paste(destdir,"/",gse_acc,"/fastqc",sep=""), sep=""))
		
		setwd(paste(destdir,"/",gse_acc,"/fastqc",sep=""))
		system(paste0("cp \"multiqc_report.html\" \"multiqc_report_salmon.html\""))
	} else {
		cat("Multiqc report already exists")
	}
}