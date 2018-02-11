
processGEO_rnaseq <- function(gse_acc, destdir, index_dir, source_code_dir, n_thread, dir_size) {
	
	# Generate fastq files, QC report, trim the fastq reads, and finally produce counts table.
	#
	# Args:
	#   gse_acc: GSE accession number.
	#   destdir: Directory to save the results.
	#   index_dir: Directory of the species indexing files.
	#	source_code_dir: Directory of the source codes used in this function.
	#	n_thread: Number of processors to use.
	#
	# Returns:
	#   Fastqc report, alignment report, counts table, eset, exploratory plots, and log file.

	if (!file.exists(paste(destdir,"/",gse_acc,"/", "process_completed",sep="")) && 
		!dir.exists(paste("/opt/raid10/genomics/naim/Thesis/geo/datasets/",gse_acc,"/salmon",sep=""))) {
	
			if (file.exists(paste(destdir,"/",gse_acc,"/", "download_completed",sep=""))) {
			
				library(GEOquery)
				library(Biobase)
				library(parallel)
				library(tximport)
				library(readr)
				library(biomaRt)
				library(EnsDb.Hsapiens.v86)
				library(EnsDb.Rnorvegicus.v79)
				library(EnsDb.Mmusculus.v79)
				library(AnnotationDbi)
				library(org.Hs.eg.db)
				library(org.Mm.eg.db)
				library(org.Rn.eg.db)
				library(reshape2)

				if (file.exists(paste(destdir,"/",gse_acc,"/", "log.Rout",sep=""))) {
					system("rm -f log.Rout")
				}
				
				if (file.exists(paste(destdir,"/",gse_acc,"/", gse_acc,"_series_matrix.txt.gz",sep=""))) {
					system(paste0("rm *_series_matrix.txt.gz"))
					closeAllConnections()
				}
				
				load("/data/geo_QuerynBrowser_merged.RData")
				sra_study_acc <- geo_QuerynBrowser_merged[which(geo_QuerynBrowser_merged$gse_id==gse_acc), "study"]
				species <- geo_QuerynBrowser_merged[which(geo_QuerynBrowser_merged$gse_id==gse_acc), "taxon"]
				
				### Download metadata
				source(paste0(source_code_dir,"/get_metadata.R"))
				setwd(paste0(destdir,"/",gse_acc))
				if (!file.exists("metadata.csv")) {
					metadata <- get_metadata(gse_acc, sra_study_acc,destdir)
				} else {
					metadata <- read.csv(file="metadata.csv", header=TRUE, row.names=1)
				}

				setwd(paste0(destdir,"/",gse_acc))
				con = file("log.Rout", open="wt")
				sink(con, type="output",  append=TRUE,split = FALSE)
				sink(con, type="message", append=TRUE,split = FALSE)	
					
				### get fastq files
				cat(paste("STEP 1: Generate fastq files from the sra files... ",Sys.time(),"\n",sep=""))
				source(paste0(source_code_dir,"/get_fastq.R"))
				get_fastq(gse_acc,sra_run_acc=metadata$Run, library_layout=metadata$LibraryLayout, destdir, n_thread, dir_size)
					
				setwd(paste0(destdir,"/",gse_acc))
				n_fastq <- (length(which(metadata$LibraryLayout=="PAIRED"))*2)+length(which(metadata$LibraryLayout=="SINGLE"))
				fastq_dumped <- length(list.files(pattern = "\\.fastq$",recursive=T,full.names=F))
				
				if(n_fastq == fastq_dumped){	
					### fastqc and Quality trimming
					cat("STEP 2: Quality check and trimming reads\n")
						
					## Trimmomatic
					cat(paste("STEP 3: Filtering reads Based on Quality of reads... ",Sys.time(),"\n",sep=""))
					source(paste0(source_code_dir,"/trim_fastq.R"))
					trim_fastq(gse_acc,sra_run_acc=metadata$Run,instrument_model= metadata$instrument_model,library_layout=metadata$LibraryLayout, 
					destdir,n_thread)

					## FastQC and MultiQC
					cat(paste("STEP 4: Quality check of the fastq files... ",Sys.time(),"\n",sep=""))
					source(paste0(source_code_dir,"/fastqc.R"))
					fastqc(gse_acc,sra_run_acc=metadata$Run, destdir, n_thread)
						
					### Salmon
					cat(paste("STEP 5: Mapping reads with Salmon... ",Sys.time(),"\n",sep=""))
					source(paste0(source_code_dir,"/salmon.R"))
					salmon(gse_acc,species,sra_run_acc=metadata$Run,library_layout=metadata$LibraryLayout,
					avg_length=metadata$avgLength, index_dir,destdir, n_thread)

					all_files <- file.exists(paste(destdir, "/",gse_acc, "/salmon/", rownames(metadata),"_transcripts_quant/",rownames(metadata), "_quant_new.sf", sep=""))
					files <- file.path(paste(destdir, "/",gse_acc, "/salmon/", rownames(metadata)[all_files],"_transcripts_quant/",rownames(metadata)[all_files], "_quant_new.sf", sep=""))

					if(length(files) == dim(metadata)[1]){

						cat(paste("STEP 6: Gene level counts from transcripts... ",Sys.time(),"\n",sep=""))
						source(paste0(source_code_dir,"/tximport2counts.R"))
						tximport2counts(gse_acc,species,destdir)
						
						### Multiqc report
						cat(paste("STEP 6: MultiQC report... ",Sys.time(),"\n",sep=""))
						source(paste0(source_code_dir,"/multiqc.R"))
						multiqc(gse_acc, destdir)
						
						### Create eset
						cat(paste("STEP 7: Create eset ... ",Sys.time(),"\n",sep=""))
						source(paste0(source_code_dir,"/create_eset.R"))
						create_eset(gse_acc,species, destdir)
							
						### Remove fastq files
						cat(paste("Remove fastq and sra files... ",Sys.time(),"\n",sep=""))
						setwd(paste(destdir,"/",gse_acc,sep=""))
						system(paste("rm -rf SRR*", sep=""))
							
						cat(paste("Process completed. ",Sys.time(),"\n","\n",sep=""))
						print(sessionInfo())
						sink(type="output")
						sink(type="message")
						
						err=grep("Error", readLines("log.Rout"), value = TRUE)
						if (length(err)==0){
							setwd(paste(destdir,"/",gse_acc,sep=""))
							system(paste0("touch process_completed"));
						} else {
							return(warning("something is wrong. Check the log file"))
						}
						closeAllConnections()
					} else {
						setwd(destdir)
						system(paste0("touch ",gse_acc,"_missing_salmon_files"))
						warning('missing salmon files')
					}
				} else {
					system(paste0("touch ",gse_acc,"_fastq_download_incomplete"))
					warning('incomplete fastq download')
				}
			} else {
				return(warning(paste0(gse_acc," :incomplete download. Processing the next study...")))
			}
	} else {
		return(warning(paste0(gse_acc," has been processed successfully. Moving to the next one...")))
	}
}

	
