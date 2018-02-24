
salmon <- function(gse_acc,species,sra_run_acc=metadata$Run,library_layout=metadata$LibraryLayout,
							avg_length=metadata$avgLength, index_dir,destdir, n_thread ) 
{

	# Quantify transcripts based on a target transcript file (index files) 
	#
	# Args:
	#   gse_acc: GSE accession number.
	#	species: Either Homo sapiens, Mus musculus, or Rattus norvegicus.
	#	sra_run_acc: SRA run accession number.
	#	library_layout: If the run is single end or paired.
	#	avg_length: Average length of the fragment.
	#   index_dir: Directory of the species indexing files.
	#   destdir: Directory to save the results.
	#	n_thread: Number of processors to use.
	#
	# Returns:
	#   1. quant.sf: plain-text, tab-separated quantification file that contains 5 column: Name,Length,EffectiveLength,TPM, and NumReads.
	#	2. cmd_info.json: A JSON format file that records the main command line parameters with which Salmon was invoked for the run 
	#		that produced the output in this directory.
	#	3. aux_info: This directory will have a number of files (and subfolders) depending on how salmon was invoked.
	#	4. meta_info.json: A JSON file that contains meta information about the run, including stats such as the number of observed 
	#		and mapped fragments, details of the bias modeling etc. 
	#	5. ambig_info.tsv: This file contains information about the number of uniquely-mapping reads as well as the total number of 
	#		ambiguously-mapping reads for each transcript. 
	#	6. lib_format_counts.json: This JSON file reports the number of fragments that had at least one mapping compatible with the 
	#		designated library format, as well as the number that didn’t. 
	#	7. libParams: The auxiliary directory will contain a text file called flenDist.txt. This file contains an approximation of the 
	#		observed fragment length distribution.


	setwd(paste0(destdir,"/",gse_acc))
	if(!dir.exists("salmon")){
		system(paste("mkdir ",paste(destdir,"/",gse_acc,sep=""),"/salmon",sep=""))
	}

	metadata <- read.csv(file="metadata.csv", header=TRUE, row.names=1)
	
	transcript_idx_redo= function(species) {
		if (species == "Homo sapiens") {
			return(paste0(index_dir,"/human_transcripts_rel87_index_kmer17"))
		}
		else if (species == "Rattus norvegicus") {
			return(paste0(index_dir,"/rat_transcripts_rel87_index_kmer17"))
		}
		else if (species == "Mus musculus") {
			return(paste0(index_dir,"/mouse_transcripts_rel87_index_kmer17"))
		}
		else {
			return(NULL)
		}
	}

	for(j in 1: length(sra_run_acc)) {		
		setwd(paste(destdir,"/",gse_acc,"/salmon",sep=""))
		print(j)
		
		transcript_idx= function(species) {
			if (species == "Homo sapiens") {
				if(avg_length[j]<=51){
					return(paste0(index_dir,"/human_transcripts_rel87_index_kmer17"))
				} else {
					return(paste0(index_dir,"/human_transcripts_rel87_index"))
				}
			}
			else if (species == "Rattus norvegicus") {
				if(avg_length[j]<=51){
					return(paste0(index_dir,"/rat_transcripts_rel87_index_kmer17"))
				} else {
					return(paste0(index_dir,"/rat_transcripts_rel87_index"))
				}
			}
			else if (species == "Mus musculus") {
				if(avg_length[j]<=51){
					return(paste0(index_dir,"/mouse_transcripts_rel87_index_kmer17"))
				} else {
					return(paste0(index_dir,"/mouse_transcripts_rel87_index"))
				}
			}
			else {
				return(NULL)
			}
		}

		if(!file.exists(paste0(rownames(metadata)[j],"_transcripts_quant/", rownames(metadata)[j],"_quant_new.sf"))){
			if (library_layout[j]=="SINGLE") {
				system(paste("salmon quant -i ",transcript_idx(species), " -p ", n_thread, " --fldMean ", 
				avg_length[j], " --fldSD 25 -l A -r ", destdir,"/",gse_acc,"/",sra_run_acc[j],"/", 
				sra_run_acc[j], "_pass.fastq -o ", rownames(metadata)[j],"_transcripts_quant", sep=""))				
			} else {
				system(paste("salmon quant -i ",transcript_idx(species), " -p ", n_thread, " -l A -1 ",
				destdir,"/",gse_acc,"/",sra_run_acc[j],"/", sra_run_acc[j], "_pass_1.fastq ", "-2 ",
				destdir,"/",gse_acc,"/",sra_run_acc[j],"/", sra_run_acc[j], "_pass_2.fastq -o ", 
				rownames(metadata)[j],"_transcripts_quant", sep=""))				
			}
			setwd(paste(destdir,"/",gse_acc,"/salmon/",rownames(metadata)[j],"_transcripts_quant",sep=""))
			file.rename(paste0("quant.sf"), paste0(rownames(metadata)[j],"_quant.sf"))
			if (file.exists(paste0(rownames(metadata)[j],"_quant.sf"))) {
				system(paste("cat ",rownames(metadata)[j],"_quant.sf| sed -E 's/\\.[0-9]+//' > ",rownames(metadata)[j],"_quant_new.sf", sep=""))
			} else {
				cat(paste0(rownames(metadata)[j],"_quant.sf doesn't exist. Processing next sample."))
			}
		} else {
			cat('file exists')
		}
	}
}