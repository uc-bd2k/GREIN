
sub_metadata <- function(metadata) {
	toDelete <- vector()
	for (i in seq_along(metadata)) {
		if (length(unique(metadata[,i])) == 1) {
			toDelete <- c(toDelete, colnames(metadata)[i])
		}
	}
	pdata <- metadata[, !colnames(metadata) %in% toDelete]
	pdata <- pdata[, !names(pdata) %in% c(grep("relation", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("file", names(pdata), value = TRUE, ignore.case=TRUE),	grep("description", names(pdata), value = TRUE, ignore.case=TRUE), 
	grep("spot", names(pdata), value = TRUE, ignore.case=TRUE), grep("process", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("date", names(pdata), value = TRUE, ignore.case=TRUE), grep("protocol", names(pdata), value = TRUE, ignore.case=TRUE), 
	grep("accession", names(pdata), value = TRUE, ignore.case=TRUE),
	#grep("title", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("instrument", names(pdata), value = TRUE, ignore.case=TRUE), grep("library", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("size", names(pdata), value = TRUE, ignore.case=TRUE),grep("platform", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("paired", names(pdata), value = TRUE, ignore.case=TRUE),grep("length", names(pdata), value = TRUE, ignore.case=TRUE),
	"Run","avgLength","spots","Model","platform_id","bases","size_MB","download_path","Experiment", "Sample","BioSample",	
	"SampleName","RunHash","ReadHash"), drop=FALSE]

	if(is.null(dim(pdata))) {
		return(message("Either the variables do not contain meaningful information or they contain single value. Please see the full metadata table."))
	} else {
		prop <- which(apply(pdata, 2, function(x) any(grepl(": ", x))))
		x <- sapply(names(prop), function(x) length(unique(gsub(":.*$", "", pdata[,x]))))
		if(any(x>1)){
			pdata <- pdata[, !(names(pdata) %in% names(prop[which(x>1)])), drop=F]
		}
		prop2 <- which(apply(pdata, 2, function(x) any(grepl(": ", x))))
		for (i in seq_along(prop2)) {
			colnames(pdata)[prop2[i]] <- unique(gsub(":.*$", "", pdata[,prop2[i]]))
			pdata[,prop2[i]] <- gsub(".*: ", "", pdata[,prop2[i]])
			pdata[,prop2[i]] <- 
				if (class(pdata[,prop2[i]])=="integer") {
					as.integer(pdata[,prop2[i]])
				} else if (class(pdata[,prop2[i]])=="numeric") {
					as.integer(pdata[,prop2[i]])						
				} else {
					as.factor(pdata[,prop2[i]])
				}
		}
	}
	if(dim(pdata)[2]==1) {
		pdata <- pdata ;
		colnames(pdata) <- "characteristics"
	} else if (dim(pdata)[2]==0) {
		pdata <- metadata[, colnames(metadata) %in% grep("title", colnames(metadata), value = TRUE), drop=FALSE] ;
		colnames(pdata) <- "characteristics"
	} else {
		#pdata <- pdata[, !names(pdata) %in% "source_name_ch1", drop=FALSE]
		pdata <- pdata
	}
	
	pdata_final <- pdata[, !names(pdata) %in% c(grep("relation", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("file", names(pdata), value = TRUE, ignore.case=TRUE),	grep("description", names(pdata), value = TRUE, ignore.case=TRUE), 
	grep("spot", names(pdata), value = TRUE, ignore.case=TRUE), grep("process", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("date", names(pdata), value = TRUE, ignore.case=TRUE), grep("protocol", names(pdata), value = TRUE, ignore.case=TRUE), 
	grep("accession", names(pdata), value = TRUE, ignore.case=TRUE),
	#grep("title", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("instrument", names(pdata), value = TRUE, ignore.case=TRUE), grep("library", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("size", names(pdata), value = TRUE, ignore.case=TRUE),grep("platform", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("paired", names(pdata), value = TRUE, ignore.case=TRUE),grep("length", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("provider", names(pdata), value = TRUE, ignore.case=TRUE), grep("^id$", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("status", names(pdata), value = TRUE, ignore.case=TRUE),
	grep("name", names(pdata), fixed = TRUE),grep("sex", names(pdata), value = TRUE, ignore.case=TRUE),
	"Run","avgLength","spots","Model","platform_id","bases","size_MB","download_path","Experiment", "Sample","BioSample",	
	"SampleName","RunHash","ReadHash"), drop=FALSE]
	
	if (dim(pdata_final)[2]>=1 & "source_name_ch1" %in% names(pdata_final)) {
		pdata_final <- pdata_final[, -which(colnames(pdata_final)=="title"), drop=F]
		colnames(pdata_final)[which(names(pdata_final) %in% "source_name_ch1")] <- "characteristics"
		if(ncol(pdata_final[!duplicated(as.list(toupper(pdata_final)))])!= ncol(pdata_final)){
			pdata_final <- pdata_final[,-which(colnames(pdata_final)=="characteristics"), drop=F]
		} else {
			pdata_final <- pdata_final
		}
	} else {
		pdata_final <- pdata_final
	}
	
	if ("title" %in% names(pdata_final) & "characteristics" %in% names(pdata_final)){
		pdata_final <- pdata_final[, -which(colnames(pdata_final)=="title"), drop=F]
	} else if("title" %in% names(pdata_final) & !("characteristics" %in% names(pdata_final))) {
		colnames(pdata_final)[which(names(pdata_final) %in% "title")] <- "characteristics"
		pdata_final <- pdata_final
	} else {
		pdata_final <- pdata_final
	}
	
	if ("title" %in% names(pdata_final) & "characteristics" %in% names(pdata_final)){
		pdata_final <- pdata_final[, -which(colnames(pdata_final)=="title"), drop=F]
	} else {
		pdata_final <- pdata_final
	}
	
	return(pdata_final)
}
