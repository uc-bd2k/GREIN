
sub_metadata <- function(metadata) {
	toDelete <- vector()
	for (i in seq_along(metadata)) {
		if (length(unique(metadata[,i])) == 1) {
			toDelete <- c(toDelete, colnames(metadata)[i])
		}
	}
	pdata <- metadata[, !colnames(metadata) %in% toDelete]
	pdata <- pdata[, !names(pdata) %in% c(grep("relation", names(pdata), value = TRUE),grep("file", names(pdata), value = TRUE),
	grep("description", names(pdata), value = TRUE), grep("spot", names(pdata), value = TRUE), grep("data_processing", names(pdata), value = TRUE),
	grep("Date", names(pdata), value = TRUE), grep("protocol", names(pdata), value = TRUE), grep("geo_accession", names(pdata), value = TRUE),
	grep("date", names(pdata), value = TRUE), grep("title", names(pdata), value = TRUE),grep("instrument", names(pdata), value = TRUE),
	"Run","avgLength","spots","Model","platform_id","bases","size_MB","download_path","Experiment", "Sample","BioSample",	
	"SampleName","RunHash","ReadHash"), drop=FALSE]

	if(is.null(dim(pdata))) {
		return(message("Either the variables do not contain meaningful information or they contain single value. Please see the full metadata table."))
	} else {
		prop <- which(apply(pdata, 2, function(x) any(grepl(": ", x))))

		for (i in seq_along(prop)) {
			if(length(unique(gsub(":.*$", "", pdata[,prop[i]])))>1) {
				pdata <- pdata[,-prop[i], drop=FALSE]
			} else {
				colnames(pdata)[prop[i]] <- unique(gsub(":.*$", "", pdata[,prop[i]]))
				pdata[,prop[i]] <- gsub(".*: ", "", pdata[,prop[i]])
				pdata[,prop[i]] <- if (class(pdata[,prop[i]])=="numeric") {as.integer(pdata[,prop[i]])} else {as.factor(pdata[,prop[i]])}
			}
		}
	}
	if(dim(pdata)[2]==1) {
		pdata <- pdata
	} else if (dim(pdata)[2]==0) {
		pdata <- metadata[, colnames(metadata) %in% grep("title", colnames(metadata), value = TRUE), drop=FALSE]
		colnames(pdata) <- "property"
	} else {
		pdata <- pdata[, !names(pdata) %in% "source_name_ch1", drop=FALSE]
	}
	
	return(pdata)
}
