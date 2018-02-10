
datatable_to_use <- function() {
	load("/opt/raid10/genomics/naim/Thesis/geo/geo_QuerynBrowser_merged_jan28_2018.RData")
	library(Biobase)
	geo=list.files("/opt/raid10/genomics/naim/Thesis/geo/datasets/")
	datatable <- matrix(, nrow = length(geo), ncol = 6)
	for ( i in seq_along(geo)) {
		if (file.exists(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",geo[i],"/metadata.csv"))) {
			metadata <- read.csv(file=paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",geo[i],"/metadata.csv"), header=TRUE, row.names=1)
		} else {
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",geo[i],"/eset.RData"))
			metadata <- pData(eset)
		}
			
		datatable[i,] <- cbind(
		geo[i],
		#paste0(
		#'<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
		#geo[i], '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'geoinfo\', 1)">', geo[i], '</a>'),
		#length(unique(unlist(lapply(colnames(exprs(eset)), function(x) unlist(strsplit(x,"_"))[2])))),
		length(unique(unlist(lapply(rownames(metadata), function(x) unlist(strsplit(x,"_"))[2])))),
		geo_QuerynBrowser_merged_jan28_2018[geo_QuerynBrowser_merged_jan28_2018$gse_id %in% geo[i], 4],
		geo_QuerynBrowser_merged_jan28_2018[geo_QuerynBrowser_merged_jan28_2018$gse_id %in% geo[i], 5], 
		geo_QuerynBrowser_merged_jan28_2018[geo_QuerynBrowser_merged_jan28_2018$gse_id %in% geo[i], 6],
		#nrow(pData(eset)))
		as.integer(nrow(metadata)))
	}
	datatable_now <- as.data.frame(datatable[complete.cases(datatable),], stringsAsFactors=FALSE)
	datatable_to_use2 <- data.frame(datatable_now[,1], as.integer(datatable_now[,2]), datatable_now[,3],
		datatable_now[,4],datatable_now[,5],datatable_now[,6], stringsAsFactors=FALSE)
	datatable_to_use <- datatable_to_use2[-which(datatable_to_use2[,2]==1),]
	colnames(datatable_to_use) <- c("GEO accession", "Number of samples", "Species", "Title", "Study summary", "sra_runs")
	#save(datatable_to_use, file=paste0("/opt/raid10/genomics/naim/Thesis/geo/datatable_to_use.RData"), compress = FALSE)
	save(datatable_to_use, file=paste0("/opt/raid10/genomics/naim/myGithub/GRIN_dev/datatable_to_use.RData"), compress = FALSE)
}
datatable_to_use()	#Run this function every week to update the list
