
datatable_init <- function() {
	geo=list.files("/opt/raid10/genomics/naim/Thesis/geo/datasets/")
	datatable <- matrix(, nrow = length(geo), ncol = 5)
	for ( i in seq_along(geo)) {
		if (file.exists(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",geo[i],"/eset.RData"))) {
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",geo[i],"/eset.RData"))
			datatable[i,] <- cbind(
			geo[i],
			#paste0(
			#'<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
			#geo[i], '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'geoinfo\', 1)">', geo[i], '</a>'),
			nrow(pData(eset)),geo_study[geo_study$gse_id %in% geo[i], 4],
			geo_study[geo_study$gse_id %in% geo[i], 5], geo_study[geo_study$gse_id %in% geo[i], 6])
		} else {
			datatable[i,] <- cbind(rep(NA, 5))
		}
	}
	datatable_now <- as.data.frame(datatable[complete.cases(datatable),], stringsAsFactors=FALSE)
	datatable_to_use <- data.frame(datatable_now[,1], as.integer(datatable_now[,2]), datatable_now[,3],
		datatable_now[,4],datatable_now[,5], stringsAsFactors=FALSE)
	colnames(datatable_to_use) <- c("GEO Accession", "Number of samples", "Species", "Title", "Study summary")
	save(datatable_to_use, file=paste0("/opt/raid10/genomics/naim/Thesis/geo/datatable_to_use.RData"))
}
#load("data/geo_study.RData")
#datatable_init()	#Run this function every week to update the list
