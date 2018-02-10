#' Complex Heatmap function
#'
#' This function allows user to draw complex heatmap with options: "Fit in Screen" and "Scrollable" along with different clustering methods using ComplexHeatmap package.
#' @param data_mat Data matrix where genes are row names and samples are column names. 
#' @param metadata Metadata associated with the data matrix where samples are the row names. The row names of the metadata need to be matched to the column names of the data matrix. 
#' @param property Property or the variable of interest from the metadata.
#' @param clustering_method Which method to use for clustering. Methods are: "precomp_GIMM", "Clustering by groups", "Pearson Correlation", "Kendall", "Spearman", and "Euclidean". 
#' "precomp_GIMM" is the precomputed GIMM clustering that is saved in the working directory. Default clustering method is "Pearson Correlation".
#' @param type Type of heatmap. There are two types: Fit in Screen and Scrollable. Default is "Fit in Screen". 
#' @export 
#' @examples
#' load("data/EDS-1013eset.RData"); assign('eset', get("EDS-1013eset"));
#' complex_heatmap(property="ER",data_mat=exprs(eset), metadata=pData(eset), clustering_method="Pearson Correlation", type="Fit in Screen")

	complex_heatmap <- function (data_mat, metadata, property, clustering_method="Pearson Correlation", type="Fit in Screen") {

		if(!is.null(property)) {
			
			exps= as.matrix(data_mat) - rowMeans(as.matrix(data_mat))
			medAbsDev <-apply(exps,1,function(x) median(abs(x)))
			topGenes= function(exps, medAbsDev) {
				if (dim(exps)[1]>= 1000 & dim(exps)[1]<=1500) {
					topGenes<-order(medAbsDev,decreasing=T)[1:1500]
				} else if (dim(exps)[1]> 1500 ) {
					topGenes<-order(medAbsDev,decreasing=T)[1:1000]
				} else if (dim(exps)[1]< 1000) {
					topGenes<-order(medAbsDev,decreasing=T)
				} else {
					topGenes=medAbsDev
				}
				return(topGenes)
			}
			topGenes=topGenes(exps, medAbsDev)
			topGenes= topGenes[!is.na(topGenes)]

			expr_data_in<- data.matrix(exps[topGenes,])
			expr_data_in[expr_data_in>2]<-2
			expr_data_in[expr_data_in<(-2)]<-(-2)
			
			metadata <- metadata[,property,drop=FALSE]
			metadata[is.na(metadata)] <- "NA"
			metadata[] <- lapply( metadata, factor)
			colnames(metadata) <- property
			x=sapply(metadata, function(x) length(unique(x)))
			nrow = round(max(x)/2)+1

			load("data/allColors.rda")
			colr= vector(mode="list", length=ncol(metadata))
			names(colr)= colnames(metadata)
			for(i in seq_len(ncol(metadata))) {
				if (length(unique(as.character(metadata[,i])))>105) {
					colr[[i]]=setNames(colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(as.character(metadata[,i])))), unique(as.character(metadata[,i])))
				} else {
					colr[[i]]=setNames(allColors[seq_along(unique(as.character(metadata[,i])))], unique(as.character(metadata[,i])))
				}
			}

			rside_length= function(property) {
				if(length(property)==1 && nchar(property)<=3 ) {
					return(15)
				}
				else if (length(property)==1 && nchar(property)>3 ) {
					return(nchar(property)*4)
				}
				else if (length(property)>1 && max(nchar(property))<=3 ) {
					return(max(nchar(property))*10)
				}
				else if (length(property)>1 && max(nchar(property))>3 ) {
					return(max(nchar(property))*4)
				}
				else {
					return(NULL)
				}
			}
			
			heat <- function(data_mat,metadata,cluster_columns,cluster_rows, top_annotation,col_fontsize,clustering_distance_columns="pearson", clustering_distance_rows="pearson") {
					
				if (type=="Fit in Screen") {
					ht <- Heatmap(data_mat,  name = "Expression",  col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
						cluster_columns = cluster_columns, cluster_rows = cluster_rows , 
						clustering_distance_columns = clustering_distance_columns, clustering_distance_rows = clustering_distance_rows,
						show_row_names = FALSE, show_column_names = if(length(colnames(data_mat))<100) {TRUE} else {FALSE},
						column_names_max_height=unit((4/10)*max(nchar(colnames(data_mat))), "cm"),column_names_gp= gpar(fontsize = col_fontsize),
						row_dend_reorder=FALSE, top_annotation = top_annotation, top_annotation_height= unit(0.7*ncol(metadata), "cm"),
						heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
						legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht,padding = unit(c(1,10,4,rside_length(property)), "mm"), heatmap_legend_side = "top", annotation_legend_side = "bottom")
				} else {
					ht <- Heatmap(data_mat,  name = "Expression",  col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
						cluster_columns = cluster_columns, cluster_rows = cluster_rows , 
						clustering_distance_columns = clustering_distance_columns, clustering_distance_rows = clustering_distance_rows,
						show_row_names = TRUE, row_names_side = "right",row_names_gp = gpar(fontsize = 10), row_names_max_width = unit(8, "cm"),
						show_column_names = if(length(colnames(data_mat))<100) {TRUE} else {FALSE},
						column_names_max_height=unit((4/10)*max(nchar(colnames(data_mat))), "cm"),column_names_gp= gpar(fontsize = col_fontsize),
						row_dend_reorder=FALSE, top_annotation = top_annotation, top_annotation_height= unit(0.7*ncol(metadata), "cm"),
						heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
						legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht,heatmap_legend_side = "top", annotation_legend_side = "right")
				}
				for(an in colnames(metadata)) {
					decorate_annotation(an, {
						# annotation names on the right
						grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 14))
					})
				}
			}
			
			gimm_rda <- list.files(pattern = "*.rda$", recursive = TRUE)
			if (clustering_method=="precomp_GIMM" & length(gimm_rda)==2) {
			
				load(gimm_rda[1])
				load(gimm_rda[2])
				toCluster_R <-gimmOutRow$clustData

				forHeatmap_pgimm <-data.matrix(toCluster_R[,-(1:2)])
				forHeatmap_pgimm[forHeatmap_pgimm>2]<-2
				forHeatmap_pgimm[forHeatmap_pgimm<(-2)]<-(-2)
				
				if(abs(12-(ncol(forHeatmap_pgimm)/30))<2) {
					col_fontsize <- abs(12-(ncol(forHeatmap_pgimm)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(forHeatmap_pgimm)/30))
				}
				data_mat <- forHeatmap_pgimm
				cluster_columns = as.dendrogram(gimmOutCol$hGClustData)
				cluster_rows = as.dendrogram(gimmOutRow$hGClustData)
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1}))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation,col_fontsize,
						clustering_distance_columns="pearson", clustering_distance_rows="pearson")
				
			} else if (clustering_method=="Clustering by groups") {
				metadata_cg <- metadata[do.call(order, c(data.frame(metadata[,1:ncol(metadata)]))),,drop=FALSE]
				forHeatmap_cg <- expr_data_in[,rownames(metadata_cg)]

				top_annotation <- HeatmapAnnotation(metadata_cg, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1} )) 
				
				if(abs(12-(ncol(forHeatmap_cg)/30))<2) {
					col_fontsize <- abs(12-(ncol(forHeatmap_cg)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(forHeatmap_cg)/30))
				}
				data_mat <- forHeatmap_cg
				cluster_columns = FALSE
				cluster_rows = TRUE
									
				heat(data_mat,metadata_cg,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="pearson", clustering_distance_rows="pearson")

			} else if (clustering_method=="Pearson Correlation") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1}))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="pearson", clustering_distance_rows="pearson")
				
			} else if (clustering_method=="Spearman") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1}))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="spearman", clustering_distance_rows="spearman")
			
			} else if (clustering_method=="Kendall") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1}))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="kendall", clustering_distance_rows="kendall")
						
			} else if (clustering_method=="Euclidean") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"),
				col = colr, annotation_legend_param=list(title_gp = gpar(fontsize = 14), 
				nrow= if(type=="Fit in Screen") {nrow} else {NULL}, ncol= if(type=="Fit in Screen") {NULL} else {1}))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="euclidean", clustering_distance_rows="euclidean")
						
			} else {
				warning("specify an available clustering method")
			}
		}
		else {
			warning("select a property first")
		}
	}

