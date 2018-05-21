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

	complex_heatmap_sig <- function (data_mat, metadata, property, n_genes, clustering_method, type="Fit in Screen") {
		
		if(!is.null(property)) {
			
			dge <- DGEList(counts=data_mat)
			dge_cpm <- cpm(dge, log=TRUE)
			exps= as.matrix(dge_cpm) - rowMeans(as.matrix(dge_cpm), na.rm=T)
			expr_data<- data.matrix(exps)
			expr_data_in <- expr_data[1:n_genes,,drop=F]
			
			if(n_genes>1) {
				if(length(which(apply(expr_data_in,1, function(x) length(unique(x)))==1))>0) {
					expr_data_in <- expr_data_in[-which(apply(expr_data_in,1, function(x) length(unique(x)))==1),, drop=F]
				} else if (length(which(apply(expr_data_in,2, function(x) length(unique(x)))==1))>0) {
					expr_data_in <- expr_data_in[,-which(apply(expr_data_in,2, function(x) length(unique(x)))==1),drop=F]
				} else {
					expr_data_in <- expr_data_in
				}
			} else {
				expr_data_in <- expr_data_in
			}
			
			if((sum(apply(expr_data_in,2,function(x) sum(x>-2 & x<2)))/length(expr_data_in))>=0.5) {
				expr_data_in[expr_data_in>2]<-2
				expr_data_in[expr_data_in<(-2)]<-(-2)
			} else {
				expr_data_in[expr_data_in>4]<-4
				expr_data_in[expr_data_in<(-4)]<-(-4)
			}
			
			metadata <- metadata[,property,drop=FALSE]
			metadata[is.na(metadata)] <- "NA"
			metadata[] <- lapply( metadata, factor)
			colnames(metadata) <- property
			for(i in 1:ncol(metadata)){
				metadata[,i] <- relevel(metadata[,i], unique(as.character(metadata[,i])[1]))
			}
			x=sapply(metadata, function(x) length(unique(x)))
			nrow = if(max(x)>10) {round(max(x)/2)+1} else {max(x)}
			n_chr = max(apply(metadata,2,nchar))
			
			load("data/allColors.rda")
			colr= vector(mode="list", length=ncol(metadata))
			names(colr)= colnames(metadata)
			allColors = unique(allColors)[-c(6,7,11,13,15,16,19,20,22,24,25,28,35)]
			for(i in seq_len(ncol(metadata))) {

				if (length(unique(as.character(metadata[,i])))>22) {
					if("NA" %in% unique(as.character(metadata[,i]))) {
						mis_num = grep("NA", unique(as.character(metadata[,i])))
						colr[[i]]=setNames(c(colorRampPalette(brewer.pal(9,"Set1"))(length(unique(as.character(metadata[,i]))[-mis_num])), "gray"), 
											c(unique(as.character(metadata[,i]))[-mis_num], "NA"))
					} else {
						colr[[i]]=setNames(colorRampPalette(brewer.pal(9,"Set1"))(length(levels(as.factor(metadata[,i])))), 
											levels(as.factor(metadata[,i])))
					}
				} else {
					if("NA" %in% unique(as.character(metadata[,i]))) {
						mis_num = grep("NA", unique(as.character(metadata[,i])))
						#colr[[i]][mis_num] = allColors[1]
						colr[[i]]=setNames(c(allColors[seq_along(unique(as.character(metadata[,i]))[-mis_num])], "gray"), 
											c(unique(as.character(metadata[,i]))[-mis_num], "NA"))
					} else {
						#colr[[i]]=setNames(allColors[seq_along(sort(unique(as.character(metadata[,i]))))], sort(unique(as.character(metadata[,i]))))
						colr[[i]]=setNames(allColors[seq_along(levels(as.factor(metadata[,i])))], levels(as.factor(metadata[,i])))
					}
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
					return(max(nchar(property))*3)
				}
				else {
					return(NULL)
				}
			}
			
			heat <- function(data_mat,metadata,cluster_columns,cluster_rows, top_annotation,col_fontsize,clustering_distance_columns="pearson", clustering_distance_rows="pearson") {
					
				if (type=="Fit in Screen") {
					ht <- Heatmap(data_mat,  name = "Expression",  col = colorRamp2(c(min(data_mat),0,max(data_mat)), c("blue", "black","yellow")),
						cluster_columns = cluster_columns, cluster_rows = cluster_rows , column_dend_reorder = FALSE,
						clustering_distance_columns = clustering_distance_columns, clustering_distance_rows = clustering_distance_rows,
						show_row_names = FALSE, show_column_names = if(length(colnames(data_mat))<100) {TRUE} else {FALSE},
						column_names_max_height=unit((4/10)*max(nchar(colnames(data_mat))), "cm"),column_names_gp= gpar(fontsize = col_fontsize),
						row_dend_reorder=FALSE, top_annotation = top_annotation, top_annotation_height= unit(0.7*ncol(metadata), "cm"),
						heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
						legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht,padding = unit(c(1,10,4,rside_length(property)), "mm"), heatmap_legend_side = "top", annotation_legend_side = "top")
				} else {
					ht <- Heatmap(data_mat,  name = "Expression",  col = colorRamp2(c(min(data_mat),0,max(data_mat)), c("blue", "black","yellow")),
						cluster_columns = cluster_columns, cluster_rows = cluster_rows , column_dend_reorder = FALSE,
						clustering_distance_columns = clustering_distance_columns, clustering_distance_rows = clustering_distance_rows,
						show_row_names = TRUE, row_names_side = "right",row_names_gp = gpar(fontsize = 10), row_names_max_width = unit(8, "cm"),
						show_column_names = if(length(colnames(data_mat))<100) {TRUE} else {FALSE},
						column_names_max_height=unit((4/10)*max(nchar(colnames(data_mat))), "cm"),column_names_gp= gpar(fontsize = col_fontsize),
						row_dend_reorder=FALSE, top_annotation = top_annotation, top_annotation_height= unit(0.7*ncol(metadata), "cm"),
						heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
						legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht,heatmap_legend_side = "top", annotation_legend_side = "top")
				}
				for(an in colnames(metadata)) {
					decorate_annotation(an, {
						# annotation names on the right
						grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 12))
					})
				}
			}
			
			if (clustering_method=="Group by properties") {
				metadata_cg <- metadata[do.call(order, c(data.frame(metadata[,1:ncol(metadata)]))),,drop=FALSE]
				forHeatmap_cg <- expr_data_in[,rownames(metadata_cg), drop=F]

				top_annotation <- HeatmapAnnotation(metadata_cg, which="column", width = unit(1,"mm"),	col = colr, 
				annotation_legend_param=list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 9),
				title_position = "topcenter", width = unit(5, "mm"),
				#nrow= if(type=="Fit in Screen") {nrow} else {nrow}, ncol= if(type=="Fit in Screen") {NULL} else {1} 
				nrow= if(n_chr>30) {NULL} else {nrow}, ncol= if(n_chr>30) {1} else {NULL}				
				)) 
				
				if(abs(12-(ncol(forHeatmap_cg)/30))<2) {
					col_fontsize <- abs(12-(ncol(forHeatmap_cg)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(forHeatmap_cg)/30))
				}
				data_mat <- forHeatmap_cg
				cluster_columns = FALSE
				cluster_rows = TRUE
									
				heat(data_mat,metadata_cg,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns=FALSE, clustering_distance_rows="pearson")

			} else if (clustering_method=="Pearson correlation") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"), col = colr, 
				annotation_legend_param=list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 9), 
				title_position = "topcenter", width = unit(5, "mm"), 
				#nrow= if(type=="Fit in Screen") {nrow} else {nrow}, ncol= if(type=="Fit in Screen") {NULL} else {1}
				nrow= if(n_chr>30) {NULL} else {nrow}, ncol= if(n_chr>30) {1} else {NULL}
				))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns="pearson", clustering_distance_rows="pearson")
				
			} else if (clustering_method=="Euclidean distance") {
			
				if(abs(12-(ncol(expr_data_in)/30))<2) {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(expr_data_in)/30))
				}
				data_mat <- expr_data_in
				cluster_columns = TRUE
				cluster_rows = TRUE
				
				top_annotation <- HeatmapAnnotation(metadata, which="column", width = unit(1,"mm"), col = colr, 
				annotation_legend_param=list(title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 9), 
				title_position = "topcenter", width = unit(5, "mm"), 
				#nrow= if(type=="Fit in Screen") {nrow} else {nrow}, ncol= if(type=="Fit in Screen") {NULL} else {1}
				nrow= if(n_chr>30) {NULL} else {nrow}, ncol= if(n_chr>30) {1} else {NULL}
				))

				heat(data_mat,metadata,cluster_columns,cluster_rows, top_annotation, col_fontsize,
						clustering_distance_columns=if(cluster_columns) {"euclidean"} else {FALSE}, clustering_distance_rows="euclidean")
						
			} else {
				warning("specify an available clustering method")
			}
		}
		else {
			warning("select a property first")
		}
	}

