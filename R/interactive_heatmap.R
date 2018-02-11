
interactive_heatmap <- function (data_mat, metadata, property, n_genes, clustering_method) {

	if(!is.null(property)) {
		dge <- DGEList(counts=data_mat)
		dge <- calcNormFactors(dge)
		dge_cpm <- cpm(dge, log=TRUE)
		exps= as.matrix(dge_cpm) - rowMeans(as.matrix(dge_cpm), na.rm=T)
		medAbsDev <-apply(exps,1,function(x) median(abs(x)))
		topGenes <- order(medAbsDev,decreasing=T)
		topGenes= topGenes[!is.na(topGenes)]

		expr_data<- data.matrix(exps[topGenes,,drop=F])
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
		
		expr_data_in[expr_data_in>2]<-2
		expr_data_in[expr_data_in<(-2)]<-(-2)
		
		metadata <- metadata[,property,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- property
		for(i in 1:ncol(metadata)){
			metadata[,i] <- relevel(metadata[,i], unique(as.character(metadata[,i])[1]))
		}
		x=sapply(metadata, function(x) length(unique(x)))
		nrow = round(max(x)/2)+1
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
		
		if (clustering_method=="Pearson correlation") {
			hr <- hclust(as.dist(1-cor(t(expr_data_in), method="pearson")))
			hc <- hclust(as.dist(1-cor(expr_data_in, method="pearson")))
				
			main_heatmap(expr_data_in, name = "Expression", colors=c("blue","black","yellow")) %>%
			#add_row_clustering(method = "hclust", clust_dist=hr(expr_data_in)) %>%
			add_col_annotation(metadata,  side = "top", colors=colr) %>%
			#add_col_clustering(method = "hclust", clust_dist=hc(expr_data_in), side = "top") %>%
			add_col_labels(side="bottom") %>%
			add_row_dendro(hr) %>%
			add_col_dendro(hc)

		} else if(clustering_method=="Euclidean distance"){
		
			main_heatmap(expr_data_in, name = "Expression", colors=c("blue","black","yellow")) %>%
			add_row_clustering(method = "hclust") %>%
			add_col_annotation(metadata,  side = "top", colors=colr) %>%
			add_col_clustering(method = "hclust", side = "top") %>%
			add_col_labels(side="bottom")

		} else if(clustering_method=="Group by properties"){
			metadata_cg <- metadata[do.call(order, c(data.frame(metadata[,1:ncol(metadata)]))),,drop=FALSE]
			forHeatmap_cg <- expr_data_in[,rownames(metadata_cg), drop=F]

			main_heatmap(forHeatmap_cg, name = "Expression", colors=c("blue","black","yellow")) %>%
			add_row_clustering(method = "hclust") %>%
			add_col_annotation(metadata_cg,  side = "top", colors=colr) %>%
			add_col_labels(side="bottom")
		
		} else {
			warning("specify an available clustering method")
		}
	
	} else {
		warning("select a property first")
	}
}