
########### Differential expression analysis using edgeR for two group comparison without covariate ###########

	analysis_2g_wocov <- function(counts, metadata, genes, property, group1, group2) {
		
		analysis_metadata <- metadata[which(metadata[,property]==group1 | metadata[,property]==group2), property,drop=FALSE]
		analysis_counts <- counts[, rownames(analysis_metadata), drop=F]
		rownames(analysis_counts) <- sapply(strsplit(rownames(analysis_counts)," : "), `[`, 1)
		group <- as.factor(analysis_metadata[,1])
		design <-model.matrix(~group)

		y <- DGEList(counts=analysis_counts, genes= genes, group=group)
		o <- order(rowSums(y$counts))
		y <- y[o,]
		keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
		y <- y[keep, , keep.lib.sizes=FALSE]
		y$samples$lib.size <- colSums(y$counts)
		y <- calcNormFactors(y) 
		
		if (ncol(analysis_counts)>2) {
			y <- estimateCommonDisp(y)
			y <- estimateTagwiseDisp(y)
			fit <- exactTest(y, pair=levels(y$samples$group))	
			top_degs <- topTags(fit, n=nrow(fit$table))
			signaturesData_final <- top_degs$table[, -c(4,5)]
			#colnames(signaturesData_final)[c(3,4)] <- c("Value_LogDiffExp", "Significance_pvalue")
			colnames(signaturesData_final)[c(3,4)] <- c("Log_FoldChange", "Adjusted_pvalue")
			signaturesData_final$Log_FoldChange <- round(signaturesData_final$Log_FoldChange,3)
			#signaturesData_final$Adjusted_pvalue <- format(signaturesData_final$Adjusted_pvalue, scientific=T)
			signaturesData_final$Adjusted_pvalue <- round(signaturesData_final$Adjusted_pvalue, 6)
			signaturesData_final[is.na(signaturesData_final)] <- "NA"
			rownames(signaturesData_final)=NULL
			counts_sig <- y$counts[match(signaturesData_final[,1], rownames(y$counts)),]
			rownames(counts) <- sapply(strsplit(rownames(counts)," : "), `[`, 1)
			counts_all <- counts[match(signaturesData_final[,1], rownames(counts)),]
			return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts_all=counts_all, counts_sig=counts_sig))
		} else {
			fit <- exactTest(y, dispersion=0.2^2, pair=levels(y$samples$group))	
			top_degs <- topTags(fit, n=nrow(fit$table))
			signaturesData_final <- top_degs$table[, -c(4,5)]
			colnames(signaturesData_final)[c(3,4)] <- c("Log_FoldChange", "Adjusted_pvalue")
			signaturesData_final$Log_FoldChange <- round(signaturesData_final$Log_FoldChange,3)
			#signaturesData_final$Adjusted_pvalue <- format(signaturesData_final$Adjusted_pvalue, scientific=T)
			signaturesData_final$Adjusted_pvalue <- round(signaturesData_final$Adjusted_pvalue, 6)
			signaturesData_final[is.na(signaturesData_final)] <- "NA"
			rownames(signaturesData_final)=NULL
			counts_sig <- y$counts[match(signaturesData_final[,1], rownames(y$counts)),]
			rownames(counts) <- sapply(strsplit(rownames(counts)," : "), `[`, 1)
			counts_all <- counts[match(signaturesData_final[,1], rownames(counts)),]
			return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts_all=counts_all, counts_sig=counts_sig))
		}
	}
	

########### Differential expression analysis using edgeR for two group comparison with covariate ###########

	analysis_2g_withcov <- function(counts, metadata, genes, property, covariate, group1, group2) {
		
		if(is.null(covariate)){
			warning("Please select a covariate")
		} else {
			analysis_metadata <- metadata[which(metadata[,property]==group1 | metadata[,property]==group2),c(property,covariate),drop=FALSE]
			analysis_counts <- counts[, rownames(analysis_metadata), drop=F]
			rownames(analysis_counts) <- sapply(strsplit(rownames(analysis_counts)," : "), `[`, 1)
			group <- analysis_metadata[,property, drop=F]
			covar <- analysis_metadata[,covariate, drop=F]
			design_data <- cbind(covar, group)
			design_data[] <- lapply(design_data, factor)
			design <- model.matrix(~., data=design_data)
			#nonEstimable(design)
			#is.fullrank(design)
			ne <- nonEstimable(design)
			if(!is.null(ne)) {
				stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n", paste(ne, collapse = " ")))
			} else {
				y <- DGEList(counts=analysis_counts, genes= genes)
				o <- order(rowSums(y$counts))
				y <- y[o,]
				keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
				y <- y[keep, , keep.lib.sizes=FALSE]
				y$samples$lib.size <- colSums(y$counts)
				y <- calcNormFactors(y) 
				
				y <- estimateGLMCommonDisp(y,design)
				y <- estimateGLMTrendedDisp(y,design)
				y <- estimateGLMTagwiseDisp(y,design)
				fit <- glmFit(y, design)
				lrt <-  glmLRT(fit)

				# y <- estimateDisp(y, design, robust=TRUE)
				# fit <- glmQLFit(y, design, robust=TRUE)
				# lrt <- glmQLFTest(fit)

				# contrast=rep(0, length(unique(analysis_metadata[,property])))
				# contrast[which(group1==unique(analysis_metadata[,property]))]=1
				# contrast[which(group2==unique(analysis_metadata[,property]))]=-1
				# contrast=c(0,contrast)
								
				top_degs <- topTags(lrt, n=nrow(lrt$table))
				signaturesData_final <- top_degs$table[, -c(4,5,6)]
				colnames(signaturesData_final)[c(3,4)] <- c("Log_FoldChange", "Adjusted_pvalue")
				signaturesData_final$Log_FoldChange <- round(signaturesData_final$Log_FoldChange,3)
				#signaturesData_final$Adjusted_pvalue <- format(signaturesData_final$Adjusted_pvalue, scientific=T)
				signaturesData_final$Adjusted_pvalue <- round(signaturesData_final$Adjusted_pvalue, 6)
				signaturesData_final[is.na(signaturesData_final)] <- "NA"
				rownames(signaturesData_final)=NULL				
				counts_sig <- y$counts[match(signaturesData_final[,1], rownames(y$counts)),]
				rownames(counts) <- sapply(strsplit(rownames(counts)," : "), `[`, 1)
				counts_all <- counts[match(signaturesData_final[,1], rownames(counts)),]
				return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts_all=counts_all, counts_sig=counts_sig))
			}
		}
	}


########### Differential expression analysis using edgeR for multi group comparison without covariate ###########

	analysis_mg_wocov <- function(counts, metadata, genes, property) {
			
		analysis_metadata <- metadata[,property,drop=FALSE]
		analysis_counts <- counts[, rownames(analysis_metadata), drop=F]
		rownames(analysis_counts) <- sapply(strsplit(rownames(analysis_counts)," : "), `[`, 1)
		group <- as.factor(analysis_metadata[,1])
		design <- model.matrix(~group)

		y <- DGEList(counts=analysis_counts, genes= genes, group=group)
		o <- order(rowSums(y$counts))
		y <- y[o,]
		keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
		y <- y[keep, , keep.lib.sizes=FALSE]
		y$samples$lib.size <- colSums(y$counts)
		y <- calcNormFactors(y) 
		
		y <- estimateGLMCommonDisp(y,design)
		y <- estimateGLMTrendedDisp(y,design)
		y <- estimateGLMTagwiseDisp(y,design)
		fit <- glmFit(y,design)
		lrt <- glmLRT(fit,coef=2:length(levels(group)))
		
		# y <- estimateDisp(y, design, robust=TRUE)
		# fit <- glmQLFit(y, design, robust=TRUE)
		# lrt <- glmQLFTest(fit, coef=2:length(levels(group)))

		top_degs <- topTags(lrt, n=nrow(lrt$table))
		signaturesData_final1 <- top_degs$table[, !(names(top_degs$table) %in% c('logCPM', 'LR', 'PValue'))]
		colnames(signaturesData_final1)[length(colnames(signaturesData_final1))] <- c("Adjusted_pvalue")
		x <- round(signaturesData_final1[!(names(signaturesData_final1) %in% c("ID_probe", "Name_GeneSymbol", "Adjusted_pvalue"))], 3)
		signaturesData_final <- data.frame("ID_probe"=signaturesData_final1$ID_probe, "Name_GeneSymbol"=signaturesData_final1$Name_GeneSymbol,x, "Adjusted_pvalue"=signaturesData_final1$Adjusted_pvalue)
		#signaturesData_final$Adjusted_pvalue <- format(signaturesData_final$Adjusted_pvalue, scientific=T)
		signaturesData_final$Adjusted_pvalue <- round(signaturesData_final$Adjusted_pvalue, 6)
		signaturesData_final[is.na(signaturesData_final)] <- "NA"
		rownames(signaturesData_final)=NULL
		counts_sig <- y$counts[match(signaturesData_final[,1], rownames(y$counts)),]
		rownames(counts) <- sapply(strsplit(rownames(counts)," : "), `[`, 1)
		counts_all <- counts[match(signaturesData_final[,1], rownames(counts)),]
		return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts_all=counts_all, counts_sig=counts_sig))
	}
	
########### Differential expression analysis using edgeR for multi group comparison with covariate ###########

	analysis_mg_withcov <- function(counts, metadata, genes, property, covariate) {
			
		if(is.null(covariate)){
			warning("Please select a covariate")
		} else {
			analysis_metadata <- metadata[,c(property,covariate),drop=FALSE]
			analysis_counts <- counts[, rownames(analysis_metadata), drop=F]
			rownames(analysis_counts) <- sapply(strsplit(rownames(analysis_counts)," : "), `[`, 1)
			group <- analysis_metadata[,property, drop=F]
			covar <- analysis_metadata[,covariate, drop=F]
			design_data <- cbind(covar, group)
			design_data[] <- lapply(design_data, factor)
			design <- model.matrix(~., data=design_data)
			ne <- nonEstimable(design)
			if(!is.null(ne)) {
				stop(paste("Design matrix not of full rank. The following coefficients not estimable:\n", paste(ne, collapse = " ")))
			} else {
				y <- DGEList(counts=analysis_counts, genes= genes)
				o <- order(rowSums(y$counts))
				y <- y[o,]
				keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
				y <- y[keep, , keep.lib.sizes=FALSE]
				y$samples$lib.size <- colSums(y$counts)
				y <- calcNormFactors(y) 

				y <- estimateGLMCommonDisp(y,design)
				y <- estimateGLMTrendedDisp(y,design)
				y <- estimateGLMTagwiseDisp(y,design)
				fit <- glmFit(y,design)
				lrt <- glmLRT(fit,coef=2:length(levels(group)))
				
				# y <- estimateDisp(y, design, robust=TRUE)
				# fit <- glmQLFit(y, design, robust=TRUE)
				# lrt <- glmQLFTest(fit, coef=2:length(levels(group)))

				top_degs <- topTags(lrt, n=nrow(lrt$table))
				signaturesData_final1 <- top_degs$table[, !(names(top_degs$table) %in% c('logCPM', 'LR', 'PValue'))]
				colnames(signaturesData_final1)[length(colnames(signaturesData_final1))] <- c("Adjusted_pvalue")
				x <- round(signaturesData_final1[!(names(signaturesData_final1) %in% c("ID_probe", "Name_GeneSymbol", "Adjusted_pvalue"))], 3)
				signaturesData_final <- data.frame("ID_probe"=signaturesData_final1$ID_probe, "Name_GeneSymbol"=signaturesData_final1$Name_GeneSymbol,x, "Adjusted_pvalue"=signaturesData_final1$Adjusted_pvalue)
				#signaturesData_final$Adjusted_pvalue <- format(signaturesData_final$Adjusted_pvalue, scientific=T)
				signaturesData_final$Adjusted_pvalue <- round(signaturesData_final$Adjusted_pvalue, 6)
				signaturesData_final[is.na(signaturesData_final)] <- "NA"
				rownames(signaturesData_final)=NULL
				counts_sig <- y$counts[match(signaturesData_final[,1], rownames(y$counts)),]
				rownames(counts) <- sapply(strsplit(rownames(counts)," : "), `[`, 1)
				counts_all <- counts[match(signaturesData_final[,1], rownames(counts)),]
				return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts_all=counts_all, counts_sig=counts_sig))
			}
		}
	}