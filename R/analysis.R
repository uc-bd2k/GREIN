
########### Differential expression analysis using edgeR for two group comparison ###########
	analysis <- function(counts, metadata, genes, property, covariate=NULL, group1, group2) {
		
		if(group1==''|group2=='') {
			stop ("select the other group")
		} else {
			if(is.null(covariate)) {
				analysis_metadata <- metadata[which(metadata[,property]==group1 | metadata[,property]==group2),c(property,covariate),drop=FALSE]
				analysis_counts <- counts[, rownames(analysis_metadata)]
				group <- as.factor(analysis_metadata[,1])
				design <-model.matrix(~group)

				y <- DGEList(counts=analysis_counts, genes= genes, group=group)
				o <- order(rowSums(y$counts))
				y <- y[o,]
				keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
				y <- y[keep,]
				y$samples$lib.size <- colSums(y$counts)
				y <- calcNormFactors(y) 
				
				if (ncol(analysis_counts)>2) {
					y <- estimateCommonDisp(y)
					y <- estimateTagwiseDisp(y)
					fit <- exactTest(y, pair=levels(y$samples$group))	
					top_degs <- topTags(fit, n=nrow(fit$table))
					signaturesData_final <- top_degs$table[, -c(4,5)]
					colnames(signaturesData_final)[c(3,4)] <- c("Value_LogDiffExp", "Significance_pvalue")
					signaturesData_final$Value_LogDiffExp <- round(signaturesData_final$Value_LogDiffExp,3)
					signaturesData_final$Significance_pvalue <- round(signaturesData_final$Significance_pvalue,3)
					signaturesData_final[is.na(signaturesData_final)] <- "NA"
					rownames(signaturesData_final)=NULL
					return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts=y$counts))
				} else {
					fit <- exactTest(y, dispersion=0.2^2, pair=levels(y$samples$group))	
					top_degs <- topTags(fit, n=nrow(fit$table))
					signaturesData_final <- top_degs$table[, -c(4,5)]
					colnames(signaturesData_final)[c(3,4)] <- c("Value_LogDiffExp", "Significance_pvalue")
					signaturesData_final$Value_LogDiffExp <- round(signaturesData_final$Value_LogDiffExp,3)
					signaturesData_final$Significance_pvalue <- round(signaturesData_final$Significance_pvalue,3)
					signaturesData_final[is.na(signaturesData_final)] <- "NA"
					rownames(signaturesData_final)=NULL
					return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts=y$counts))
				}
			} else {
				group <- metadata[,property]
				covar <- metadata[,covariate]
				design <-model.matrix(~covar+group)

				y <- DGEList(counts=counts, genes= genes, group=group)
				o <- order(rowSums(y$counts))
				y <- y[o,]
				keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
				y <- y[keep,]
				y$samples$lib.size <- colSums(y$counts)
				y <- calcNormFactors(y) 
				y <- estimateDisp(y, design, robust=TRUE)
				
				contrast=rep(0, length(unique(metadata[,property])))
				contrast[which(group1==unique(metadata[,property]))]=1
				contrast[which(group2==unique(metadata[,property]))]=-1
				contrast=c(0,contrast)
				
				fit <- glmFit(y, design)
				lrt <-  glmLRT(fit, contrast=contrast)
				top_degs <- topTags(lrt, n=nrow(lrt$table))
				signaturesData_final <- top_degs$table[, -c(4,5)]
				colnames(signaturesData_final)[c(3,4)] <- c("Value_LogDiffExp", "Significance_pvalue")
				signaturesData_final$Value_LogDiffExp <- round(signaturesData_final$Value_LogDiffExp,3)
				signaturesData_final$Significance_pvalue <- round(signaturesData_final$Significance_pvalue,3)
				signaturesData_final[is.na(signaturesData_final)] <- "NA"
				rownames(signaturesData_final)=NULL
				return(list(signaturesData_final=signaturesData_final, disp=y$common.dispersion, counts=y$counts))


			
			}
		}
		
	}