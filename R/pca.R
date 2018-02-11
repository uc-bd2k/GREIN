
	pca <- function (data_mat, metadata, property, type) {

		exps <- data_mat[complete.cases(data_mat),]
		pseudoCount <- as.matrix(exps)
		groups <- metadata[,which(colnames(metadata)==property)] 
		b.PCA <- prcomp(t(pseudoCount),retx=TRUE,center=TRUE)
		pca.mat <- if (ncol(b.PCA$x)<5) {
						 b.PCA$x
					} else {
						cbind(b.PCA$x[,1],b.PCA$x[,2],b.PCA$x[,3],b.PCA$x[,4],b.PCA$x[,5])
					}
		load("data/allColors.rda")
		if (length(unique(groups))>105) {
			colcols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(groups)))
		} else {
			colcols <- allColors[seq_along(unique(groups))]
		}

		if (type=="2D") {
			pairsD3(pca.mat, group = groups, labels =c("PC1","PC2","PC3","PC4","PC5"), cex = 3,
			width = 1200, col = colcols[factor(as.character(groups),exclude=NULL)], theme = "colour",opacity = 1)
		} else {
			par3d(windowRect = c(10, 10, 800, 800),cex=1.2)
			if (ncol(b.PCA$x)<3) {
				plot3d(b.PCA$x[,1],b.PCA$x[,2],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="PC1", ylab = "PC2 ", size=1.5,
				zlab = "PC3 ", asp=.5, main=paste("3D PCA plot of the ", property, sep=""))
				rglwidget()
			} else {
				plot3d(b.PCA$x[,1],b.PCA$x[,2],b.PCA$x[,3],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="PC1", ylab = "PC2 ", size=1.5,
				zlab = "PC3 ", asp=.5, main=paste("3D PCA plot of the ", property, sep=""))
				rglwidget()

			}

		}
	}
	
