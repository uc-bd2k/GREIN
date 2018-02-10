#' PCA function
#'
#' This function allows user to draw 2D and 3D PCA plot.
#' @param data_mat Data matrix where genes are row names and samples are column names. 
#' @param metadata Metadata associated with the data matrix where samples are the row names. The row names of the metadata need to be matched to the column names of the data matrix. 
#' @param property Property or the variable of interest from the metadata.
#' @export
#' @examples
#' load("data/EDS-1013eset.RData"); assign('eset', get("EDS-1013eset"));
#' PCA(data_mat=exprs(eset), metadata=pData(eset), property="ER", type="2D")


	pca <- function (data_mat, metadata, property, type) {

		exps <- data_mat[complete.cases(data_mat),]
		pseudoCount <- as.matrix(exps)
		groups <- metadata[,which(colnames(metadata)==property)] 
		#pseudoCount <- removeBatchEffect(as.matrix(exps),batch=groups)
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
	
