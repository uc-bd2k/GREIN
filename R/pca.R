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


	pca <- function (b.PCA, metadata, property, type) {

		# exps <- data_mat[complete.cases(data_mat),]
		# pseudoCount <- exps[abs(rowSums(exps))>0,]

		# #pseudoCount <- as.matrix(dat)
		groups <- metadata[,which(colnames(metadata)==property)] 
		groups <- relevel(groups, unique(as.character(groups)[1]))
		# #pseudoCount <- removeBatchEffect(as.matrix(exps),batch=groups)
		# b.PCA <- prcomp(t(pseudoCount),retx=TRUE,center=TRUE, scale=TRUE)
		# pca.mat <- if (ncol(b.PCA$x)<5) {
						 # b.PCA$x
					# } else {
						# cbind(b.PCA$x[,1],b.PCA$x[,2],b.PCA$x[,3],b.PCA$x[,4],b.PCA$x[,5])
					# }
		load("data/allColors.rda")
		allColors = unique(allColors)[-c(6,7,11,13,15,16,19,20,22,24,25,28,35)]
		if (length(unique(groups))>22) {
			colcols <- colorRampPalette(brewer.pal(9,"Set1"))(length(levels(groups)))
			#colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(groups)))
		} else {
			colcols <- allColors[seq_along(levels(groups))]
		}
		d_tsne <- data.frame(b.PCA, group = groups)
		colnames(d_tsne) <- c("dim1", "dim2", "dim3", "group")

		if (type=="pca2D") {
			if (ncol(b.PCA)<5) {
				pairsD3(b.PCA, group = groups, labels =c("PC1","PC2","PC3","PC4","PC5"), cex = 3,
				width = 1200, col = colcols[factor(as.character(groups),exclude=NULL)], theme = "colour",opacity = 1)
			} else {
				pairsD3(cbind(b.PCA[,1],b.PCA[,2],b.PCA[,3],b.PCA[,4],b.PCA[,5]), group = groups, labels =c("PC1","PC2","PC3","PC4","PC5"), cex = 3,
				width = 1200, col = colcols, theme = "colour",opacity = 1)			
			}
		} else if (type=="pca3D") {
			par3d(windowRect = c(10, 10, 800, 800),cex=1.2)
			if (ncol(b.PCA)<3) {
				plot3d(b.PCA[,1],b.PCA[,2],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="PC1", ylab = "PC2 ", size=1.15,
				zlab = "PC3 ", asp=.5, main=paste("3D PCA plot of '", property,"'", sep=""))
				rglwidget()
			} else {
				plot3d(b.PCA[,1],b.PCA[,2],b.PCA[,3],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="PC1", ylab = "PC2 ", size=1.15,
				zlab = "PC3 ", asp=.5, main=paste0("3D PCA plot of '", property,"'"))
				rglwidget()

			}
		} else if (type=="tsne3D") {
			if (ncol(b.PCA)<3) {
				plot3d(b.PCA[,1],b.PCA[,2],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="dim1", ylab = "dim2 ", size=1.15,
				zlab = "dim3 ", asp=.5, main=paste("3D t-SNE plot of '", property,"'", sep=""))
				rglwidget()
			} else {
				plot_ly(d_tsne, x = ~dim1, y = ~dim2, z = ~dim3, color = ~group,  width = 900, colors = colcols) %>% add_markers() %>%
					layout(title = paste0("t-SNE plot of <span style='color:#337ab7'>", property,"</span>"),
							#margin = list(l = 50, r = 50, b = 50, t = 50),
							scene = list(aspectmode='cube',xaxis = list(title = 'dim1', showgrid = T),
								yaxis = list(title = 'dim2', showgrid = T),
								zaxis = list(title = 'dim3', showgrid = T),
								camera = list(eye = list(x = -1.65, y = 1.65, z = 1.65))
							)
						)			
			}
		} else if(type=="tsne2D"){
			plot_ly(d_tsne[,-3], x = ~dim1, y = ~dim2, color = ~group,  width = 900, colors = colcols) %>% add_markers() %>%
					layout(title = paste0("t-SNE plot of <span style='color:#337ab7'>", property,"</span>"))		
		} else {
			return()
		}	
	}
	
