


	plotHeatFit <- function(){
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		rownames(exprs(eset)) <- paste0(fData(eset)[,2], " : ", fData(eset)[,1])
		complex_heatmap (exprs(eset), sub_metadata(pData(eset)), input$property, input$clustering_method, type="Fit in Screen")
	}
	
	output$fitscreen = renderUI({
		validate(
			need(!is.null(input$property), "please select a property")
		)
		plotOutput("HeatFit", height="auto")
	})
	output$HeatFit = renderPlot({
		validate(
			need(!is.null(input$property), "please select a property")
		)  
		print(plotHeatFit())
		
	}, 	#width=1200, 
		#height=1000
		height = reactive({
			if(!is.null(input$property)) {
				df <- metadata()[,input$property,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- input$property
				x=sapply(df, function(x) length(unique(x)))
				nrow = max(x)/2
				return((length(input$property)*25) + 650+(nrow*10))
			} else {
				return(600)
			}
		})
	)
	
	output$downloadFit = downloadHandler(
		filename = function() {
			datasetname <- sub(" :.*$", "", input$geo_acc)
			if (length(input$property)==1) {
				return(paste(datasetname,"_",input$property,"_heatmap_fit.png",sep=""))
			} else {
				return(paste(datasetname,"_","(",paste(input$property, collapse="+"),")","_heatmap_fit.png",sep=""))
			}
		},
		content = function(file) {
			png(file, width = 18, height = 16, units = "in",res=200)
			plotHeatFit()
			dev.off()
		}	
	)
	
	
	### Scrollable
	plotHeatScroll <- function(){
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		rownames(exprs(eset)) <- paste0(fData(eset)[,2], " : ", fData(eset)[,1])
		complex_heatmap (exprs(eset), sub_metadata(pData(eset)), input$property, input$clustering_method, type="Scrollable")
	}
	output$scroll = renderPlot({
		validate(
			need(!is.null(input$property), "please select a property")
		)
		print(plotHeatScroll())
	},  #width=1600, #height=12000
		height = reactive({
			exps= as.matrix(exprseset()) - rowMeans(as.matrix(exprseset()))
			medAbsDev<-apply(exps,1,function(x) median(abs(x)))
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
			perGene<-250/1000
			value <- max((length(input$property)*25)+(4*96),(length(input$property)*25)+(96*perGene*length(topGenes))) #1 inch=96 pixel
			print(value)
			return(value)
				
		})
	)
	
	output$downloadScroll = downloadHandler(
		filename = function() {
			datasetname <- sub(" :.*$", "", input$geo_acc)
			if (length(input$property)==1) {
				return(paste(datasetname,"_",input$property,"_heatmap_scroll.png",sep=""))
			} else {
				return(paste(datasetname,"_","(",paste(input$property, collapse="+"),")","_heatmap_scroll.png",sep=""))
			}
		},
		content = function(file) {
			height = reactive({
				exps= as.matrix(exprseset()) - rowMeans(as.matrix(exprseset()))
				medAbsDev<-apply(exps,1,function(x) median(abs(x)))
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

				perGene<-300/1000
				value <- max((length(input$property)*25)+(4*96),(length(input$property)*25)+(96*perGene*length(topGenes))) #1 inch=96 pixel
				print(value)
				return(value)
			})
			png(file, width = 3000, height = height(), res=175, units = "px")
			plotHeatScroll()
			dev.off()
		}	
	)