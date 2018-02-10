
load("data/geo_study.RData")
source("R/sub_metadata.R", local = TRUE)

shinyServer(function(input, output, session) {
	
	load("/opt/raid10/genomics/naim/Thesis/geo/datatable_to_use.RData")
	datafrm <- data.frame(datatable_to_use ,check.names = FALSE)
	datafr <- datafrm[order(datafrm$Species),]
	#apply(data.matrix(datatable_to_use[,1]), 1,function(x) unlist(strsplit(unlist(strsplit(x, split='>', fixed=TRUE))[2], split='<', fixed=TRUE))[1])
	#Analyze= paste0('Analyze : ',datatable_to_use[,1][grep("GSE", datatable_to_use[,1])])
	#<button><b>Analyze</b></button>
	output$datatable <- DT::renderDataTable (
		datafr,
		selection=list(mode="single", target="row"),
		escape = which(colnames(datafr) %in% c("GEO Accession", "Number of samples", "Species", "Title", "Study summary")),
		rownames = FALSE, filter = 'top', 
		options = list(
			columnDefs = list(
				list(width = '45%', targets =4
					#render=JS(
					# "function(data, type, full) {",
					# "return '<a class=\"btn btn-primary btn-sm\" onclick= href=#/' + '>' + 'Analyze' + '</a>';",
					# "}"  
					#)
				)
			),
			pageLength = 10
		)
	# %>% formatStyle(0, cursor = 'pointer') 
    )
	
	observe({
		hide(selector = "#explorgeo li a[data-value=explore]")
	})
	observeEvent(input$datatable_rows_selected,{
		if (!is.null(input$datatable_rows_selected)) {
			shinyjs::show(selector = "#explorgeo li a[data-value=explore]")
		} else {
			shinyjs::hide(selector = "#explorgeo li a[data-value=explore]")
		}
	})

	observeEvent(input$datatable_rows_selected, {
		info = input$datatable_rows_selected
		if (is.null(info)) return()
		updateTabsetPanel(session, "explorgeo", selected = "explore")
		updateTextInput(session, 'geo_acc', 
		#value = apply(data.matrix(datafr[info,1]), 1,
		#	function(x) unlist(strsplit(unlist(strsplit(x, split='>', fixed=TRUE))[2], split='<', fixed=TRUE))[1]))
		value= datafr[info,1])
	})
	
	#sub(" :.*$", "", input$geo_acc)
	
	############### organize eset data for following use ##############
	exprseset <- reactive({
		#load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		assign('eset', get(load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))))
		eset <- get('eset')
		rownames(exprs(eset)) <- paste0(fData(eset)[,2], " : ", fData(eset)[,1])
		return(exprs(eset))
	})
	metadata <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		return(sub_metadata(pData(eset)))
	})
	
	######################## Description ########################
	rows <- reactive ({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		return(ncol(exprs(eset)))
		ncol(exprseset())
	})
	observeEvent(input$geo_acc,{
		output$geo_summary <- renderText(
			paste(
				paste0(h4("Study link: "), paste0('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
					input$geo_acc, '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'geoinfo\', 1)">', input$geo_acc, '</a>')),
				h4(paste0("No. of samples : ",rows())), 
				paste(h4("Title: "),geo_study[which(geo_study$gse_id==input$geo_acc),5], sep="\n"),
				paste(h4("Summary: "),geo_study[which(geo_study$gse_id==input$geo_acc),6], sep="\n"),
			sep="\n")
		)
	})
  
	
	######################## metadata column ########################


	output$metadata_column <- renderUI({	
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		validate(
			need(any(input$tab2=="metadata"),""),
			need(!is.null(dim(sub_metadata(pData(eset)))), "Either the variables do not contain meaningful information or they contain single value. Please download the metadata table to take a look.")
		)
		metadata <- sub_metadata(pData(eset))
		col_names <- colnames(metadata)
		checkboxGroupInput("property", h4("Select properties"), choices  = col_names, 
							selected=col_names[(length(col_names)-1):length(col_names)])
	})
	
	output$full_metacol <- renderUI({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		validate(
			need(any(input$tab2=="metadata"),"")
		)
		#checkboxInput("full_meta", label = h4("Show full metadata"), value = FALSE)
		radioButtons("full_meta", label = h4("Show full metadata"), c("yes","no"), selected="no")
	})
	
	output$metadata <- DT::renderDataTable({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		if(input$full_meta=="yes"){
			datatable(pData(eset),selection = 'none',options = list(pageLength = if(nrow(pData(eset))>=10) {10} else {nrow(pData(eset))},
								scrollX = TRUE, dom = 't'),
				rownames= 
				if(ncol(pData(eset))<5) {
					apply(data.matrix(rownames(pData(eset))),1, function(x) paste0(unlist(strsplit(x,"_"))[1], "_",unlist(strsplit(x,"_"))[2]))
				} else {
					FALSE
				}
			)
		} else {
			validate(
				need(!is.null(dim(sub_metadata(pData(eset)))), "Either the variables do not contain meaningful information or they contain single value. Please download the metadata table to take a look.")
			)
			datatable(sub_metadata(pData(eset))[, input$property, drop = FALSE],selection = 'none',options = list(pageLength = if(nrow(sub_metadata(pData(eset)))>=10) {10} else {nrow(sub_metadata(pData(eset)))}),
				rownames= 
				if(ncol(sub_metadata(pData(eset)))<5) {
					apply(data.matrix(rownames(sub_metadata(pData(eset)))),1, function(x) paste0(unlist(strsplit(x,"_"))[1], "_",unlist(strsplit(x,"_"))[2]))
				} else {
					FALSE
				}
			)
		}
		
	})
	output$downloadmeta <- downloadHandler(
		filename = function() { 
			if(input$full_meta=="yes"){
				paste(input$geo_acc, '_full_metadata.csv', sep='') 
			} else {
				paste(input$geo_acc, '_filtered_metadata.csv', sep='') 
			}
		},
		content = function(file) {
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
			if(input$full_meta=="yes"){
				write.csv(pData(eset), file)
			} else {
				write.csv(sub_metadata(pData(eset)), file)
			}
		}
	)	
	observeEvent(input$tab2,{
		output$metadatatext <- renderUI({
			if (input$tab2=="metadata") {
				includeMarkdown("www/metadatatext.Rmd")
			} else {
				return()
			}
		})
	})
	observeEvent(input$tab2,{
		if (input$tab2=="metadata") {
			shinyjs::show("property")
			shinyjs::show("full_meta")
		} else {
			shinyjs::hide("property")
			shinyjs::hide("full_meta")
		}
	})
	observeEvent(input$full_meta,{
		if (input$full_meta=="yes") {
			shinyjs::hide("property")
		} else {
			shinyjs::show("property")
		}
	})
	

	######################## Counts table ########################
	counts_table <- reactive({
        load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		if(ncol(exprs(eset))>5) {
			exprs(eset)[,1:5]
		} else {
			exprs(eset)[,1:ncol(exprs(eset))]
		}
    })
	output$counts <- DT::renderDataTable({
		datatable(
			round(counts_table(), 3)
			,colnames = apply(data.matrix(colnames(counts_table())),1, function(x) paste(unlist(strsplit(x,"_"))[1], unlist(strsplit(x,"_"))[2], sep="\n"))
			,selection = 'none'
			,options = list(pageLength = 10)
		)
	})
	output$downloadcounts <- downloadHandler(
		filename = function() { 
			paste(input$geo_acc, '_counts.csv', sep='') 
		},
		content = function(file) {
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
			write.csv(exprs(eset), file)
		}
	)
	observeEvent(input$tab2,{
		output$countstext <- renderUI({
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
			if (input$tab2=="countsTable" && ncol(exprs(eset))>5) {
				helpText(h5(paste0("First 5 of ",ncol(exprs(eset))," columns are printed because of table width limitation")))
			} else {
				return()
			}
		})
	})
	
	######################## QC report ########################
	
	addResourcePath("datasets", "/opt/raid10/genomics/naim/Thesis/geo/datasets")
	output$multiqc <- renderUI({
		tags$iframe(
			seamless="seamless",
			src=paste0("datasets/", input$geo_acc,"/fastqc/multiqc_report.html"), height=1200, width='100%'
		)
	})

	output$downloadqc <- downloadHandler(
		filename = function() { 
			paste(input$geo_acc, '_QCreport.html', sep='') 
		},
		content = function(file) {
			file.copy(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/fastqc/multiqc_report.html"), file,overwrite = TRUE)
			file.remove("multiqc_report.html")
		},contentType = "text/html"
	)
		
	
	############################### Precomputed plots ###############################
	
	############ Sample Correlation ############
	output$corplot <- renderPlotly({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		corr= cor(exprs(eset),method = "spearman")
		
		m <- list(
		  l = 200,
		  r = 10,
		  b = 100,
		  t = 30,
		  pad = 2
		)
		showticklabels <- ifelse (nrow(corr)>50, FALSE, TRUE)
		plot_ly(x = rownames(corr), y = colnames(corr), z = corr, 
				key = corr, type = "heatmap", source = "heatplot", colors = c("blue","black","yellow")) %>%
		  layout(xaxis= list(showticklabels = showticklabels), 
				 yaxis = list(showticklabels = showticklabels),
				 margin = m, showlegend=FALSE)
	})	

	############ Density plot ############
	output$densityp <- renderPlotly({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		pseudoCount = log2(exprs(eset) + 1)
		df = reshape2::melt(pseudoCount)
		df2 = data.frame(df, Samples = df[,2])
		
		p <- ggplot(df2, aes(x = value, colour = Samples)) + ylim(c(0, 0.25)) +
		geom_density()+  theme(legend.position = "right",legend.text=element_text(size=6), legend.key.size=unit(0.3,"cm"))+
		xlab("log2(counts + 1)")
	
		#df2 = stack(as.data.frame(pseudoCount))
		#p <- ggplot(df2, aes(x = value)) + 
		#  geom_density(aes(fill = Samples), alpha = 0.5) + xlab(expression(log[2](counts + 1)))
		ggplotly(p)	
	
	})	

	############ Heatmap ############
	# exprseset <- reactive({
		# load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		# rownames(exprs(eset)) <- paste0(fData(eset)[,2], " : ", fData(eset)[,1])
		# return(exprs(eset))
	# })
	# metadata <- function(){
		# load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		# return(sub_metadata(pData(eset)))
	# }
	source("R/complex_heatmap.R", local = TRUE)
	#source("R/heatmap_server.R", local = TRUE)
	
	output$heat_prop <- renderUI({	
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		col_names <- colnames(sub_metadata(pData(eset)))
		checkboxGroupInput("heat_property", h4("Select properties"), choices  = col_names, 
							selected=col_names[(length(col_names)-1):length(col_names)])
	})
	
	ngenes <- function(data_mat) {
		data_mat <- data_mat[complete.cases(data_mat), ]
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
		expr_data<- data.matrix(exps[topGenes,])
		return(nrow(expr_data))
	}
	
	output$cluster_genes <- renderUI({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		numericInput("cluster_ngenes", label=h4('Number of highly variable genes'), 
			value= if(ngenes(exprs(eset))>1000) {
				1000
			} else {
				ngenes(exprs(eset))
			}, 
			min = 3, max = ngenes(exprs(eset))
		)
		
	})

	plotHeatFit <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		rownames(exprs(eset)) <- paste0(fData(eset)[,2], " : ", fData(eset)[,1])
		complex_heatmap (data_mat=exprs(eset), metadata=sub_metadata(pData(eset)), property=input$heat_property,
		input$cluster_ngenes, clustering_method=input$clustering_method, type="Fit in Screen")
	})
	
	output$fitscreen = renderUI({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)
		plotOutput("HeatFit", height="auto")
	})
	output$HeatFit = renderPlot({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)  
		print(plotHeatFit())
		
	}
	,height = 700
		# reactive({
			# if(!is.null(input$heat_property)) {
				# df <- sub_metadata(pData(eset))[,input$heat_property,drop=FALSE]
				# df[is.na(df)] <- "NA"
				# df[] <- lapply( df, factor)
				# colnames(df) <- input$heat_property
				# x=sapply(df, function(x) length(unique(x)))
				# nrow = max(x)/2
				# return((length(input$heat_property)*25) + 650+(nrow*10))
			# } else {
				# return(600)
			# }
		# })
	)
	
	output$downloadFit = downloadHandler(
		filename = function() {
			datasetname <- sub(" :.*$", "", input$geo_acc)
			if (length(input$heat_property)==1) {
				return(paste(datasetname,"_",input$heat_property,"_heatmap_fit.png",sep=""))
			} else {
				return(paste(datasetname,"_","(",paste(input$heat_property, collapse="+"),")","_heatmap_fit.png",sep=""))
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
		complex_heatmap (data_mat=exprs(eset), metadata=sub_metadata(pData(eset)), property=input$heat_property, 
		input$cluster_ngenes, clustering_method=input$clustering_method, type="Scrollable")
	}
	output$scroll = renderPlot({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)
		print(plotHeatScroll())
	},  #width=1600, #height=12000
		height = reactive({
			if(input$cluster_ngenes>=1000) {
				value <- max((length(input$heat_property)*25)+(4*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24)) #1 inch=96 pixel
			} else {
				value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24))
			}
			return(value)				
		})
	)
	
	output$downloadScroll = downloadHandler(
		filename = function() {
			datasetname <- sub(" :.*$", "", input$geo_acc)
			if (length(input$heat_property)==1) {
				return(paste(datasetname,"_",input$heat_property,"_heatmap_scroll.png",sep=""))
			} else {
				return(paste(datasetname,"_","(",paste(input$heat_property, collapse="+"),")","_heatmap_scroll.png",sep=""))
			}
		},
		content = function(file) {
			height = reactive({
				if(input$cluster_ngenes>=1000) {
					value <- max((length(input$heat_property)*25)+(4*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24)) #1 inch=96 pixel
				} else {
					value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24))
				}
				return(value)
			})
			png(file, width = 3000, height = height(), res=175, units = "px")
			plotHeatScroll()
			dev.off()
		}	
	)
	
	observeEvent(input$tab2,{
		if (input$tab2=="plots" && input$preplots=="heatmap") {
			shinyjs::show("clustering_method")
		} else {
			shinyjs::hide("clustering_method")
		}
	})
	observeEvent(input$preplots,{
		if (input$tab2=="plots" && input$preplots=="heatmap") {
			shinyjs::show("clustering_method")
		} else {
			shinyjs::hide("clustering_method")
		}
	})
	
	observeEvent(input$tab2,{
		if (input$tab2=="plots" && input$preplots=="heatmap") {
			shinyjs::show("cluster_ngenes")
		} else {
			shinyjs::hide("cluster_ngenes")
		}
	})
	
	observeEvent(input$tab2,{
		output$heatmaptext <- renderUI({
			if (input$tab2=="plots" && input$preplots=="heatmap") {
				includeMarkdown("www/heatmaptext.Rmd")
			} else {
				return()
			}
		})
	})
	
	
	############ PCA ############
	source("R/pca.R", local = TRUE)
	#source("./pca.R")
	
	output$pca_prop <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="pca","")
		)
		selectInput(inputId="pca_property",label=h4('Select a property'),
			choices=colnames(metadata()), selected=colnames(metadata())[length(colnames(metadata()))]
		)
	})
	output$legend = renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="pca","")
		)
		plotOutput("legendplot",  height="1600px")
	
	})
	output$legendplot = renderPlot({
		load("data/allColors.rda")
		groups <- metadata()[,which(colnames(metadata())==input$pca_property)]		
		groups[is.na(groups)] <- "NA"
		colcols <- allColors[seq_along(unique(groups))]

		if (length(unique(groups))>=50) {
			ncolm=2
		} else {
			ncolm=1
		} 

		plot(1, type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i")
		legend("top",title = paste(input$pca_property),bty = "n",cex = 1,ncol=ncolm,
		fill=colcols,
		legend= 
		sapply(unique(groups), function(x) levels(factor(as.character(substring(x,1,25)),exclude=NULL))))
		#dev.off()
	},	height = function(){
			groups <- metadata()[,which(colnames(metadata())==input$pca_property)]
			if (length(unique(groups))>=60) {
				return(1600)
			} 
			if (length(unique(groups))<60 & length(unique(groups))>=40) {
				return(1000)
			}
			if (length(unique(groups))<40 & length(unique(groups))>=15) {
				return(700)
			}
			if (length(unique(groups))<15) {
				return(400)
			}
		}
	)
	plot2D <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		pca(data_mat=exprs(eset), metadata=sub_metadata(pData(eset)), property=input$pca_property, type="2D")

    })
	output$pca2Dplot = renderUI({
		pairsD3Output("p2D",width = "auto",height=800)
	})
	output$p2D <- renderPairsD3({
		print(plot2D())
	})
	output$pca3Dplot <- renderRglwidget({
		try(rgl.close())
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		pca(data_mat=exprs(eset), metadata=sub_metadata(pData(eset)), property=input$pca_property, type="3D")
		rglwidget()
	})
	output$downloadPCA2D = downloadHandler(
		filename = function() {
			datasetname <- input$geo_acc
			return(paste(datasetname,"_",input$pca_property,"_pca_2D.html",sep=""))
		},
		content = function(file) {
			savePairs(plot2D(), file)
		}	
	)
	output$downloadPCA3D = downloadHandler(
		filename = function() {
			datasetname <- input$geo_acc
			return(paste(datasetname,"_",input$pca_property,"_pca_3D.html",sep=""))
		},
		content = function(file) {
			writeWebGL(dir = "webGL",filename=file,snapshot = FALSE, commonParts = TRUE)
		}	
	)
	observeEvent(input$preplots,{
		if (input$preplots=="pca") {
			shinyjs::show("pca_property")
		} else {
			shinyjs::hide("pca_property")
		}
	})

	################################################ Analysis nav page ################################################
	
	source("R/analysis.R", local = TRUE)
	observe({
		hide(selector = "#explorgeo li a[data-value=analyze]")
	})
	observeEvent(input$analysisbtn,{
		if (!is.null(input$analysisbtn)) {
			shinyjs::show(selector = "#explorgeo li a[data-value=analyze]")
		} else {
			shinyjs::hide(selector = "#explorgeo li a[data-value=analyze]")
		}
	})
	
	observeEvent(input$analysisbtn, {
		if (input$analysisbtn) {
			updateTabsetPanel(session, "explorgeo", selected = "analyze")
			updateTextInput(session, 'geo_acc2', value = input$geo_acc)
		}
    })
	
	observeEvent(input$geo_acc2, {
		output$metadata_analysis <- renderUI ({
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc2,"/eset.RData"))
			validate(
				need(!is.null(dim(sub_metadata(pData(eset)))), "Either the variables do not contain meaningful information or they contain single value. Please select another dataset to analyze."),
				need(input$signature=="sig_data","")
			)
			metadata <- sub_metadata(pData(eset))
			col_names <- colnames(metadata)
			selectInput(inputId="analysis_property",label=h4('Select a property'),
				choices=c('',col_names), selected=NULL
			)
		})
	})
	
	metadata_small <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc2,"/eset.RData"))
		metadata <- sub_metadata(pData(eset))
		metadata <- metadata[, input$analysis_property,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- input$analysis_property
		return(metadata)
	})
	counts <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc2,"/eset.RData"))
		return(exprs(eset))
	})

	output$group1 <- renderUI ({
		validate(
			need(input$analysis_property!='', ""),
			need(input$signature=="sig_data","")
		)		
		selectInput(inputId="input_group1", label='Group 1',
		choices=c('',levels(metadata_small()[,1])), selected = NULL)
	})
	output$group2 <- renderUI ({
		validate(
			need(any(input$analysis_property!=''),""),
			need(input$signature=="sig_data","")
		)		
		available <- levels(metadata_small()[,1])[which(levels(metadata_small()[,1])!=input$input_group1)]
		selectInput( inputId="input_group2", label='Group 2',
		choices=c('',available), selected = NULL)
	})
	
	
	sigdata <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc2,"/eset.RData"))
		analysis_genes <- cbind(ID_geneid=fData(eset)[,1], Name_GeneSymbol=fData(eset)[,2])
		metadata <- sub_metadata(pData(eset))
		signature_data= analysis(counts=counts(), metadata=metadata, genes=analysis_genes, property=input$analysis_property,
			covariate=NULL, group1=input$input_group1, group2=input$input_group2)
		return(signature_data)
	})
	
	output$sig_data <- DT::renderDataTable({
		validate(
			need(input$signature=="sig_data",""),
			need(input$analysis_property!='', ""),
			need(input$input_group1!='', ""),
			need(input$input_group2!='', "")
		)
			datatable(
				sigdata()$signaturesData_final
				,selection = 'none'
				,rownames = FALSE
				, filter = 'top'
				,options = list(pageLength = 10)
			)

	})	
	
	output$downloadsigdata <- downloadHandler(
		filename = function() { 
			paste(input$geo_acc, '_signatureData.csv', sep='') 
		},
		content = function(file) {
			write.csv(sigdata()$signaturesData_final, file, row.names=FALSE )
		}
	)	

	observeEvent(input$signature,{
		output$sig_method <- renderUI({
			if (input$signature=='sig_data') {
				includeMarkdown("www/signature_method.Rmd")
			} else {
				return()
			}
		})
	})
	observeEvent(input$analyze,{
		if (!is.null(input$input_group1)) {
			shinyjs::show("sig_data")
		} else {
			shinyjs::hide("sig_data")
		}
	})
	
	observe({
		if (!is.null(input$input_group1)) {
			shinyjs::show("downloadsigdata")
		} else {
			shinyjs::hide("downloadsigdata")
		}
	})
	
	
	############ Power curve ############
	
	pow_dat <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc2,"/eset.RData"))
		y <- DGEList(counts=exprs(eset))
		o <- order(rowSums(y$counts))
		y <- y[o,]
		keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
		y <- y[keep,]
		y$samples$lib.size <- colSums(y$counts)
		y <- calcNormFactors(y) 
		y <- estimateCommonDisp(y)
		y$depth <- min(rowMeans(y$counts))
		y$cv <- sqrt(y$common.dispersion)
		return(y)
	})
		
	output$alpha <- renderUI({
		validate(
			need(any(input$signature=="power"),"")
		)
		sliderInput("sig_lev", label=h4('Level of significance (alpha)') ,min = 0, max = 0.2, value = 0.05, step= 0.01)
	})
	
	output$effect <- renderUI({
		validate(
			need(any(input$signature=="power"),"")
		)
		sliderInput("fc", label=h4('Effect size (fold change)') ,min = 0, max = 5, value = 2, step= 0.1)
	})
	
	output$samples <- renderUI({
		validate(
			need(any(input$signature=="power"),"")
		)
		val = round(rnapower(power=0.999, depth=pow_dat()$depth, cv= pow_dat()$cv, alpha=input$sig_lev, effect=input$fc))
		numericInput("n_samples", label=h4('No. of samples') ,min = 2, max= val+5, step= 2, value = val)
	})
		
	output$powerplot <- renderPlotly({	

		power <- rnapower(n=c(1:input$n_samples), depth=pow_dat()$depth, cv= pow_dat()$cv, alpha=input$sig_lev, effect=input$fc)
		df.power <- as.data.frame(cbind(samples=1:input$n_samples, power))
		if(input$n_samples>10) {
			p.line <- ggplot(df.power) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ xlim(1,input$n_samples)	+labs(x="No. of samples", y="Power") #+scale_x_discrete(limit = c(1:input$n_samples))
		} else {
			p.line <- ggplot(df.power) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ labs(x="No. of samples", y="Power") +scale_x_discrete(limit = c(1:input$n_samples))		
		}
		ggplotly(p.line)
	})	
	

})

#sudo nano /var/log/shiny-server/GRSN-shiny-20170705-134227-52677.log