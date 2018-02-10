
# load("data/geo_study.RData")
#load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))

shinyServer(function(input, output, session) {

	load("data/geo_QuerynBrowser_merged.RData")
	load("data/geo_QuerynBrowser_merged_jan28_2018.RData")
	load("data/datatable_to_use.RData")
	load("data/allColors.rda")
	source("R/sub_metadata.R", local = TRUE)
	source("R/dropdownButton.R", local = TRUE)
	source("R/complex_heatmap.R", local = TRUE)
	source("R/interactive_heatmap.R", local = TRUE)
	source("R/pca.R", local = TRUE)
	source("R/analysis.R", local = TRUE)
	source("R/complex_heatmap_sig.R", local = TRUE)

	datafrm <- data.frame(datatable_to_use[,-6] ,check.names = FALSE)
	#datafrm <- datafrm[-which(datafrm[,2]==1),]
	datafr2 <- datafrm[order(datafrm$Species),]
	datafr2 <- data.frame('GEO accession'= paste0('<button id="blah" type="button" class="btn btn-primary action-button" value=', datafr2[,1], '>', datafr2[,1], '</button>'),
					datafr2[,-1], stringsAsFactors=F, check.names=F)
	datafr <- reactiveValues()
	
	observe({
		datafr$x <- if(input$user_geo==""){ 
					datafr2
				} else if (input$user_geo!=""){
					info2 = unlist(lapply(datafr2[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
					datafr2[grep(toupper(input$user_geo),info2,ignore.case=TRUE),,drop=F]
					#datafr2[grep(toupper(input$user_geo),datafr2[,1],ignore.case=TRUE),,drop=F]
				} else {
					return()
				}	
	})

	filenames <- function(){
		details <- file.info(list.files(path="/opt/raid10/genomics/naim/Thesis/geo/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
		files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
		x=gsub("^.*/", "", sub(".txt","",files))
		dir_name <- list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request", full.names = F, recursive = FALSE)
		return(c(dir_name,x))
	}
	observeEvent(input$user_geo, {
		if(!is.null(input$start_process)){
			updateSelectInput(session, 'geo_next_lp', label="User requested datasets in the processing queue", choices=filenames(), selected=filenames()[1])
		}
	})
	
	output$study_stat <- renderPlotly({
		x <- geo_QuerynBrowser_merged[,c(1,4)]
		x$study_status <- ifelse(geo_QuerynBrowser_merged[,1] %in% datatable_to_use[,1], "Processed", "In progress")
		
		status <- c('Processed','In progress')
		#processed  <- x[which(x[,3]=='Processed'),c(2,3)]
		processed  <- datatable_to_use[,3,drop=F]
		waiting  <- x[which(x[,3]=='In progress'),c(2,3)]
		hs <- c(nrow(processed[which(processed[,1]=='Homo sapiens'),,drop=F]), nrow(waiting[which(waiting[,1]=='Homo sapiens'),,drop=F])) 
		mm <- c(nrow(processed[which(processed[,1]=='Mus musculus'),,drop=F]), nrow(waiting[which(waiting[,1]=='Mus musculus'),,drop=F])) 
		rn <- c(nrow(processed[which(processed[,1]=='Rattus norvegicus'),,drop=F]), nrow(waiting[which(waiting[,1]=='Rattus norvegicus'),,drop=F])) 
		dat_study <- data.frame("Status"=as.factor(status), "Homo_sapiens"=hs,'Mus_musculus'=mm, 'Rattus_norvegicus'=rn, check.names=F )

		plot_ly(dat_study, x = ~Status, y = ~Homo_sapiens, type = 'bar', name = 'Homo sapiens') %>%
			add_trace(y = ~Mus_musculus, name = 'Mus musculus') %>%
			add_trace(y = ~Rattus_norvegicus, name = 'Rattus norvegicus') %>%
		layout(title= 'Processing status of GEO RNA-seq datasets', yaxis = list(title = 'Count'), barmode = 'stack')
	})

	output$sample_stat <- renderPlotly({
	
		x=datatable_to_use[,c(2,3)]
		species <- c('Homo sapiens','Mus musculus', 'Rattus norvegicus') 
		samples <- c(sum(x[which(x[,2]=='Homo sapiens'),1]), sum(x[which(x[,2]=='Mus musculus'),1]), sum(x[which(x[,2]=='Rattus norvegicus'),1]))
		dat_sam <- data.frame("Species"=as.factor(species), "Samples"=as.numeric(samples))
		
		plot_ly(dat_sam, labels = ~Species, values = ~samples, type = 'pie',
			textposition = 'inside',
			textinfo = 'label+text',
			insidetextfont = list(color = '#FFFFFF'),
			hoverinfo = 'percent',
			text = ~paste(samples),
			showlegend = FALSE) %>%
		layout(title = 'Number of processed GEO samples',
			xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
			yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
	})
	
	output$datatable <- DT::renderDataTable ({
		datatable(datafr$x,
		selection = 'none',
		escape = FALSE,
		rownames = FALSE, filter = 'top', 
		options = list(searchHighlight = TRUE,
			columnDefs = list(
				list(width = '45%', targets =4),
				list(className = 'dt-center', targets =0:2
				)
			),
			pageLength = 10
		)) %>%formatStyle(1, valueColumns=1,color = '#337ab7', cursor = 'pointer')		
    })
	
	observe({
		hide(selector = "#grin li a[data-value=explore]")
	})
	observeEvent(input$datatable_cell_clicked,{
		info = input$datatable_cell_clicked
		if(!is.numeric(info$value)){
			info2 = unlist(lapply(info$value, function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		} else {
			info2=NULL
		}
		#if (!is.null(info$value) && info$col == 0) {
		if (!is.null(info2) && info$col == 0) {
			shinyjs::show(selector = "#grin li a[data-value=explore]")
		} 
	})
		
	observeEvent(input$datatable_cell_clicked, {
		info = input$datatable_cell_clicked
		if(!is.numeric(info$value)){
			info2 = unlist(lapply(info$value, function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		} else {
			info2=NULL
		}
		#if (is.null(info$value) || info$col != 0) return()
		if (is.null(info2) || info$col != 0) return()
		updateTabsetPanel(session, "grin", selected = "explore")
		#updateTextInput(session, 'geo_acc', value= info$value)
		updateTextInput(session, 'geo_acc', value= info2)
	})
	
	
	############### organize eset as master data for all purpose ##############
	
	collapse_data <- reactive({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",toupper(input$geo_acc),"/eset.RData"))
		id <- unlist(lapply(colnames(exprs(eset)), function(x) unlist(strsplit(x, "_"))[2]))
		pdata <- pData(eset)
		fdata <- fData(eset)
		if("geo_accession" %in% names(pdata)) {
			if(length(unique(id))!= nrow(pdata)){
				pdata$ids <- id
				colData <- data.frame(subject=pdata$ids) 
				dds <-DESeqDataSetFromMatrix(round(exprs(eset)),colData,design=~subject)
				ddsCollapsed <- collapseReplicates(dds,groupby = dds$subject,renameCols = F)
				countsTable <- DESeq2::counts(ddsCollapsed)
				pdata <- pdata[colnames(countsTable), -which(colnames(pdata)=="ids")]
				colnames(countsTable) <- as.character(colData(ddsCollapsed)$subject)
				
				genes <- fdata[complete.cases(fdata), 2, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% rownames(genes)),]
				rownames(countsTable) <- paste0(rownames(countsTable), " : ", genes[,1])
				#rownames(countsTable) <- paste0(fData(eset)[,1], " : ", fData(eset)[,2])
				rownames(pdata) <- colnames(countsTable)
				sub_pdata <- sub_metadata(pdata)
				collapse_data <- list(counts=countsTable, pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
				return(collapse_data)
			} else {
				countsTable <- round(exprs(eset))
				colnames(countsTable) <- id
				genes <- fdata[complete.cases(fdata), 2, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% rownames(genes)),]
				rownames(countsTable) <- paste0(rownames(countsTable), " : ", genes[,1])
				rownames(pdata) <- id
				sub_pdata <- sub_metadata(pdata)
				collapse_data <- list(counts=countsTable, pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
				return(collapse_data)
			}
		} else {
			countsTable <- round(exprs(eset))
			colnames(countsTable) <- id
			genes <- fdata[complete.cases(fdata), 2, drop=F]
			countsTable <- countsTable[which(rownames(countsTable) %in% rownames(genes)),]
			rownames(countsTable) <- paste0(rownames(countsTable), " : ", genes[,1])
			rownames(pdata) <- colnames(countsTable)
			sub_pdata <- sub_metadata(pdata)
			dat <- list(counts=countsTable,pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
			return(dat)
		}
	})

	######################## Description ########################
	
	geo_rows <- reactive ({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		return(length(unique(unlist(lapply(colnames(exprs(eset)), function(x) unlist(strsplit(x,"_"))[2])))))
	})
	srr_rows <- reactive ({
		load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
		return(ncol(exprs(eset)))
	})

	observeEvent(input$geo_acc,{
		output$geo_summary <- renderText(
			paste(
				paste0(h4("Study link: "), 
				paste0('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
					toupper(input$geo_acc), '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'geoinfo\', 1)" target="_blank">', toupper(input$geo_acc), '</a>')),
				h4(paste0("No. of GEO samples : ",geo_rows())),
				h4(paste0("No. of SRA runs : ",srr_rows())),
				paste(h4("Title: "), geo_QuerynBrowser_merged_jan28_2018[which(geo_QuerynBrowser_merged_jan28_2018$gse_id==input$geo_acc),5], sep="\n"),
				paste(h4("Summary: "),geo_QuerynBrowser_merged_jan28_2018[which(geo_QuerynBrowser_merged_jan28_2018$gse_id==input$geo_acc),6], sep="\n"),
			sep="\n")
		)
	})
  
	
	######################## metadata column ########################
	
	output$metadata_column <- renderUI({	
		validate(
			need(any(input$tab2=="metadata"),"")
		)
		#metadata <- sub_metadata(collapse_data()$pdata)
		metadata <- collapse_data()$sub_pdata
		col_names <- names(metadata)
		checkboxGroupInput("property", h4("Select properties"), choices  = col_names, selected=col_names[1:length(col_names)])
	})
	
	output$full_metacol <- renderUI({
		validate(
			need(any(input$tab2=="metadata"),"")
		)
		radioButtons("full_meta", label = h4("Show full metadata"), c("yes","no"), selected="no")
	})
		
	output$metadata_full <- DT::renderDataTable({
		validate(
			need(any(input$tab2=="metadata"),"")
		)
		datatable(collapse_data()$pdata, selection = 'none',options = list(searchHighlight = TRUE, pageLength = if(nrow(collapse_data()$pdata)>=10) {10} else {nrow(collapse_data()$pdata)},
					scrollX = TRUE), rownames= TRUE
		)
	})
	output$metadata <- DT::renderDataTable({
		validate(
			need(any(input$tab2=="metadata"),"")
		)
		datatable(collapse_data()$sub_pdata[, colnames(collapse_data()$sub_pdata) %in% input$property, drop = FALSE],
		selection = 'none',
		options = list(searchHighlight = TRUE, pageLength = if(nrow(collapse_data()$sub_pdata)>=10) {10} else {nrow(collapse_data()$sub_pdata)}, scrollX = TRUE),
			rownames= TRUE
		)
	})
		
	output$downloadmeta <- downloadHandler(
		filename = function() { 
			if(input$full_meta=="yes"){
				paste(toupper(input$geo_acc), '_full_metadata.csv', sep='') 
			} else {
				paste(toupper(input$geo_acc), '_filtered_metadata.csv', sep='') 
			}
		},
		content = function(file) {
			load(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc,"/eset.RData"))
			if(input$full_meta=="yes"){
				write.csv(collapse_data()$pdata, file)
			} else {
				write.csv(collapse_data()$sub_pdata, file)
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
	observeEvent(input$tab2,{
		if (input$tab2=="metadata") {
			shinyjs::show("downloadmeta")
		} else {
			shinyjs::hide("downloadmeta")
		}
	})
	observeEvent(input$full_meta,{
		if (input$full_meta=="yes") {
			shinyjs::show("metadata_full")
			shinyjs::hide("metadata")
		} else {
			shinyjs::hide("metadata_full")
			shinyjs::show("metadata")
		}
	})

	######################## Counts table ########################
	
	output$counts <- DT::renderDataTable({	
		ct <- data.frame(gene_symbol= unlist(lapply(strsplit(rownames(collapse_data()$counts), ' : ', fixed = TRUE), '[', 2)), collapse_data()$counts, stringsAsFactors=F)
		rownames(ct) <- unlist(lapply(strsplit(rownames(collapse_data()$counts), ' : ', fixed = TRUE), '[', 1))
		datatable(
			ct,
			#collapse_data()$counts, 
			selection = 'none',options = list(pageLength = if(nrow(collapse_data()$counts)>=10) {10} else {nrow(collapse_data()$counts)}, scrollX = TRUE)
		)
	})
	output$downloadcounts <- downloadHandler(
		filename = function() { 
			if(input$download_type_counts == "Gene level") {
				paste(toupper(input$geo_acc), '_GeneLevel_counts.csv', sep='') 
			} else {
				paste(toupper(input$geo_acc), '_TranscriptLevel_counts.csv', sep='') 
			}
		},
		content = function(file) {
			if(input$download_type_counts == "Gene level") {
				ct <- data.frame(gene_symbol= unlist(lapply(strsplit(rownames(collapse_data()$counts), ' : ', fixed = TRUE), '[', 2)), collapse_data()$counts, stringsAsFactors=F)
				rownames(ct) <- unlist(lapply(strsplit(rownames(collapse_data()$counts), ' : ', fixed = TRUE), '[', 1))
				write.csv(ct, file)
			} else {
				setwd(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",input$geo_acc))
				tr <- read.csv(file="transcripts_counts.csv", header=TRUE, row.names=1)
				write.csv(tr, file)
			}
		}
	)
	observeEvent(input$tab2,{
		if (input$tab2=="countsTable") {
			shinyjs::show("downloadcounts")
			shinyjs::show("download_type_counts")
		} else {
			shinyjs::hide("downloadcounts")
			shinyjs::hide("download_type_counts")
		}
	})

	
	######################## QC report ########################
	
	addResourcePath("datasets", "/opt/raid10/genomics/naim/Thesis/geo/datasets")
	output$multiqc <- renderUI({
		tags$iframe(
			seamless="seamless",
			src=paste0("datasets/", toupper(input$geo_acc),"/fastqc/multiqc_report.html"), height=1200, width='100%'
		)
	})
	output$downloadqc <- downloadHandler(
		filename = function() { 
			paste(toupper(input$geo_acc), '_QCreport.html', sep='') 
		},
		content = function(file) {
			file.copy(paste0("/opt/raid10/genomics/naim/Thesis/geo/datasets/",toupper(input$geo_acc),"/fastqc/multiqc_report.html"), file,overwrite = TRUE)
			file.remove("multiqc_report.html")
		},contentType = "text/html"
	)
		
	############################### Visualization ###############################
	
	############ Sample Correlation ############
	
	output$corplot <- renderPlotly({
		corr= cor(collapse_data()$counts,method = "spearman")
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
		pseudoCount = cpm(collapse_data()$counts, normalized.lib.sizes=FALSE, log=TRUE)
		df = reshape2::melt(pseudoCount)
		df2 = data.frame(df, Samples = df[,2])
		
		p <- ggplot(df2, aes(x = value, colour = Samples)) + ylim(c(0, 0.25)) +
		geom_density()+  theme(legend.position = "right",legend.text=element_text(size=6), legend.key.size=unit(0.3,"cm"))+
		xlab("log2(counts per million)")
	
		ggplotly(p)	
	})	

	############ Heatmap ############

	output$heat_prop <- renderUI({	
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		col_names <- colnames(collapse_data()$sub_pdata)
		checkboxGroupInput("heat_property", h4("Select properties"), choices  = col_names, 
		selected= if(length(col_names)>2) {
					col_names[(length(col_names)-1):length(col_names)]
				} else {
					col_names[1:length(col_names)]
				}
		)
	})
	data_dl <- function(data_mat, n_genes) {
		dge <- DGEList(counts=data_mat)
		dge <- calcNormFactors(dge)
		dge_cpm <- cpm(dge, log=TRUE)
		#data_mat <- cpm(data_mat, normalized.lib.sizes=TRUE, log=TRUE)
		exps= as.matrix(dge_cpm) - rowMeans(as.matrix(dge_cpm), na.rm=T)
		medAbsDev <-apply(exps,1,function(x) median(abs(x)))
		topGenes <-order(medAbsDev,decreasing=T)
		topGenes= topGenes[!is.na(topGenes)]

		expr_data<- data.matrix(dge_cpm[topGenes,])
		expr_data_in <- expr_data[1:n_genes,]
		gene_sym <- unlist(lapply(rownames(expr_data_in), function(x) unlist(strsplit(x, ' : '))[2]))
		dat2 <- data.frame(GeneSymbol=gene_sym, expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	}
	
	output$cluster_genes <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		numericInput("cluster_ngenes", label=h4('Number of highly variable genes'), 
			value= if(nrow(collapse_data()$counts)>100) {
				100
			} else {
				nrow(collapse_data()$counts)
			}, 
			min = 1, max = nrow(collapse_data()$counts)
		)
		
	})
	output$cluster_meth <- renderUI({	
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		selectInput("clustering_method", label = h4("Grouping genes and samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected=
			#if(any(apply(collapse_data()$counts,1,sum, na.rm=TRUE)==0)) {
			if(input$cluster_ngenes<10){
				"Euclidean distance"
			} else {
				"Pearson correlation"
			}
		)
	})

	plotHeatFit <- function(){
		if(!is(try(complex_heatmap (data_mat=collapse_data()$counts, metadata=collapse_data()$sub_pdata, property=input$heat_property,
		n_genes=input$cluster_ngenes, clustering_method=input$clustering_method, type="Fit in Screen"), silent=TRUE), 'try-error')) warnings("Please try another")
	}
	
	output$fitscreen = renderUI({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)
		plotOutput("HeatFit", height="auto") %>% withSpinner(type=5)
	})
	output$HeatFit = renderPlot({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)  
		print(plotHeatFit())
		
	}
		,height = reactive({
			if(!is.null(input$heat_property)) {
				df <- collapse_data()$sub_pdata[,input$heat_property,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- input$heat_property
				x=sapply(df, function(x) length(unique(x)))
				nrow = max(x)/2
				if(input$cluster_ngenes<=50)
					return((length(input$heat_property)*25) + (200+ (input$cluster_ngenes*5.5))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				else {
					return((length(input$heat_property)*25) + 550 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				}
			} 
		})
	)
	output$downloadFit <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_fit=="png"){
				if (length(input$heat_property)==1) {
					return(paste(datasetname,"_(",input$heat_property,")_",input$cluster_ngenes,"genes_heatmap_fit.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property, collapse="+"),")","_",input$cluster_ngenes,"genes_heatmap_fit.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes,"Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_fit == "png") {
				png(file, width = 18, height = 16, units = "in",res=200)
				plotHeatFit()
				dev.off()
			} else if(input$download_type_fit == "csv"){
				dat2 <- data_dl(collapse_data()$counts, input$cluster_ngenes)
				write.csv(dat2, file)
			} else {
				return()
			}
		}
	)

	### Scrollable

	plotHeatScroll <- function(){
		if(!is(try(complex_heatmap (data_mat=collapse_data()$counts, metadata=collapse_data()$sub_pdata, property=input$heat_property, 
		input$cluster_ngenes, clustering_method=input$clustering_method, type="Scrollable"), silent=TRUE), 'try-error')) warnings("Please try another")
	}
	output$scroll = renderUI({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)
		plotOutput("HeatScroll", height="auto") %>% withSpinner(type=5)
	})
	
	output$HeatScroll = renderPlot({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)  
		print(plotHeatScroll())
	}
		,height = reactive({
			if(input$cluster_ngenes>=100) {
				value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24))
				#value <- max((length(input$heat_property)*25)+(4*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24)) #1 inch=96 pixel
			} else {
				#value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24))
				value <- (length(input$heat_property)*25) + 400+(input$cluster_ngenes*20)
			}
			return(value)				
		})
	)
	
	output$downloadScroll <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scroll=="png"){
				if (length(input$heat_property)==1) {
					return(paste(datasetname,"_","(",input$heat_property,")","_",input$cluster_ngenes,"genes_heatmap_scroll.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property, collapse="+"),")","_",input$cluster_ngenes,"genes_heatmap_scroll.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes,"Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_scroll == "png") {
				height = reactive({
					if(input$cluster_ngenes>=1000) {
						value <- max((length(input$heat_property)*25)+(4*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				plotHeatScroll()
				dev.off()
			} else{
				dat2 <- data_dl(collapse_data()$counts, input$cluster_ngenes)
				write.csv(dat2, file)
			}
		}
	)

	### Interactive
	
	height_in = reactive({
		if(!is.null(input$heat_property)) {
			df <- collapse_data()$sub_pdata[,input$heat_property,drop=FALSE]
			df[is.na(df)] <- "NA"
			df[] <- lapply( df, factor)
			colnames(df) <- input$heat_property
			x=sapply(df, function(x) length(unique(x)))
			nrow = max(x)/2
			if(input$cluster_ngenes<=50)
				return((length(input$heat_property)*25) + (200+ (input$cluster_ngenes*5.5))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
			else {
				return((length(input$heat_property)*25) + 550 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
			}
		} 
	})
	plotInteractiveHeat <- function(){
		#if(!is(try(interactive_heatmap (data_mat=collapse_data()$counts, metadata=collapse_data()$sub_pdata, property=input$heat_property, 
		#n_genes=input$cluster_ngenes, clustering_method=input$clustering_method), silent=TRUE), 'try-error')) warnings("Please try another")
		interactive_heatmap (data_mat=collapse_data()$counts, metadata=collapse_data()$sub_pdata, property=input$heat_property, 
		n_genes=input$cluster_ngenes, clustering_method=input$clustering_method)
	}

	output$iheatmapui = renderUI({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)
		iheatmaprOutput("iheatmap", width="auto", height = height_in()) %>% withSpinner(type=5)
	})
	output$iheatmap = renderIheatmap({
		validate(
			need(!is.null(input$heat_property), "please select a property")
		)  
		#plotInteractiveHeat()
		interactive_heatmap (data_mat=collapse_data()$counts, metadata=collapse_data()$sub_pdata, property=input$heat_property, 
		n_genes=input$cluster_ngenes, clustering_method=input$clustering_method)

	})

	
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
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="fitscrn"){
			shinyjs::show("downloadFit")
			shinyjs::show("download_type_fit")
		} else {
			shinyjs::hide("downloadFit")
			shinyjs::hide("download_type_fit")
		}
	})
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="scrolb"){
			shinyjs::show("downloadScroll")
			shinyjs::show("download_type_scroll")
		} else {
			shinyjs::hide("downloadScroll")
			shinyjs::hide("download_type_scroll")
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
	
	output$pca_prop <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="pca","")
		)
		selectInput(inputId="pca_property",label=h4('Select a property'),
			choices=colnames(collapse_data()$sub_pdata), selected=colnames(collapse_data()$sub_pdata)[length(colnames(collapse_data()$sub_pdata))]
		)
	})
	output$legend = renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="pca","")
		)
		#tags$div(class = "container_pca",
			plotOutput("legendplot", height='1600px')
		#)
	})
	
	output$legendplot = renderPlot({
		groups <- collapse_data()$sub_pdata[,which(colnames(collapse_data()$sub_pdata)==input$pca_property)]		
		groups[is.na(groups)] <- "NA"
		colcols <- allColors[seq_along(unique(groups))]

		if (length(unique(groups))>=50) {
			ncolm=2
		} else {
			ncolm=1
		} 

		plot(1, type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i")
		legend("top",title = input$pca_property
		,bty = "n",cex = 1, ncol=ncolm,
		fill=colcols, 
		legend= 
		sapply(as.character(unique(groups)), function(x) paste(strsplit(x, "(?<=.{20})", perl = TRUE)[[1]], collapse="\n")))
		#sapply(unique(groups), function(x) levels(factor(as.character(x),exclude=NULL))))
	}	
		,height = function(){
			groups <- collapse_data()$sub_pdata[,which(colnames(collapse_data()$sub_pdata)==input$pca_property)]
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
		data_mat <- cpm(collapse_data()$counts, normalized.lib.sizes=FALSE, log=TRUE)
		pca(data_mat=data_mat, metadata=collapse_data()$sub_pdata, property=input$pca_property, type="2D")
    })
	output$pca2Dplot = renderUI({
		pairsD3Output("p2D",width = "auto",height=800) %>% withSpinner(type=5)
	})
	output$p2D <- renderPairsD3({
		print(plot2D())
	})
	output$pca3Dplot <- renderRglwidget({
		try(rgl.close())
		data_mat <- cpm(collapse_data()$counts, normalized.lib.sizes=FALSE, log=TRUE)
		pca(data_mat=data_mat, metadata=collapse_data()$sub_pdata, property=input$pca_property, type="3D")
		rglwidget()
	})
	output$downloadPCA2D = downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			return(paste(datasetname,"_",input$pca_property,"_pca_2D.html",sep=""))
		},
		content = function(file) {
			savePairs(plot2D(), file)
		}	
	)
	output$downloadPCA3D = downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			return(paste(datasetname,"_",input$pca_property,"_pca_3D.html",sep=""))
		},
		content = function(file) {
			writeWebGL(dir = ".",filename=file,snapshot = FALSE, commonParts = TRUE)
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
	
	observe({
		hide(selector = "#grin li a[data-value=analyze]")
	})
	observeEvent(input$analysisbtn,{
		if (!is.null(input$analysisbtn)) {
			shinyjs::show(selector = "#grin li a[data-value=analyze]")
		} else {
			shinyjs::hide(selector = "#grin li a[data-value=analyze]")
		}
	})
	observeEvent(input$analysisbtn, {
		if (input$analysisbtn) {
			updateTabsetPanel(session, "grin", selected = "analyze")
		}
    })

	# observeEvent(input$geo_acc, {
		# output$geo_acc2_ui <- renderUI({	
			# info = input$datatable_cell_clicked
			# validate(
				# need(any(c(!is.null(info$value),input$analysisbtn)), "")
			# )
			# selectInput('geo_acc2', h3('Selected study'), choices= c('',datafr$x[,1]),selected= input$geo_acc)
		# })
	# })
	# observeEvent(input$geo_acc2, {
		# updateSelectInput(session, 'geo_acc',choices= c('',datafr$x[,1]),selected= input$geo_acc2)
	# })

	observe({
		info = input$datatable_cell_clicked
		if (!is.null(info$value) || input$analysisbtn) {
			updateTextInput(session, 'geo_acc2', value = toupper(input$geo_acc))
		}
    })
	observeEvent(input$geo_acc2, {
		updateTextInput(session, 'geo_acc', value = toupper(input$geo_acc2))
	})

	################## Power Analysis ##################
	
	####### Power curve #######
	
	pow_dat <- reactive({
		group <- as.factor(collapse_data()$sub_pdata[,1])
		if(any(as.numeric(table(group))<2)){
			y <- DGEList(counts=collapse_data()$counts, genes= collapse_data()$fdata)
			o <- order(rowSums(y$counts))
			y <- y[o,]
			keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
			y <- y[keep, , keep.lib.sizes=FALSE]
			y$samples$lib.size <- colSums(y$counts)
			y <- calcNormFactors(y) 
			y <- estimateCommonDisp(y)
			y <- estimateTagwiseDisp(y)
			y$depth <- y$pseudo.lib.size/1000000
			y$cv <- sqrt(y$common.dispersion)
			y$bcov_tagwise <- sqrt(y$tagwise.dispersion)
			return(y)
		} else {
			y <- DGEList(counts=collapse_data()$counts, genes= collapse_data()$fdata, group=group)
			o <- order(rowSums(y$counts))
			y <- y[o,]
			keep <- rowSums(cpm(y)>1) >= min(summary(factor(y$samples$group)))
			y <- y[keep, , keep.lib.sizes=FALSE]
			y$samples$lib.size <- colSums(y$counts)
			y <- calcNormFactors(y) 
			y <- estimateCommonDisp(y)
			y <- estimateTagwiseDisp(y)
			y$depth <- y$pseudo.lib.size/1000000
			y$cv <- sqrt(y$common.dispersion)
			y$bcov_tagwise <- sqrt(y$tagwise.dispersion)
			return(y)
		}
	})
		
	output$alpha <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(any(c(input$sigtab=="power_curve", input$sigtab=="detec_gene")),"")
		)
		numericInput("sig_lev", label=h4('Level of significance (alpha)') ,min = 0, max = 0.1, value = 0.05, step= 0.01)
	})
	
	output$effect <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(any(c(input$sigtab=="power_curve", input$sigtab=="detec_gene")),"")
		)
		sliderInput("fc", label=h4('Effect size (fold change)') ,min = 1, max = 5, value = 2, step= 0.1)
	})
	
	output$samples <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(any(c(input$sigtab=="power_curve", input$sigtab=="detec_gene")),"")
		)
		val = ncol(collapse_data()$counts)
		numericInput("n_samples", label=h4('No. of samples') ,min = 2, max= val+5, step= 2, value = val)
	})

	output$depth_ui <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(input$sigtab=="power_curve","")
		)
		val = round(pow_dat()$depth,2)
		numericInput("depth", label=h4('Sequencing depth (in million)') ,min = 2, max= val+5, step= 0.5, value = val)
	})
	
	df.power <- reactive ({
		#power <- rnapower(n=c(1:input$n_samples), depth=pow_dat()$depth, cv= pow_dat()$cv, alpha=input$sig_lev, effect=input$fc)
		power <- rnapower(n=c(1:input$n_samples), depth=input$depth, cv= pow_dat()$cv, alpha=input$sig_lev, effect=input$fc)
		df.power <- as.data.frame(cbind(samples=1:input$n_samples, power))
		return(df.power)
	})
	output$powerplot <- renderPlotly({	
		if(input$n_samples>10) {
			p.line <- ggplot(df.power()) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ xlim(1,input$n_samples)	+labs(x="No. of samples", y="Power") #+scale_x_discrete(limit = c(1:input$n_samples))
		} else {
			p.line <- ggplot(df.power()) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ labs(x="No. of samples", y="Power") +scale_x_discrete(limit = c(1:input$n_samples))		
		}
		ggplotly(p.line)
	})	
	
	observeEvent(input$signature,{
		output$pow_method <- renderUI({
			#if (input$signature=='power') {
			if(input$signature=='power' & input$sigtab=="power_curve"){
				includeMarkdown("www/power_method.Rmd")
			} else {
				return()
			}
		})
	})

	####### gene detect #######
	
	output$power_ui <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(input$sigtab=="detec_gene","")
		)
		val = 0.8
		numericInput("power_level", label=h4('Power') ,min = 0.1, max= 1, step= 0.1, value = val)
	})

	output$search_gene_ui <- renderUI({
		validate(
			need(input$signature=="power","")
			,need(input$sigtab=="detec_gene","")
		)
		selectInput(inputId="search_gene",label=h4('Search for gene of interest'), 
			choices= toupper(pow_dat()$genes[,1]), multiple = TRUE, selected= NULL
		)
	})
	
	to_vary <- reactive ({
		logCPM_range <- if(min(pow_dat()$AveLogCPM)>0) {
							0:(max(pow_dat()$AveLogCPM)*100) /100
						} else {
							(min(pow_dat()$AveLogCPM)*100):(max(pow_dat()$AveLogCPM)*100) /100
						}
		depth_range <- 2^logCPM_range * pow_dat()$depth
		
		effect_vary <- NULL
		effect_vary <- c(effect_vary,input$fc)
		alpha_vary <- NULL
		alpha_vary <- c(alpha_vary,input$sig_lev)
		samples_vary <- NULL
		samples_vary <- c(samples_vary,input$n_samples)
		power_vary <- NULL
		power_vary <- c(power_vary,input$power_level)
		
		bcov_range <- as.matrix(rnapower(depth=depth_range, n=samples_vary, alpha=alpha_vary, power=power_vary, effect=effect_vary))
		colnames(bcov_range) <- "line of detectability"
		
		df_logCPM_bcov_rnapower <- melt(data.frame(logCPM=logCPM_range, bcov=bcov_range), id="logCPM")
		df_logCPM_bcov <- data.frame(logCPM=pow_dat()$AveLogCPM, bcov=pow_dat()$bcov_tagwise, genes=pow_dat()$genes[,1])
		df_logCPM_bcov_rnapower_overlap <- data.frame(logCPM=df_logCPM_bcov$logCPM, variable="genes", value=df_logCPM_bcov$bcov, genes=df_logCPM_bcov$genes)
		df_logCPM_bcov_rnapower$variable <- gsub(c("bcov.line.of.detectability|line.of.detectability"), paste0("line of detectability"), df_logCPM_bcov_rnapower$variable)
		df_logCPM_bcov_rnapower$genes <- paste0("NA", 1:nrow(df_logCPM_bcov_rnapower))
		df_logCPM_bcov_rnapower_overlap2 <- rbind(df_logCPM_bcov_rnapower_overlap, df_logCPM_bcov_rnapower)
		return(df_logCPM_bcov_rnapower_overlap2)
	})
			
	output$detect_gene <- renderPlotly({
		dataset <- to_vary()
		p <- ggplot(dataset, aes(x=logCPM, y=value, colour=variable)) + 
			geom_point(size=1, alpha=0.6, aes(text=paste0("gene symbol: ",to_vary()$genes))) +
			labs(x="Average log(counts per million)", y="Biological coefficient of variation") + 
			guides(colour=guide_legend(override.aes = list(size=10, alpha=0.7), title.position="top", title="")) + 
			theme(legend.text=element_text(size = 14), axis.text=element_text(size=10), axis.title=element_text(size=12))
		gene_det <- p + geom_point(data=to_vary()[which(toupper(to_vary()$genes) %in% input$search_gene),],color="blue",size=1.5)
		
		ggplotly(gene_det, source="gene_det_source", tooltip=c("text","x","y"))
	})
	
	observeEvent(input$signature,{
		output$gene_detect_ui <- renderUI({
			if(input$signature=='power' & input$sigtab=="detec_gene"){
				includeMarkdown("www/gene_detetctability.Rmd")
			} else {
				return()
			}
		})
	})

	################## Signature ##################
		
	observeEvent(input$geo_acc2, {
		output$metadata_analysis <- renderUI ({
			validate(
				need(!is.null(dim(collapse_data()$sub_pdata)), "Either the variables do not contain meaningful information or they contain single value. Please select another dataset to analyze."),
				need(input$signature=="sig_data","")
			)
			metadata <- collapse_data()$sub_pdata
			col_names <- colnames(metadata)
			selectInput(inputId="analysis_property",label=h4('Variable of interest'),
				choices=c('',col_names), selected=NULL
			)
		})
	})
	
	metadata_small <- reactive({
		metadata <- collapse_data()$sub_pdata
		metadata <- metadata[, input$analysis_property,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- input$analysis_property
		return(metadata)
	})
	
	output$ana_type <- renderUI({
		validate(
			need(!is.null(dim(collapse_data()$sub_pdata)), "Either the variables do not contain meaningful information or they contain single value. Please select another dataset to analyze."),
			need(input$signature=="sig_data",""),
			need(input$analysis_property!='',"")
		)
		selectInput(inputId="analysis_type",label=h4('Type of comparison'),
			choices= 
				if(length(levels(metadata_small()[,1]))==2 && ncol(collapse_data()$sub_pdata)>1){
					c('Two group without covariate','Two group with covariate')
				} else if(length(levels(metadata_small()[,1]))==2 && ncol(collapse_data()$sub_pdata)==1 && any(as.numeric(table(metadata_small()[,1]))<2)){
					'Two group without covariate'
				} else if(length(levels(metadata_small()[,1]))>2 && ncol(collapse_data()$sub_pdata)==1 && all(as.numeric(table(metadata_small()[,1]))>=2)){
					c('Two group without covariate','Multi group without covariate')
				} else if(length(levels(metadata_small()[,1]))>2 && ncol(collapse_data()$sub_pdata)>1 && all(as.numeric(table(metadata_small()[,1]))>=2)){
					#c('Two group without covariate','Two group with covariate', 'Multi group without covariate','Multi group with covariate')
					c('Two group without covariate','Two group with covariate', 'Multi group without covariate')
				} else {
					'Two group without covariate'
				},
			selected= 'Two group without covariate'
		)
	})
	output$group1 <- renderUI ({
		validate(
			need(input$signature=="sig_data",""),
			need(any(c(input$analysis_type=="Two group without covariate",input$analysis_type=="Two group with covariate")),""),
			need(input$analysis_property!='',"")
		)		
		selectInput(inputId="input_group1", label=h4('Experimental group'),
		choices=c('',levels(metadata_small()[,1])), selected = NULL)
	})
	output$group2 <- renderUI ({
		validate(
			need(input$signature=="sig_data",""),
			need(any(c(input$analysis_type=="Two group without covariate",input$analysis_type=="Two group with covariate")),""),
			need(input$analysis_property!='',"")
		)		
		available <- levels(metadata_small()[,1])[which(levels(metadata_small()[,1])!=input$input_group1)]
		selectInput( inputId="input_group2", label=h4('Control group'),
		choices=c('',available), selected = NULL)
	})	
	output$cov <- renderUI({
		validate(
			need(input$signature=="sig_data",""),
			need(input$analysis_property!='',""),
			need(any(c(input$analysis_type=="Multi group with covariate",input$analysis_type=="Two group with covariate")),"")
		)
		selectInput(inputId="covar",label=h4('Select covariates'), 
			choices= colnames(collapse_data()$sub_pdata)[which(colnames(collapse_data()$sub_pdata)!=input$analysis_property)]
			, multiple = TRUE, selected= NULL
		)
	})
	
	design_data <- reactive({
		metad <- collapse_data()$sub_pdata
		analysis_metadata <- metad[which(metad[,input$analysis_property]==input$input_group1 | metad[,input$analysis_property]==input$input_group2),c(input$analysis_property,input$covar),drop=FALSE]
		group <- analysis_metadata[,input$analysis_property, drop=F]
		covar <- analysis_metadata[,input$covar, drop=F]
		design_data <- cbind(covar, group)
		design_data[] <- lapply(design_data, factor)
		return(design_data)
	})
	design_rank <- reactive({
		design <- model.matrix(~., data=design_data())
		ne <- nonEstimable(design)	
		return(ne)
	})
	
	output$rank_warn <- renderText({
		validate(
			need(any(input$analysis_type=='Two group with covariate',input$analysis_type=='Multi group with covariate'),"")
			,need(any(!is.null(design_rank()), any(as.numeric(sapply(design_data()[,sapply(design_data(), is.factor)], nlevels))<2)), "")
		)
		if(any(as.numeric(sapply(design_data()[,sapply(design_data(), is.factor)], nlevels))<2)){
			paste("All the variables need to have 2 or more levels")
		} else {
			paste("Design matrix is not of full rank.  The following coefficients are not estimable:\n", paste(design_rank(), collapse = " "))
		}
	})
	
	sigdata <- reactive({
		filepath <- paste('/mnt/raid/tmp/iLincs/signatures/',toupper(input$geo_acc2), '_signatureData.txt', sep='')
		analysis_genes <- cbind(Ensembl_ID=rownames(collapse_data()$fdata), Gene_symbol=collapse_data()$fdata[,1])
		metadata <- collapse_data()$sub_pdata
		if(input$analysis_type=='Two group without covariate' && input$input_group2!=''){
			signature_data= analysis_2g_wocov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property=input$analysis_property,
				group1=input$input_group1, group2=input$input_group2)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,3], Significance_pvalue=signature_data$signaturesData_final[,4], stringsAsFactors=F)
			#write.table(signature_data$signaturesData_final[,-1],filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$analysis_type=='Two group with covariate' && input$covar!=''){
			signature_data= analysis_2g_withcov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property=input$analysis_property,
				covariate=input$covar, group1=input$input_group1, group2=input$input_group2)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,3], Significance_pvalue=signature_data$signaturesData_final[,4], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$analysis_type=='Multi group without covariate') {
			signature_data= analysis_mg_wocov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property=input$analysis_property)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,3], Significance_pvalue=signature_data$signaturesData_final[,4], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		# } else if(input$analysis_type=='Multi group with covariate' && input$covar!='') {
			# signature_data= analysis_mg_withcov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property=input$analysis_property,
				# covariate=input$covar)
			# write.table(signature_data$signaturesData_final[,-1],filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else {
			return()
		}
		return(signature_data)
	})
	
	output$sig_data <- DT::renderDataTable({
		datatable(
			sigdata()$signaturesData_final
			,selection = 'none'
			,rownames = FALSE
			, filter = 'top' 
			,options = list(pageLength = 10, scrollX = TRUE
				,columnDefs = list(list(className = 'dt-center', targets ="_all"))
			)
		)
	})	
	output$downloadsigdata <- downloadHandler(
		filename = function() { 
			paste(toupper(input$geo_acc2), '_signatureData.csv', sep='') 
		},
		content = function(file) {
			write.csv(sigdata()$signaturesData_final[,-1], file, row.names=FALSE )
		}
	)	
	observeEvent(input$signature,{
		output$sig_method <- renderUI({
			if (input$signature=='sig_data') {
				includeMarkdown("www/signature_method.Rmd")
			}
		})
	})
	observeEvent(input$analysis_type,{
		if (input$analysis_type=='Multi group without covariate') {
			shinyjs::show("downloadsigdata")
		}
	})
	observeEvent(input$analysis_type,{
		if (input$analysis_type=='Two group with covariate') {
			shinyjs::hide("downloadsigdata")
		}
	})
	observeEvent(input$input_group2,{
		if (input$input_group2=='' && input$analysis_type=='Two group without covariate') {
			shinyjs::hide("downloadsigdata")
		} else if (input$input_group2!='' && input$analysis_type!='Two group without covariate'){
			shinyjs::hide("downloadsigdata")
		} else {
			shinyjs::show("downloadsigdata")
		}
	})
	observeEvent(input$covar,{
		if (is.null(input$covar) && input$analysis_type=='Two group with covariate') {
			shinyjs::hide("downloadsigdata")
		} else {
			shinyjs::show("downloadsigdata")
		}
	})
	observeEvent(input$covar,{
		if(!is.null(input$covar) && !is.null(design_rank())){
			shinyjs::hide("downloadsigdata")
		} 
	})
	observeEvent(input$analysis_property,{
		if (input$analysis_property!='') {
			shinyjs::show("downloadsigdata")
		} else {
			shinyjs::hide("downloadsigdata")
		}
	})

	##### show to ilincs
	
	# output$dynamiclink <- renderUI({
		# if (!is.null(input$analysis_property)&&!is.null(input$input_group1)&&!is.null(input$input_group2)) {
			# filepath<-paste('/mnt/raid/tmp/iLincs/signatures/',input$geo_acc2, '_signatureData.txt', sep='')
			# filename<-paste(input$geo_acc2, '_signatureData.txt', sep='')
			# write.table(sigdata()$signaturesData_final[,-1],filepath, row.names=FALSE, quote=F, sep="\t" )
			# link<-paste("http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",input$geo_acc2,"_signatureData.txt",sep='')
			# tags$a(href=link,paste0("Upload signatures to iLincs"),target="_blank",class="btn btn-primary", icon("paper-plane"))
			# )
		# } else {
			# return()
		# }
	# })
	
	output$dynamiclink <- renderUI({
		if(input$analysis_type=='Two group without covariate'){
			if (input$input_group2!='' && input$analysis_type=='Two group without covariate') {
				actionButton(inputId='ul2ilincs', label="Upload signatures to iLincs", icon = icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			} 
		} else if(input$analysis_type=='Multi group without covariate'){
			actionButton(inputId='ul2ilincs', label="Upload signatures to iLincs", icon = icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
		} else if(input$analysis_type=='Two group with covariate'){
			if (is.null(design_rank())) {
				actionButton(inputId='ul2ilincs', label="Upload signatures to iLincs", icon = icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			}
		}
	})

	observe({
		if (input$signature=='sig_data') {
			shinyjs::show("dynamiclink")
		} else {
			shinyjs::hide("dynamiclink")
		}
	})
	observeEvent(input$analysis_property,{
		if(input$analysis_property == ''){
			shinyjs::hide("dynamiclink")
		} else {
			shinyjs::show("dynamiclink")
		}
	})
	observeEvent(input$input_group2,{
		if (input$input_group2!='' && input$analysis_type!='Two group without covariate'){
			shinyjs::hide("dynamiclink")
		} else {
			shinyjs::show("dynamiclink")
		}
	})
	observeEvent(input$covar,{
		if (is.null(input$covar) && input$analysis_type=='Two group with covariate') {
			shinyjs::hide("dynamiclink")
		} else {
			shinyjs::show("dynamiclink")
		}
	})
	observeEvent(input$analysis_type,{
		if (input$analysis_type=='Multi group without covariate') {
			shinyjs::hide("dynamiclink")
		}
	})
	
	####### modal Heatmap in regular signature tab #######

	output$heat_typeui_mod_regsig <- renderUI({	
		radioButtons("heat_type_mod_regsig", h4("Type of heatmap"), c("Heatmap across all the samples"="all","Heatmap across the comparison samples"="notall"),
			selected="notall")
	})
	output$heat_prop_mod_regsig <- renderUI({	
		#col_names <- colnames(design_tab())[-1]
		col_names <- colnames(collapse_data()$sub_pdata)
		checkboxGroupInput("heat_property_mod_regsig", h4("Select properties"), choices  = col_names, 
			selected=col_names[which(col_names %in% input$analysis_property)])
	})
	output$cluster_meth_mod_regsig <- renderUI({	
		selectInput("clustering_method_mod_regsig", label = h4("Grouping samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected=
			if(input$cluster_ngenes_mod_regsig<10){
				"Euclidean distance"
			} else {
				"Pearson correlation"
			}
		)
	})	
	output$cluster_genes_mod_regsig <- renderUI({
		numericInput("cluster_ngenes_mod_regsig", label=h4('Number of top DE genes')
			,value= 100
			,min = 1, max = nrow(sigdata()$counts_all)
		)
	})

	plot_counts_mod_regsig <- reactive({
		counts <- if(input$heat_type_mod_regsig=="all"){
					sigdata()$counts_all
				} else {
					sigdata()$counts_sig
				}
		gene_sym <- collapse_data()$fdata[match(rownames(counts), rownames(collapse_data()$fdata)),,drop=F]
		rownames(counts) <- paste0(rownames(counts), " : ", gene_sym[,1])
		return(counts)
	})
	plotHeatFit_mod_regsig <- function(){
		complex_heatmap_sig (data_mat=plot_counts_mod_regsig(), metadata=collapse_data()$sub_pdata[colnames(plot_counts_mod_regsig()),,drop=F]
		, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
		, type="Fit in Screen")
	}
	data_dl_mod_regsig <- reactive({
		dge <- if(input$heat_type_mod_regsig=="all"){
					DGEList(counts=sigdata()$counts_all)
				} else {
					DGEList(counts=sigdata()$counts_sig)
				}
		dge_cpm <- cpm(dge, log=TRUE)
		expr_data<- data.matrix(dge_cpm)
		expr_data_in <- expr_data[1:input$cluster_ngenes_mod_regsig,]
		gene_sym <- collapse_data()$fdata[match(rownames(expr_data_in), rownames(collapse_data()$fdata)),,drop=F]
		dat2 <- data.frame(GeneSymbol=gene_sym[,1], expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	})

	output$fitscreen_mod_regsig = renderUI({
		validate(
			need(!is.null(input$heat_property_mod_regsig), "please select a property")
		)
		plotOutput("HeatFit_mod_regsig2", height="auto") %>% withSpinner(type=5)
	})
	output$HeatFit_mod_regsig2 = renderPlot({
		validate(
			need(!is.null(input$heat_property_mod_regsig), "please select a property")
		)  
		print(plotHeatFit_mod_regsig())
		
	}
		,height = reactive({
			if(!is.null(input$heat_property_mod_regsig)) {
				df <- collapse_data()$sub_pdata[,input$heat_property_mod_regsig,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- input$heat_property_mod_regsig
				x=sapply(df, function(x) length(unique(x)))
				nrow = max(x)/2
				if(input$cluster_ngenes_mod_regsig<=50)
					return((length(input$heat_property_mod_regsig)*25) + (160+ (input$cluster_ngenes_mod_regsig*5.50))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				else {
					return((length(input$heat_property_mod_regsig)*25) + 450 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				}
			} 
		})
	)
	output$downloadFit_mod_regsig <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_fit_mod_regsig=="png"){
				if (length(input$heat_property_mod_regsig)==1) {
					return(paste(datasetname,"_(",input$heat_property_mod_regsig,")_top ",input$cluster_ngenes_mod_regsig," DE genes_heatmap_fit.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property_mod_regsig, collapse="+"),")","_top ",input$cluster_ngenes_mod_regsig," DE genes_heatmap_fit.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes_mod_regsig," DE Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_fit_mod_regsig == "png") {
				png(file, width = 18, height = 16, units = "in",res=200)
				plotHeatFit_mod_regsig()
				dev.off()
			} else if(input$download_type_fit_mod_regsig == "csv"){
				dat2 <- data_dl_mod_regsig()
				write.csv(dat2, file)
			} else {
				return()
			}
		}
	)

	#### Scrollable

	plotHeatScroll_mod_regsig <- function(){
		complex_heatmap_sig (data_mat=plot_counts_mod_regsig(), metadata=collapse_data()$sub_pdata[colnames(plot_counts_mod_regsig()),,drop=F]
		, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
		, type="Scrollable")
	}
	output$scroll_mod_regsig = renderUI({
		validate(
			need(!is.null(input$heat_property_mod_regsig), "please select a property")
		)
		plotOutput("HeatScroll_mod_regsig", height="auto") %>% withSpinner(type=5)
	})
	output$HeatScroll_mod_regsig = renderPlot({
		validate(
			need(!is.null(input$heat_property_mod_regsig), "please select a property")
		)  
		print(plotHeatScroll_mod_regsig())
	}
		,height = reactive({
			if(input$cluster_ngenes_mod_regsig>=100) {
				value <- max((length(input$heat_property_mod_regsig)*25)+(10*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*24))
			} else {
				value <- (length(input$heat_property_mod_regsig)*25) + 400+(input$cluster_ngenes_mod_regsig*20)
			}
			return(value)				
		})
	)
	
	output$downloadScroll_mod_regsig <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scroll_mod_regsig=="png"){
				if (length(input$heat_property_mod_regsig)==1) {
					return(paste(datasetname,"_","(",input$heat_property_mod_regsig,")","_top ",input$cluster_ngenes_mod_regsig," DE genes_heatmap_scroll.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property_mod_regsig, collapse="+"),")","_top ",input$cluster_ngenes_mod_regsig," DE genes_heatmap_scroll.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes_mod_regsig," DE Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_scroll_mod_regsig == "png") {
				height = reactive({
					if(input$cluster_ngenes_mod_regsig>=1000) {
						value <- max((length(input$heat_property_mod_regsig)*25)+(4*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property_mod_regsig)*25)+(10*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				plotHeatScroll_mod_regsig()
				dev.off()
			} else{
				dat2 <- data_dl_mod_regsig()
				write.csv(dat2, file)
			}
		}
	)

	observe({
		if (input$heatmodal_regsig=="fitscrn_mod_regsig"){
			shinyjs::show("downloadFit_mod_regsig")
			shinyjs::show("download_type_fit_mod_regsig")
		} else {
			shinyjs::hide("downloadFit_mod_regsig")
			shinyjs::hide("download_type_fit_mod_regsig")
		}
	})
	observe({
		if (input$heatmodal_regsig=="scrolb_mod_regsig"){
			shinyjs::show("downloadScroll_mod_regsig")
			shinyjs::show("download_type_scroll_mod_regsig")
		} else {
			shinyjs::hide("downloadScroll_mod_regsig")
			shinyjs::hide("download_type_scroll_mod_regsig")
		}
	})
	observeEvent(input$user_analysis_prop,{
		if(input$user_analysis_prop == ''){
			shinyjs::hide("heatmap_mod_regsig")
		} else {
			shinyjs::show("heatmap_mod_regsig")
		}
	})

	observeEvent(input$analysis_type,{
		if (input$analysis_type=='Multi group without covariate') {
			shinyjs::show("heatmap_mod_regsig")
		}
	})
	observeEvent(input$analysis_type,{
		if (input$analysis_type=='Two group with covariate') {
			shinyjs::hide("heatmap_mod_regsig")
		}
	})
	observeEvent(input$input_group2,{
		if (input$input_group2=='' && input$analysis_type=='Two group without covariate') {
			shinyjs::hide("heatmap_mod_regsig")
		} else if (input$input_group2!='' && input$analysis_type!='Two group without covariate'){
			shinyjs::hide("heatmap_mod_regsig")
		} else {
			shinyjs::show("heatmap_mod_regsig")
		}
	})
	observeEvent(input$covar,{
		if (is.null(input$covar) && input$analysis_type=='Two group with covariate') {
			shinyjs::hide("heatmap_mod_regsig")
		} else {
			shinyjs::show("heatmap_mod_regsig")
		}
	})
	observeEvent(input$covar,{
		if(!is.null(input$covar) && !is.null(design_rank())){
			shinyjs::hide("heatmap_mod_regsig")
		} 
	})
	observeEvent(input$analysis_property,{
		if (input$analysis_property!='') {
			shinyjs::show("heatmap_mod_regsig")
		} else {
			shinyjs::hide("heatmap_mod_regsig")
		}
	})

	
	############## User made experimental design ##############
	
	metadata_user_exp <- reactive({
		metadata <- collapse_data()$sub_pdata
		metadata <- metadata[, input$user_analysis_prop,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- input$user_analysis_prop
		rownames(metadata) <- paste0(rownames(metadata), " : ",metadata[,input$user_analysis_prop])
		return(metadata)
	})

	output$user_meta_prop <- renderUI ({
		validate(
			need(input$signature=="exp_desg","")
		)
		metadata <- collapse_data()$sub_pdata
		col_names <- colnames(metadata)
		selectInput(inputId="user_analysis_prop",label=h4('Variable of interest'),
			choices=c('',col_names), selected=NULL
		)
	})
	output$user_g1_ui <- renderUI ({
		validate(
			need(input$signature=="exp_desg","")
		)
		dropdownButton(label = "Select experimental samples", status = "primary", width = 390 #250
			,tags$div(class = "container", 
				checkboxGroupInput(inputId = "user_g1", label = "", 
				#choices = as.character(metadata_user_exp()[,1]))
				choices = rownames(metadata_user_exp()))
			)
		)
	})
	output$user_g2_ui <- renderUI ({
		validate(
			need(input$signature=="exp_desg",""),
			need(input$user_analysis_prop!='',"")
		)
		x <- input$user_g1
		avail <- rownames(metadata_user_exp())[!(rownames(metadata_user_exp()) %in% x)]
		dropdownButton(label = "Select control samples", status = "primary", width = 390 #220
			,tags$div(class = "container", 
				checkboxGroupInput(inputId = "user_g2", label = "", width="100%",
				choices = avail)
			)
		)
	})		
	output$user_ana_type <- renderUI({
		validate(
			need(input$signature=="exp_desg","")
			,need(input$user_analysis_prop!='',"")
		)
		selectInput(inputId="user_analysis_type",label=h4('Type of comparison'),
			choices= 
				if(ncol(collapse_data()$sub_pdata)>1){
					c('Two group without covariate','Two group with covariate')
				} else {
					'Two group without covariate'
				},
			selected= 'Two group without covariate'
		)
	})
	output$user_cov_ui <- renderUI({
		validate(
			need(input$signature=="exp_desg","")
			,need(input$user_analysis_prop!='',"")
			,need(input$user_analysis_type=="Two group with covariate","")
			#,need(input$user_g1!="","")
			#,need(input$user_g2!="","")
		)
		selectInput(inputId="user_cov",label=h4('Select covariates'), 
			choices= colnames(collapse_data()$sub_pdata)[which(colnames(collapse_data()$sub_pdata)!=input$user_analysis_prop)]
			, multiple = TRUE, selected= NULL
		)
	})
	output$gen_sig_btn_ui <- renderUI({
		validate(
			need(input$signature=="exp_desg","")
			,need(input$user_analysis_prop!='',"")
			,need(input$user_g1!="","")
			,need(input$user_g2!="","")
			,need(all(design_user_levels()>=2), "")
		)
		if(input$user_analysis_type=="Two group without covariate"){
			actionButton("gen_sig_btn", "Generate signature",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
		} else if(input$user_analysis_type=="Two group with covariate" & !is.null(input$user_cov)){
			if(is.null(design_user_rank())){
				actionButton("gen_sig_btn", "Generate signature",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
			} else {
				return()
			}
		} else {
			return()
		}
	})	
		
	metadata_user_exp_cov <- reactive({
		metadata <- collapse_data()$sub_pdata
		metadata <- metadata[, c(input$user_analysis_prop,input$user_cov),drop=FALSE]
		rownames(metadata) <- paste0(rownames(metadata), " : ",metadata[,input$user_analysis_prop])
		return(metadata)
	})
	design_tab <- reactive({
		if(input$user_analysis_type=="Two group without covariate"){
			dat1 <- metadata_user_exp()[rownames(metadata_user_exp()) %in% input$user_g1,,drop=FALSE]
			dat1$"Selected groups"[rownames(dat1)==input$user_g1] <- "experimental"
			dat2 <- metadata_user_exp()[rownames(metadata_user_exp()) %in% input$user_g2,,drop=FALSE]
			dat2$"Selected groups"[rownames(dat2)==input$user_g2] <- "control"			
			dat <- rbind(dat1,dat2)
			dat <- dat[, c("Selected groups", input$user_analysis_prop)]
			return(dat)
		} else if(input$user_analysis_type=="Two group with covariate"){
			dat1 <- metadata_user_exp_cov()[rownames(metadata_user_exp_cov()) %in% input$user_g1,,drop=FALSE]
			dat1$"Selected groups"[rownames(dat1)==input$user_g1] <- "experimental"
			dat2 <- metadata_user_exp_cov()[rownames(metadata_user_exp_cov()) %in% input$user_g2,,drop=FALSE]
			dat2$"Selected groups"[rownames(dat2)==input$user_g2] <- "control"
			dat <- rbind(dat1,dat2)
			dat <- dat[ ,c("Selected groups", input$user_analysis_prop, input$user_cov)]
			return(dat)
		} else {
			return()
		}
	})
	design_user_levels <- reactive({
		if(input$user_analysis_type=="Two group without covariate"){
			metad <- design_tab()
			#design_data <- metad[which(metad[,"Selected groups"]=="experimental" | metad[,"Selected groups"]=="control"),"Selected groups",drop=FALSE]
			design_data <- metad[,-1,drop=F]
			design_data[] <- lapply(design_data, factor)
			nlev <- as.numeric(sapply(design_data[,sapply(design_data, is.factor)], nlevels))
			return(nlev)
		} else if(input$user_analysis_type=="Two group with covariate"){
			metad <- design_tab()
			#analysis_metadata <- metad[which(metad[,"Selected groups"]=="experimental" | metad[,"Selected groups"]=="control"),c("Selected groups",input$user_cov),drop=FALSE]
			#group <- analysis_metadata[,"Selected groups", drop=F]
			#covar <- analysis_metadata[,input$user_cov, drop=F]
			#design_data <- cbind(covar, group)
			design_data <- metad[,-1,drop=F]
			design_data[] <- lapply(design_data, factor)
			nlev <- as.numeric(sapply(design_data[,sapply(design_data, is.factor)], nlevels))
			return(nlev)
		} else {
			return()
		}
	})
	design_user_rank <- reactive({
		metad <- design_tab()
		analysis_metadata <- metad[which(metad[,"Selected groups"]=="experimental" | metad[,"Selected groups"]=="control"),c("Selected groups",input$user_cov),drop=FALSE]
		group <- analysis_metadata[,"Selected groups", drop=F]
		covar <- analysis_metadata[,input$user_cov, drop=F]
		design_data <- cbind(covar, group)
		design_data[] <- lapply(design_data, factor)
		if(all(design_user_levels()>=2)){
			design <- model.matrix(~., data=design_data)
			ne <- nonEstimable(design)	
			return(ne)
		} else {
			return()
		}
	})

	output$designtable = DT::renderDataTable({
		validate(
			need(any(input$signature=="exp_desg"),"")
		)
		datatable( design_tab()
			,selection = 'none'
			,options = list(pageLength = if(nrow(design_tab())>=10) {10} else {nrow(design_tab())}, scrollX = TRUE)
			,rownames= TRUE
		)
	}, server = FALSE) 

	output$user_rank_warn <- renderText({
		validate(
			need(input$user_analysis_type=='Two group with covariate',"")
			,need(input$user_cov!="","")
			,need(all(design_user_levels()>=2), "")
			,need(!is.null(design_user_rank()), "")
		)	
		paste("Design matrix is not of full rank. The following coefficients are not estimable:\n", paste(design_user_rank(), collapse = " "))
	})
	output$user_level_warn <- renderText({
		validate(
			need(any(any(design_user_levels()<2)), "")
		)	
		paste("All the variables need to have 2 or more levels")
	})
	
	sigdata_user <- reactive({
		
		filepath <- paste('/mnt/raid/tmp/iLincs/signatures/',toupper(input$geo_acc2), '_signatureData.txt', sep='')
		analysis_genes <- cbind(Ensembl_ID=rownames(collapse_data()$fdata), Gene_symbol=collapse_data()$fdata[,1])
		metadata <- design_tab()
		rownames(metadata) <- sapply(strsplit(rownames(design_tab())," : "), `[`, 1)
		if(input$user_analysis_type=='Two group without covariate'){
			signature_data= analysis_2g_wocov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property="Selected groups",
				group1="experimental", group2="control")
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,3], Significance_pvalue=signature_data$signaturesData_final[,4], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$user_analysis_type=='Two group with covariate'){
			signature_data= analysis_2g_withcov(counts=collapse_data()$counts, metadata=metadata, genes=analysis_genes, property="Selected groups",
				covariate=input$user_cov, group1="experimental", group2="control")
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,3], Significance_pvalue=signature_data$signaturesData_final[,4], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else {
			return()
		}
		return(signature_data)
	})

	observeEvent(input$gen_sig_btn, {
		if (input$gen_sig_btn) {
			updateTabsetPanel(session, "expr_design", selected = "user_exp_sig")
		}
    })
	output$sigtable <- DT::renderDataTable({
		validate(
			need(input$gen_sig_btn,"")
		)
		datatable(
			sigdata_user()$signaturesData_final
			,selection = 'none'
			,rownames = FALSE
			, filter = 'top' 
			,options = list(pageLength = 10, scrollX = TRUE
				,columnDefs = list(list(className = 'dt-center', targets ="_all"))
			)
		)
	})		
	output$download_user_sigdata <- downloadHandler(
		filename = function() { 
			paste(toupper(input$geo_acc2), '_signatureData.csv', sep='') 
		},
		content = function(file) {
			write.csv(sigdata_user()$signaturesData_final[,-1], file, row.names=FALSE )
		}
	)	

	output$user_dynamiclink <- renderUI({
		if(input$user_analysis_type=='Two group without covariate'){
			if (input$user_g2!='' && input$user_analysis_type=='Two group without covariate') {
				actionButton(inputId='user_ul2ilincs', label="Upload signatures to iLincs", icon = icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			} else {
				return()
			}
		} else if(input$user_analysis_type=='Two group with covariate'){
			#if (is.null(design_user_rank()) & all(design_user_levels()>=2)) {
				actionButton(inputId='user_ul2ilincs', label="Upload signatures to iLincs", icon = icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			# } else {
				# return()
			# }
		}
	})
	
	observeEvent(input$user_analysis_prop,{
		if(input$user_analysis_prop == ''){
			shinyjs::hide("user_dynamiclink")
		} else {
			shinyjs::show("user_dynamiclink")
		}
	})
	observe({
		shinyjs::hide("gen_sig_btn")
	})
	observeEvent(input$user_cov,{
		if(!is.null(input$user_cov) && !is.null(design_user_rank())){
			shinyjs::hide("gen_sig_btn")
		} else {
			shinyjs::show("gen_sig_btn")
		}
	})

	observe({
		hide(selector = "#expr_design li a[data-value=user_exp_sig]")
	})
	observeEvent(input$gen_sig_btn,{
		if (!is.null(input$gen_sig_btn)) {
			shinyjs::show(selector = "#expr_design li a[data-value=user_exp_sig]")
		} else {
			shinyjs::hide(selector = "#expr_design li a[data-value=user_exp_sig]")
		}
	})
	observeEvent(input$user_cov, {
		if(input$user_cov!="" & !is.null(design_user_rank())){
			shinyjs::hide("download_user_sigdata")
		} else {
			shinyjs::show("download_user_sigdata")
		}
	})
	observeEvent(input$expr_design,{
		if(input$expr_design=="user_exp_sig"){
			shinyjs::show("download_user_sigdata")
		}	
	})
	observeEvent(input$expr_design,{
		if(input$expr_design=="user_exp_sig"){
			shinyjs::disable("user_analysis_type")
			shinyjs::disable("user_analysis_prop")
			shinyjs::disable("user_cov")
		} else {
			shinyjs::enable("user_analysis_type")
			shinyjs::enable("user_analysis_prop")
			shinyjs::enable("user_cov")
		}
	})
	observeEvent(input$expr_design,{
		output$user_method <- renderUI({
			if (input$signature=="exp_desg") {
				includeMarkdown("www/user_design.Rmd")
			}
		})
	})


############ modal Heatmap ############

	output$heat_typeui_mod <- renderUI({	
		radioButtons("heat_type_mod", h4("Type of heatmap"), c("Heatmap across all the samples"="all", "Heatmap across the comparison samples"="notall"),
			selected="notall")
	})
	output$heat_prop_mod <- renderUI({	
		col_names <- colnames(collapse_data()$sub_pdata)
		checkboxGroupInput("heat_property_mod", h4("Select properties"), choices  = col_names, selected=input$user_analysis_prop)
	})
	output$cluster_meth_mod <- renderUI({	
		selectInput("clustering_method_mod", label = h4("Grouping samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected=
			if(input$cluster_ngenes_mod<10){
				"Euclidean distance"
			} else {
				"Pearson correlation"
			}
		)
	})	
	output$cluster_genes_mod <- renderUI({
		numericInput("cluster_ngenes_mod", label=h4('Number of top DE genes')
			,value= 100
			,min = 1, max = nrow(sigdata_user()$counts_all)
		)
	})

	plot_counts_mod <- reactive({
		counts <- if(input$heat_type_mod=="all"){
			sigdata_user()$counts_all
		} else {
			sigdata_user()$counts_sig
		}
		gene_sym <- collapse_data()$fdata[match(rownames(counts), rownames(collapse_data()$fdata)),,drop=F]
		rownames(counts) <- paste0(rownames(counts), " : ", gene_sym[,1])
		return(counts)
	})
	plotHeatFit_mod <- function(){
		complex_heatmap_sig (data_mat=plot_counts_mod(), metadata=collapse_data()$sub_pdata[colnames(plot_counts_mod()),,drop=F]
		, property=input$heat_property_mod,	input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod, type="Fit in Screen")
	}
	data_dl_mod <- reactive({
		dge <- if(input$heat_type_mod=="all"){
				DGEList(counts=sigdata_user()$counts_all)
			} else {
				DGEList(counts=sigdata_user()$counts_sig)
			}
		dge_cpm <- cpm(dge, log=TRUE)
		expr_data<- data.matrix(dge_cpm)
		expr_data_in <- expr_data[1:input$cluster_ngenes_mod,]
		gene_sym <- collapse_data()$fdata[match(rownames(expr_data_in), rownames(collapse_data()$fdata)),,drop=F]
		dat2 <- data.frame(GeneSymbol=gene_sym[,1], expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	})

	output$fitscreen_mod = renderUI({
		validate(
			need(!is.null(input$heat_property_mod), "please select a property")
		)
		plotOutput("HeatFit_mod", height="auto") %>% withSpinner(type=5)
	})
	output$HeatFit_mod = renderPlot({
		validate(
			need(!is.null(input$heat_property_mod), "please select a property")
		)  
		print(plotHeatFit_mod())
		
	}
		,height = reactive({
			if(!is.null(input$heat_property_mod)) {
				df <- collapse_data()$sub_pdata[,input$heat_property_mod,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- input$heat_property_mod
				x=sapply(df, function(x) length(unique(x)))
				nrow = max(x)/2
				if(input$cluster_ngenes_mod<=50)
					return((length(input$heat_property_mod)*25) + (160+ (input$cluster_ngenes_mod*5.50))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				else {
					return((length(input$heat_property_mod)*25) + 450 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
				}
			} 
		})
	)
	output$downloadFit_mod <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_fit_mod=="png"){
				if (length(input$heat_property_mod)==1) {
					return(paste(datasetname,"_(",input$heat_property_mod,")_top ",input$cluster_ngenes_mod," DE genes_heatmap_fit.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property_mod, collapse="+"),")","_top ",input$cluster_ngenes_mod," DE genes_heatmap_fit.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes_mod," DE Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_fit_mod == "png") {
				png(file, width = 18, height = 16, units = "in",res=200)
				plotHeatFit_mod()
				dev.off()
			} else if(input$download_type_fit_mod == "csv"){
				dat2 <- data_dl_mod()
				write.csv(dat2, file)
			} else {
				return()
			}
		}
	)

	#### Scrollable

	plotHeatScroll_mod <- function(){
		complex_heatmap_sig (data_mat=plot_counts_mod(), metadata=collapse_data()$sub_pdata[colnames(plot_counts_mod()),,drop=F]
		, property=input$heat_property_mod, input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod, type="Scrollable")
	}
	output$scroll_mod = renderUI({
		validate(
			need(!is.null(input$heat_property_mod), "please select a property")
		)
		plotOutput("HeatScroll_mod", height="auto") %>% withSpinner(type=5)
	})
	output$HeatScroll_mod = renderPlot({
		validate(
			need(!is.null(input$heat_property_mod), "please select a property")
		)  
		print(plotHeatScroll_mod())
	}
		,height = reactive({
			if(input$cluster_ngenes_mod>=100) {
				value <- max((length(input$heat_property_mod)*25)+(10*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*24))
			} else {
				value <- (length(input$heat_property_mod)*25) + 400+(input$cluster_ngenes_mod*20)
			}
			return(value)				
		})
	)
	
	output$downloadScroll_mod <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scroll_mod=="png"){
				if (length(input$heat_property_mod)==1) {
					return(paste(datasetname,"_","(",input$heat_property_mod,")","_top ",input$cluster_ngenes_mod," DE genes_heatmap_scroll.png",sep=""))
				} else {
					return(paste(datasetname,"_","(",paste(input$heat_property_mod, collapse="+"),")","_top ",input$cluster_ngenes_mod," DE genes_heatmap_scroll.png",sep=""))
				}
			} else {
				return(paste(datasetname,": Top_",input$cluster_ngenes_mod," DE Genes.csv",sep=""))
			}
		},
		content = function(file) {
			if(input$download_type_scroll_mod == "png") {
				height = reactive({
					if(input$cluster_ngenes_mod>=1000) {
						value <- max((length(input$heat_property_mod)*25)+(4*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property_mod)*25)+(10*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				plotHeatScroll_mod()
				dev.off()
			} else{
				dat2 <- data_dl_mod()
				write.csv(dat2, file)
			}
		}
	)

	observe({
		if (input$heatmodal=="fitscrn_mod"){
			shinyjs::show("downloadFit_mod")
			shinyjs::show("download_type_fit_mod")
		} else {
			shinyjs::hide("downloadFit_mod")
			shinyjs::hide("download_type_fit_mod")
		}
	})
	observe({
		if (input$heatmodal=="scrolb_mod"){
			shinyjs::show("downloadScroll_mod")
			shinyjs::show("download_type_scroll_mod")
		} else {
			shinyjs::hide("downloadScroll_mod")
			shinyjs::hide("download_type_scroll_mod")
		}
	})
	observeEvent(input$user_analysis_prop,{
		if(input$user_analysis_prop == ''){
			shinyjs::hide("heatmap_mod")
		} else {
			shinyjs::show("heatmap_mod")
		}
	})



################################################ Processing users data nav page ################################################
	
	
	observeEvent(input$user_geo,{
		info2 = unlist(lapply(datafr$x[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		#if(toupper(input$user_geo)=='' || toupper(input$user_geo) %in% datafr$x[,1] || toupper(input$user_geo) %in% filenames() || !toupper(input$user_geo) %in% geo_QuerynBrowser_merged[,1] || nchar(input$user_geo)<7 || substr(toupper(input$user_geo), 1, 3) != "GSE" || !grepl("^[[:digit:]]*$", substr(toupper(input$user_geo), 4, nchar(input$user_geo)))) {
		if(toupper(input$user_geo)=='' || toupper(input$user_geo) %in% info2 || toupper(input$user_geo) %in% filenames() || !toupper(input$user_geo) %in% geo_QuerynBrowser_merged_jan28_2018[,1] || nchar(input$user_geo)<7 || substr(toupper(input$user_geo), 1, 3) != "GSE" || !grepl("^[[:digit:]]*$", substr(toupper(input$user_geo), 4, nchar(input$user_geo)))) {
			shinyjs::hide("start_process")
		} else {
			shinyjs::show("start_process")
		}
	})

	output$warn <- renderText({
		validate(
			need(toupper(input$user_geo)!='',"")
		)
		paste0("This dataset has already been processed")
	})
	output$warn2 <- renderText({
		validate(
			need(toupper(input$user_geo)!='',"")
		)	
		paste0("Please enter a valid GEO series accession")
	})
	output$warn3 <- renderText({
		validate(
			need(toupper(input$user_geo)!='',"")
		)	
		paste0("This dataset does not exist in the our database")
	})
	output$warn4 <- renderText({
		validate(
			need(toupper(input$user_geo)!='',"")
		)	
		paste0("You have already requested this dataset. To see the status please press the 'Processing console' button")
	})
	
	observeEvent(input$user_geo,{
		info2 = unlist(lapply(datafr$x[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		#if(toupper(input$user_geo) %in% datafr$x[,1]) {
		if(toupper(input$user_geo) %in% info2) {
			shinyjs::show("warn")
		} else {
			shinyjs::hide("warn")
		}
	})
	observeEvent(input$user_geo,{
		if(nchar(input$user_geo)<7 | substr(toupper(input$user_geo), 1, 3) != "GSE" | !grepl("^[[:digit:]]*$", substr(toupper(input$user_geo), 4, nchar(input$user_geo)))) {
			shinyjs::show("warn2")
		} else {
			shinyjs::hide("warn2")
		}
	})
	observeEvent(input$user_geo,{
		load("/opt/raid10/genomics/naim/Thesis/geo/geo_QuerynBrowser_merged_jan28_2018.RData")
		if(substr(toupper(input$user_geo), 1, 3) == "GSE" && nchar(input$user_geo)>=7 && grepl("^[[:digit:]]*$", substr(toupper(input$user_geo), 4, nchar(input$user_geo))) && !toupper(input$user_geo) %in% geo_QuerynBrowser_merged_jan28_2018[,1]){
			shinyjs::show("warn3")
		} else {
			shinyjs::hide("warn3")
		}
	})
	observeEvent(input$user_geo,{
		if(toupper(input$user_geo) %in% filenames()) {
			shinyjs::show("warn4")
			shinyjs::hide("warn3")
		} else {
			shinyjs::hide("warn4")
		}
	})


	observe({
		hide(selector = "#grin li a[data-value=out_cons]")
	})
	observeEvent(input$start_process,{
		if (!is.null(input$start_process) && input$user_geo!='') {
			updateTabsetPanel(session, "grin", selected = "out_cons");
			shinyjs::show(selector = "#grin li a[data-value=out_cons]")
		} else {
			shinyjs::hide(selector = "#grin li a[data-value=out_cons]")
		}
	})
	observeEvent(input$process_log,{
		if (!is.null(input$process_log)) {
			updateTabsetPanel(session, "grin", selected = "out_cons");
			shinyjs::show(selector = "#grin li a[data-value=out_cons]")
		} else {
			shinyjs::hide(selector = "#grin li a[data-value=out_cons]")
		}
	})

	observeEvent(input$start_process,{
		if (input$start_process) {
			if(input$user_geo!=''){
				updateTabsetPanel(session, "grin", selected = "out_cons")
						
				if(!file.exists(paste0("/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", toupper(input$user_geo), ".txt")) & 
					input$user_geo!='' & 
					!dir.exists(paste0("/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", toupper(input$user_geo)))){
					system(paste0("touch /opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", toupper(input$user_geo), ".txt"))
				}
				
				filenames <- function(){
					details <- file.info(list.files(path="/opt/raid10/genomics/naim/Thesis/geo/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
					files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
					x=gsub("^.*/", "", sub(".txt","",files))
					return(x)
				}
				create_dir <- function(){
					system(paste0("mkdir /opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", filenames()[1]))
					system(paste0("chmod -R 777 /opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", filenames()[1]))
					dir_name <- list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE)
					return(dir_name)
				}

				if(length(list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE))==0) {
					dir_name= create_dir()
				} else {
					dir_name= list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE)
				}

				updateTextInput(session, 'geo_current', value = dir_name)
				updateSelectInput(session, 'geo_next', label="Waiting to process", choices=filenames(), selected=filenames()[1])
									
				output$steps <- renderText({
					paste0("There are 5 steps to be completed")
				})
							
				logfilename <- paste0("/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/",dir_name,"/log.txt")
				fileReaderData <- reactiveFileReader(500, session, logfilename, readLines)
				output$process_log <- renderText({
					validate(
						need(input$geo_current!="","")
					)
					
					text <- fileReaderData()
					text[is.na(text)] <- ""
					paste(text, collapse = '\n')
				})
			}
		}
    })	
	
	observeEvent(input$process_log,{
		if (input$process_log) {
			updateTabsetPanel(session, "grin", selected = "out_cons")

			filenames <- function(){
				details <- file.info(list.files(path="/opt/raid10/genomics/naim/Thesis/geo/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
				files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
				x=gsub("^.*/", "", sub(".txt","",files))
				return(x)
			}
			dir_created <- function(){
				dir_name <- list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE)
				return(dir_name)
			}

			if(length(list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE))==0) {
				dir_name= dir_created()
			} else {
				dir_name= list.dirs(path = "/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/", full.names = F, recursive = FALSE)
			}

			updateTextInput(session, 'geo_current', value = dir_name)
			updateSelectInput(session, 'geo_next', label="Waiting to process", choices=filenames(), selected=filenames()[1])
								
			output$steps <- renderText({
				paste0("There are 5 steps to be completed")
			})
						
			logfilename <- paste0("/opt/raid10/genomics/naim/Thesis/geo/user_geo_request/",dir_name,"/log.txt")
			fileReaderData <- reactiveFileReader(500, session, logfilename, readLines)
			output$process_log <- renderText({
				validate(
					need(input$geo_current!="","")
				)
				
				text <- fileReaderData()
				text[is.na(text)] <- ""
				paste(text, collapse = '\n')
			})
		}
    })	


	
})