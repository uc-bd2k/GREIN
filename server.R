
shinyServer(function(input, output, session) {

	source("R/sub_metadata.R", local = TRUE)
	source("R/dropdownButton.R", local = TRUE)
	source("R/complex_heatmap.R", local = TRUE)
	source("R/interactive_heatmap.R", local = TRUE)
	source("R/pca.R", local = TRUE)
	source("R/analysis.R", local = TRUE)
	source("R/complex_heatmap_sig.R", local = TRUE)
	
	####======== Help tab ========####
	pdfurl <- "step-by-step.pdf"
	output$help_pdf <- renderText({
		return(paste('<iframe style="height:600px; width:100%" src="', pdfurl, '"></iframe>', sep = ""))
	})
		
	
	####======== main data for GRIN datatable and search in section 1 ========####
	datafr2 <- as.data.frame(read_feather(paste0("data/GREIN_datatable.csv")),stringsAsFactors=F)
	all_samples <- as.data.frame(read_feather(paste0("data/all_samples_feather.csv")),stringsAsFactors=F, check.names=F)
	samples_ont <- as.data.frame(read_feather(paste0("data/samples_ontology_full_feather.csv")),stringsAsFactors=F, check.names=F)
	geo_QuerynBrowser_merged_jan28_2018 <- as.data.frame(read_feather(paste0("data/geo_QuerynBrowser_merged_jan28_2018.csv")),stringsAsFactors=F)	
	datatable_to_use <- as.data.frame(read_feather(paste0("data/datatable_to_use.csv")),stringsAsFactors=F)
	
	datafr <- reactiveValues()
	observe({
		user_geo = trimws(input$user_geo, which ="both")
		datafr$x <- if(input$user_geo=="" && input$sample_ont_ser==""){ 
					datafr2
				} else if (input$user_geo!="" && input$sample_ont_ser==""){
					info2 = unlist(lapply(datafr2[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
					datafr2[grep(toupper(user_geo),info2,ignore.case=TRUE),,drop=F]
				} else if(input$user_geo=="" && input$sample_ont_ser!=""){
					info2 = unlist(lapply(datafr2[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
					user_search <- c(samples_ont[grep(input$sample_ont_ser,samples_ont[,2]),1], samples_ont[grep(input$sample_ont_ser,samples_ont[,3]),1]) 					
					geo_acc <- unique(all_samples[which(all_samples[,"Sample"] %in% user_search),"gse_acc"])
					datafr2[which(info2 %in% geo_acc),,drop=F]
				} else {
					return()
				}
	})

	####======== Section 2 pop-ups and warnings ========####
	filenames <- reactive({
		details <- file.info(list.files(path="data/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
		files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
		x=gsub("^.*/", "", sub(".txt","",files))
		dir_name <- list.dirs(path = "data/user_geo_request", full.names = F, recursive = FALSE)
		return(c(dir_name,x))
	})
	observeEvent(input$user_geo, {
		if(!is.null(input$start_process)){
			updateSelectInput(session, 'geo_next_lp', label=" ", choices=filenames(), selected=filenames()[1])
		}
	})
	observeEvent(input$user_geo,{
		info2 = unlist(lapply(datafr$x[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		user_geo = trimws(input$user_geo, which ="both")
		if(toupper(user_geo)=='' || toupper(user_geo) %in% info2 || toupper(user_geo) %in% filenames() || nchar(user_geo)<7 || substr(toupper(user_geo), 1, 3) != "GSE" || !grepl("^[[:digit:]]*$", substr(toupper(user_geo), 4, nchar(user_geo)))) {
			shinyjs::hide("start_process")
		} else {
			shinyjs::show("start_process")
		}
	})

	output$warn <- renderText({
		validate(
			need(input$user_geo!='',"")
		)
		paste0("This dataset has already been processed. Please see the following table.")
	})
	output$warn2 <- renderText({
		validate(
			need(input$user_geo!='',"")
		)	
		paste0("Please enter a valid GEO series accession.")
	})
	output$warn3 <- renderText({
		validate(
			need(input$user_geo!='',"")
		)	
		paste0("This dataset does not exist in the GEO database.")
	})
	output$warn4 <- renderText({
		validate(
			need(input$user_geo!='',"")
		)	
		paste0("This dataset is currently in the processing que. To see the status please press the 'Processing console' button")
	})
	output$warn5 <- renderText({
		validate(
			need(input$user_geo!='',"")
		)	
		paste0("This dataset is not yet processed. You can initialize the processing by clicking the 'Start processing' button.")
	})
	
	observeEvent(input$user_geo,{
		info2 = unlist(lapply(datafr$x[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		user_geo = trimws(input$user_geo, which ="both")
		if(toupper(user_geo) %in% info2) {
			shinyjs::show("warn")
		} else {
			shinyjs::hide("warn")
		}
	})
	observeEvent(input$user_geo,{
		user_geo = trimws(input$user_geo, which ="both")
		if(nchar(user_geo)<7 | substr(toupper(user_geo), 1, 3) != "GSE" | !grepl("^[[:digit:]]*$", substr(toupper(user_geo), 4, nchar(user_geo)))) {
			shinyjs::show("warn2")
		} else {
			shinyjs::hide("warn2")
		}
	})
	observeEvent(input$user_geo,{
		user_geo = trimws(input$user_geo, which ="both")
		if(substr(toupper(user_geo), 1, 3) == "GSE" && nchar(user_geo)>=7 && grepl("^[[:digit:]]*$", substr(toupper(user_geo), 4, nchar(user_geo)))){
			shinyjs::show("warn3")
		} else {
			shinyjs::hide("warn3")
		}
	})
	observeEvent(input$user_geo,{
		user_geo = trimws(input$user_geo, which ="both")
		if(toupper(user_geo) %in% filenames()) {
			shinyjs::show("warn4")
			shinyjs::hide("warn3")
		} else {
			shinyjs::hide("warn4")
		}
	})
    observeEvent(input$user_geo,{		
        info2 = unlist(lapply(datafr$x[,1], function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		user_geo = trimws(input$user_geo, which ="both")
        if(!toupper(user_geo) %in% info2 & !toupper(user_geo) %in% filenames()){
            shinyjs::show("warn5")
        } else {
            shinyjs::hide("warn5")
        }
    })
	output$warn_ont <- renderText({
		validate(
			need(input$sample_ont_ser!='',"")
			,need(nrow(datafr$x)>0,"")
		)
		paste0("Click any dataset in the following table and see ontologies in the metadata tab.")
	})
	output$warn_ont2 <- renderText({
		validate(
			need(input$sample_ont_ser!='',"")
			,need(nrow(datafr$x)<1,"")
		)
		paste0("Nothing found!")
	})

	observe({
		if (input$sample_ont_ser!="" && input$user_geo=="") {
			shinyjs::disable("user_geo")
		} else if(input$user_geo!="" && input$sample_ont_ser==""){
			shinyjs::disable("sample_ont_ser")
		} else {
			shinyjs::enable("user_geo")
			shinyjs::enable("sample_ont_ser")
		}
	})

	####======== Section 3: study status ========####
	
	output$study_stat <- renderPlotly({
		x <- geo_QuerynBrowser_merged_jan28_2018[,c(1,4)]
		x$study_status <- ifelse(geo_QuerynBrowser_merged_jan28_2018[,1] %in% datatable_to_use[,1], "Processed", "In progress")
		
		status <- c('Processed','In progress')
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

	####======== Section 3: sample status ========####
	output$sample_stat <- renderPlotly({
		x=datatable_to_use[,c(2,3)]
		species <- c('Homo sapiens','Mus musculus', 'Rattus norvegicus') 
		samples <- c(length(x[which(x[,2]=='Homo sapiens'),1]), length(x[which(x[,2]=='Mus musculus'),1]), length(x[which(x[,2]=='Rattus norvegicus'),1]))
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
	
	####======== Section 3: sample density&histogram ========####
	output$samplesize_den <- renderPlotly({

		df <- data.frame(xlog=log10(as.numeric(datafr2[,2])),x= as.numeric(datafr2[,2]))
		p <- ggplot(df, aes(x=xlog, text=paste0("no. of samples: ",x))) + 
		  geom_histogram(alpha = 0.7, fill = "darkblue") + 
		  #geom_density(fill = "white", alpha = 0.5) + 
		  theme(panel.background = element_rect(fill = '#ffffff'), plot.title = element_text(size = 12, vjust=3),
			axis.text.y = element_text(size=7)) + 
		  scale_x_continuous(name="number of samples",breaks=c(1,2,3,3.8),labels=c("10", "100", "1000","6000")) +
		  labs(title = "Sample size distribution of the processed datasets\n\n")
		ggplotly(p, tooltip=c("text","y")) %>% config(displayModeBar = F) 
	})

	####======== Section 4: main data table ========####
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
		shinyjs::hide(selector = "#grin li a[data-value=explore]")
	})
	observeEvent(input$datatable_cell_clicked,{
		info = input$datatable_cell_clicked
		if(!is.numeric(info$value)){
			info2 = unlist(lapply(info$value, function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		} else {
			info2=NULL
		}
		if (!is.null(info2) && info$col == 0) {
			shinyjs::show(selector = "#grin li a[data-value=explore]")
		} else {
			shinyjs::hide(selector = "#grin li a[data-value=explore]")
		}
	})
		
	observeEvent(input$datatable_cell_clicked, {
		info = input$datatable_cell_clicked
		if(!is.numeric(info$value)){
			info2 = unlist(lapply(info$value, function(x) unlist(strsplit(unlist(strsplit(x, ">"))[2], "</button"))[1]))
		} else {
			info2=NULL
		}
		if (is.null(info2) || info$col != 0) return()
		updateTabsetPanel(session, "grin", selected = "explore")
		updateTextInput(session, 'geo_acc', value= info2)
	})

	############### organize eset as master data for all purpose ##############
	
	collapse_data <- reactive({
		load(paste0("data/",toupper(input$geo_acc),"/eset.RData"))
		pdata <- pData(eset)
		x <- exprs(eset)
		f <- fData(eset)
		rownames(x) <- f[,1]

		id <- unlist(lapply(rownames(pdata), function(x) unlist(strsplit(x, "_"))[2]))
		if("geo_accession" %in% names(pdata)) {
			if(length(unique(id))!= nrow(pdata)){
				pdata$ids <- id
				groupby <- pdata$ids
				if (!is.factor(groupby)) 
					groupby <- factor(groupby)
				groupby <- droplevels(groupby)
				sp <- split(seq(along = groupby), groupby)
				countsTable <- sapply(sp, function(i) rowSums(round(x)[,i, drop = FALSE]))
				colsToKeep <- sapply(sp, `[`, 1)
				pdata <- pdata[colsToKeep,]
				rownames(pdata) <- pdata$ids
				colnames(countsTable) <- rownames(pdata)
				pdata <- pdata[,-which(colnames(pdata)=="ids")]	
				meta <- data.frame(rowname=rownames(pdata), pdata, stringsAsFactors=FALSE,check.names=F)
				
				sub_pdata <- sub_metadata(pdata)	
				sub_meta <- data.frame(rowname=rownames(sub_pdata), sub_pdata, stringsAsFactors=FALSE, check.names=F)
				
				genes <- f[complete.cases(f), 1:3, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
				ct <- data.frame(rowname=rownames(countsTable), gene_symbol=genes[,2], countsTable, stringsAsFactors=F)

				collapse_data <- list(counts=ct, pdata=meta, sub_pdata=sub_meta, fdata=genes)
				return(collapse_data)
			} else {
				countsTable <- round(x)
				colnames(countsTable) <- id
				genes <- f[complete.cases(f), 1:3, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
				#rownames(countsTable) <- paste0(genes[,1], " : ", genes[,2])
				ct <- data.frame(rowname=rownames(countsTable), gene_symbol=genes[,2], countsTable, stringsAsFactors=F)

				rownames(pdata) <- id
				meta <- data.frame(rowname=rownames(pdata), pdata, stringsAsFactors=FALSE, check.names=F)
				
				sub_pdata <- sub_metadata(pdata)	
				sub_meta <- data.frame(rowname=rownames(sub_pdata), sub_pdata, stringsAsFactors=FALSE, check.names=F)		
				collapse_data <- list(counts=ct, pdata=meta, sub_pdata=sub_meta, fdata=genes)
				return(collapse_data)
			}
		} else {
			countsTable <- round(x)
			colnames(countsTable) <- id
			genes <- f[complete.cases(f), 1:3, drop=F]
			countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
			ct <- data.frame(rowname=rownames(countsTable), gene_symbol=genes[,2], countsTable, stringsAsFactors=F)

			rownames(pdata) <- id
			meta <- data.frame(rowname=rownames(pdata), pdata, stringsAsFactors=FALSE, check.names=F)
			
			sub_pdata <- sub_metadata(pdata)	
			sub_meta <- data.frame(rowname=rownames(sub_pdata), sub_pdata, stringsAsFactors=FALSE, check.names=F)	
			dat <- list(counts=ct, pdata=meta, sub_pdata=sub_meta, fdata=genes)
			return(dat)
		}
	})
	
	collapse_data2 <- reactive({
		load(paste0("data/",toupper(input$geo_acc),"/eset.RData"))
		pdata <- pData(eset)
		x <- exprs(eset)
		f <- fData(eset)
		rownames(x) <- f[,1]
		id <- unlist(lapply(rownames(pdata), function(x) unlist(strsplit(x, "_"))[2]))
		
		if("geo_accession" %in% names(pdata)) {
			if(length(unique(id))!= nrow(pdata)){
				pdata$ids <- id
				groupby <- pdata$ids
				if (!is.factor(groupby)) 
					groupby <- factor(groupby)
				groupby <- droplevels(groupby)
				sp <- split(seq(along = groupby), groupby)
				countsTable <- sapply(sp, function(i) rowSums(round(x)[,i, drop = FALSE]))
				colsToKeep <- sapply(sp, `[`, 1)
				pdata <- pdata[colsToKeep,]
				rownames(pdata) <- pdata$ids
				pdata <- pdata[,-which(colnames(pdata)=="ids")]				
				genes <- f[complete.cases(f), 1:2, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
				rownames(countsTable) <- paste0(genes[,1], " : ", genes[,2])
				#rownames(pdata) <- colnames(countsTable)
				sub_pdata <- sub_metadata(pdata)
				collapse_data <- list(counts=data.matrix(countsTable), pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
				return(collapse_data)
			} else {
				countsTable <- round(x)
				colnames(countsTable) <- id
				genes <- f[complete.cases(f), 1:2, drop=F]
				countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
				rownames(countsTable) <- paste0(genes[,1], " : ", genes[,2])
				rownames(pdata) <- id
				sub_pdata <- sub_metadata(pdata)
				collapse_data <- list(counts=data.matrix(countsTable), pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
				return(collapse_data)
			}
		} else {
			countsTable <- round(x)
			colnames(countsTable) <- id
			genes <- f[complete.cases(f), 1:2, drop=F]
			countsTable <- countsTable[which(rownames(countsTable) %in% genes[,1]),]
			rownames(countsTable) <- paste0(genes[,1], " : ", genes[,2])
			rownames(pdata) <- colnames(countsTable)
			sub_pdata <- sub_metadata(pdata)
			collapse_data <- list(counts=data.matrix(countsTable),pdata=pdata, sub_pdata=sub_pdata, fdata=genes)
			return(collapse_data)
		}
	})


	######################## Description ########################
	
	des_data <- reactive ({
		meta <- collapse_data()$pdata
		geo_rows <- length(unique(meta[,1])) #length(unique(unlist(lapply(meta[,1], function(x) unlist(strsplit(x,"_"))[2]))))
		srr_rows <- nrow(meta)
		datatable_to_use <- as.data.frame(read_feather(paste0("data/datatable_to_use.csv")),stringsAsFactors=F)
		return(list(geo_rows=geo_rows, srr_rows=srr_rows, datatable_to_use=datatable_to_use))
	})
	
	output$geo_summary <- DT::renderDataTable({
		study_link=paste0('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',toupper(input$geo_acc), '" onclick="ga(\'send\', \'event\', \'click\', \'link\', \'geoinfo\', 1)" target="_blank">', toupper(input$geo_acc), '</a>')
		dat <- data.frame(v1=study_link, v2=des_data()$geo_rows, v3=des_data()$srr_rows, 
			v4=des_data()$datatable_to_use[which(des_data()$datatable_to_use$`GEO accession`==toupper(input$geo_acc)),3],
			v5=des_data()$datatable_to_use[which(des_data()$datatable_to_use$`GEO accession`==toupper(input$geo_acc)),4],
			v6=des_data()$datatable_to_use[which(des_data()$datatable_to_use$`GEO accession`==toupper(input$geo_acc)),5])
		colnames(dat) <- c("Study link","No. of GEO samples","No. of SRA runs","Species","Title","Summary")
		
		datatable(t(dat), selection = 'none', escape = FALSE, filter = 'none', colnames = ""
			,options = list(dom = 't', columnDefs = list(list(width = '17%', targets = 0)))
		) %>% formatStyle(0,  fontWeight = 'bold')
	})
	
	######################## metadata column ########################
	
	meta_dat <- reactive({
		meta <- collapse_data()$pdata
		rownames(meta) <- meta[,1]
		meta <- meta[,-1,drop=F]
		sub_meta <- collapse_data()$sub_pdata
		rownames(sub_meta) <- sub_meta[,1]
		sub_meta <- sub_meta[,-1,drop=F]
		met <- list(meta=meta, sub_meta=sub_meta)
		return(met)

	})
	
	sub_meta_dat <- reactive({
		sub_meta <- collapse_data()$sub_pdata
		rownames(sub_meta) <- sub_meta[,1]
		sub_meta <- sub_meta[,-1,drop=F]
		return(sub_meta)
	})
	
	output$metadata_column <- renderUI({	
		col_names <- names(sub_meta_dat())
		checkboxGroupInput("property", h4("Select properties"), choices  = col_names, selected=col_names[1:length(col_names)])
	})
	
	output$full_metacol <- renderUI({
		validate(
			need(input$tab2=="metadata","")
		)
		choices <- if(nrow(sam_ontology())==0){
				c("Filtered metadata","Full metadata")
			} else {
				c("Filtered metadata","Full metadata","Ontology data")
			}
		radioButtons("full_meta", label = h4("Types of metadata:"), choices = choices, selected="Filtered metadata")
	})
		
	output$metadata_full <- DT::renderDataTable({
		datatable(meta_dat()$meta, selection = 'none',options = list(searchHighlight = TRUE, pageLength = if(nrow(meta_dat()$meta)>=10) {10} else {nrow(meta_dat()$meta)},
					scrollX = TRUE), rownames= TRUE
		)
	})
	output$metadata <- DT::renderDataTable({
		sub_meta_dat <- sub_meta_dat()
		datatable(sub_meta_dat[, colnames(sub_meta_dat) %in% input$property, drop = FALSE],
		selection = 'none',
		options = list(searchHighlight = TRUE, pageLength = if(nrow(sub_meta_dat)>=10) {10} else {nrow(sub_meta_dat)}, scrollX = TRUE),
			rownames= TRUE
		)
	})
	
	sam_ontology <- reactive ({
		sample_data <- all_samples[which(all_samples[,"gse_acc"] %in% input$geo_acc), ]
		ont_data <- samples_ont[which(samples_ont[,"Sample"] %in% sample_data[,"Sample"]),]
		x <- merge(sample_data, ont_data, by.x=3, by.y=1)
		x <- x[,c(2,1,3,5,6)]
		dat <- x[order(x[,1]),]
		return(dat)
	})
	output$ontology <- DT::renderDataTable({
		
		dat <- sam_ontology()
		datatable(dat,
		selection = 'none',
		options = list(searchHighlight = TRUE, pageLength = if(nrow(dat)>=10) {10} else {nrow(dat)}, scrollX = TRUE),
			rownames= FALSE
		)
	})

	output$downloadmeta <- downloadHandler(
		filename = function() { 
			if(input$full_meta=="Full metadata"){
				paste(toupper(input$geo_acc), '_full_metadata.csv', sep='') 
			} else if(input$full_meta=="Filtered metadata") {
				paste(toupper(input$geo_acc), '_filtered_metadata.csv', sep='') 
			} else {
				paste(toupper(input$geo_acc), 'sample_ontology_metadata.csv', sep='') 
			}
		},
		content = function(file) {
			if(input$full_meta=="Full metadata"){
				x <- meta_dat()$meta
				fwrite(x, row.names=TRUE, file)
			} else if(input$full_meta=="Filtered metadata"){
				x <- meta_dat()$sub_meta
				fwrite(x, row.names=TRUE, file)
			} else {
				x <- sam_ontology()
				fwrite(x, row.names=TRUE, file)
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
			shinyjs::show("downloadmeta")
		} else {
			shinyjs::hide("property")
			shinyjs::hide("full_meta")
			shinyjs::hide("downloadmeta")
		}
	})
	observeEvent(input$full_meta,{
		if (input$full_meta=="Full metadata") {
			shinyjs::hide("property")
			shinyjs::show("metadata_full")
			shinyjs::hide("metadata")
			shinyjs::hide("ontology")
		} else if(input$full_meta=="Filtered metadata"){
			shinyjs::show("property")
			shinyjs::hide("metadata_full")
			shinyjs::show("metadata")
			shinyjs::hide("ontology")
		} else {
			shinyjs::hide("property")
			shinyjs::hide("metadata_full")
			shinyjs::hide("metadata")
			shinyjs::show("ontology")			
		}
	})

	
	######################## Counts table ########################
	
	counts_dat <- reactive ({
		ct <- collapse_data()$counts
		rownames(ct) <- ct[,1]
		ct <- ct[,-1]
		return(ct)		
	})
	metac <- reactive({
		pdata <- collapse_data()$pdata
		return(nrow(pdata))
	})
	output$counts_sample_ui <- renderUI({
		sam <- nrow(sub_meta_dat())
		value <- if(sam>500) {
			500
		} else {
			sam
		}
		numericInput("counts_sample", label=h4('Number of samples to show') ,min = 2, max =sam, value = value, step= 1)
	})

	vcounts <- reactiveValues(doPlot = FALSE)
	observeEvent(input$load_counts, {
		vcounts$doPlot <- input$load_counts
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vcounts$doPlot <- FALSE
		}
	}) 
  
	output$counts_warn <- renderText({
		validate(
			need(input$counts_sample>1000,"")
		)
		paste0("You have selected ",input$counts_sample," samples to show which might take longer to load.")
	})

	output$counts <- DT::renderDataTable({
		if (input$geo_acc!="") {
			if(vcounts$doPlot==FALSE) return()
			isolate ({
				datatable(
					counts_dat()[,1:(1+input$counts_sample),drop=F] ,selection = 'none',options = list(pageLength = 15, scrollX = TRUE)
				)
			})
		}
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
				fwrite(counts_dat(), row.names=TRUE, file)
			} else {
				tr <- read.csv(file=paste0("data/",toupper(input$geo_acc),"/transcripts_counts.csv"), header=TRUE, row.names=1)
				fwrite(tr, row.names=TRUE, file)
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
	
	addResourcePath("data", "data")
	output$multiqc <- renderUI({
		tags$iframe(
			seamless="seamless",
			src=paste0("data/", toupper(input$geo_acc),"/multiqc_report.html"), height=1200, width='100%'
		)
	})
	output$downloadqc <- downloadHandler(
		filename = function() { 
			paste(toupper(input$geo_acc), '_QCreport.html', sep='') 
		},
		content = function(file) {
			file.copy(paste0("data/",toupper(input$geo_acc),"/multiqc_report.html"), file,overwrite = TRUE)
			file.remove("multiqc_report.html")
		},contentType = "text/html"
	)
		
	############################### Visualization ###############################
	
	############ Sample Correlation ############
	
	output$cor_sample_ui <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(any(c(input$preplots=="corr",input$preplots=="dens")),"")
		)
		sam <- nrow(sub_meta_dat())
		if(input$preplots=="dens"){
			value <- if(sam>100) {
				100
			} else {
				sam
			}
		} else {
			value <- if(sam>300) {
				300
			} else {
				sam
			}			
		}
		numericInput("cor_sample", label=h4('Number of samples') ,min = 2, max =sam, value = value, step= 1)
	})
	output$draw_cor_ui <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="corr","")
		)
		actionButton("draw_cor", "Draw correlation plot",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
	})	
	output$corr_warn <- renderText({
		validate(
			need(metac()>400,""),
			need(input$draw_cor==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the correlation plot.")
	})
	observeEvent(input$draw_cor,{
		shinyjs::hide("corr_warn")
	})

	vcorr <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_cor, {
		vcorr$doPlot <- input$draw_cor
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vcorr$doPlot <- FALSE
		}
	}) 

	output$corplot <- renderPlotly({
		if (vcorr$doPlot == FALSE) return()
		isolate({
			countsTable <- collapse_data()$counts
			rownames(countsTable) <- countsTable[,1]
			countsTable <- countsTable[,-c(1,2)]
			corr= cor(data.matrix(countsTable[,which(colSums(countsTable)>0)]),method = "spearman")
			rownames(corr) <- colnames(corr)
			corr <- corr[(1:input$cor_sample),(1:input$cor_sample)]
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
	})	

	############ Density plot ############
		
	output$draw_dens_ui <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="dens","")
		)
		actionButton("draw_dens", "Draw density plot",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
	})	
	
	vdens <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_dens, {
		vdens$doPlot <- input$draw_dens
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vdens$doPlot <- FALSE
		}
	}) 

	output$dens_warn <- renderText({
		validate(
			need(metac()>200,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the density plot.")
	})
	observeEvent(input$draw_dens,{
		shinyjs::hide("dens_warn")
	})
	
	observeEvent(input$geo_acc,{
		output$densityp <- renderPlotly({
			if (vdens$doPlot == FALSE) return()
			isolate({
				countsTable <- collapse_data()$counts
				rownames(countsTable) <- countsTable[,1]
				countsTable <- countsTable[,-c(1,2)]
				dat <- data.matrix(countsTable[,which(colSums(countsTable)>0)])
				pseudoCount = cpm(dat, normalized.lib.sizes=FALSE, log=TRUE)
				df = reshape2::melt(pseudoCount)
				df2 = data.frame(df, Samples = df[,2])			
				df2 <- df2[1:(input$cor_sample*(nrow(df2)/length(unique(df2$Samples)))),]				  
				p <- ggplot(df2, aes(x = value, colour = Samples)) + ylim(c(0, 0.25)) +
				geom_density()+  theme(legend.position = "right",legend.text=element_text(size=6), legend.key.size=unit(0.3,"cm"))+
				xlab("log2(counts per million)")
				ggplotly(p)	
			})
		})	
	})

	############# Heatmap ############

	heatmap_data <- reactive({
		ct <- collapse_data()$counts
		rownames(ct) <- ct[,1]
		heat_dat <- ct
		rownames(heat_dat) <- paste0(ct[,1], " : ", ct[,2])
		heat_dat <- as.matrix(heat_dat[,-c(1,2)])
		dge <- DGEList(counts=heat_dat[,which(colSums(heat_dat)>0)])
		dge_cpm <- cpm(dge, log=TRUE)		#data_mat <- data_mat[complete.cases(data_mat), ]
		exps= as.matrix(dge_cpm) - rowMeans(as.matrix(dge_cpm), na.rm=T)
		medAbsDev <-apply(exps,1,function(x) median(abs(x)))
		topGenes <- order(medAbsDev,decreasing=T)
		topGenes= topGenes[!is.na(topGenes)]
		expr_data<- data.matrix(exps[topGenes,,drop=F])
		ht <- data.frame(rowname=rownames(expr_data), expr_data, stringsAsFactors=FALSE)
		rownames(ht) <- ht[,1]
		ht <- ht[,-1]
		dat <- data.matrix(ht)
		sub_meta <- sub_meta_dat()
		return(list(counts=dat, sub_pdata=sub_meta))
	})
	data_dl <- reactive({
		expr_data <- heatmap_data()$counts
		expr_data_in <- expr_data[1:input$cluster_ngenes,]
		gene_sym <- unlist(lapply(rownames(expr_data_in), function(x) unlist(strsplit(x, ' : '))[2]))
		dat2 <- data.frame(GeneSymbol=gene_sym, expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	})

	output$heat_prop <- renderUI({	
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		col_names <- names(sub_meta_dat())
		checkboxGroupInput("heat_property", h4("Select properties"), choices  = col_names 
		,selected= if(length(col_names)>2) {
					col_names[(length(col_names)-1):length(col_names)]
				} else {
					col_names[1:length(col_names)]
				}
		)
	})
	output$cluster_genes <- renderUI({
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		numericInput("cluster_ngenes", label=h4('Number of highly variable genes'), 
			value= 100, min = 1, max = 20000, width="230px"
		)
		
	})
	output$cluster_meth <- renderUI({	
		validate(
			need(input$tab2=="plots",""),
			need(input$preplots=="heatmap","")
		)
		selectInput("clustering_method", label = h4("Grouping genes and samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected="Pearson correlation"
		)
	})
	
	vheat <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heat, {
		vheat$doPlot <- input$draw_heat
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheat$doPlot <- FALSE
		}
	}) 

	vheatscr <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatscr, {
		vheatscr$doPlot <- input$draw_heatscr
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatscr$doPlot <- FALSE
		}
	}) 

	vheatint <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatint, {
		vheatint$doPlot <- input$draw_heatint
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatint$doPlot <- FALSE
		}
	}) 

	output$heat_warn <- renderText({
		validate(
			need(metac()>700,""),
			need(input$draw_heat==0,""),
			need(input$draw_heatscr==0,""),
			need(input$draw_heatint==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the heatmap.")
	})
	observeEvent(input$draw_heat,{
		shinyjs::hide("heat_warn")
	})
	observeEvent(input$draw_heatscr,{
		shinyjs::hide("heat_warn")
	})
	observeEvent(input$draw_heatint,{
		shinyjs::hide("heat_warn")
	})

	### Static heatmap

	output$HeatStatic = renderPlot({
		if (vheat$doPlot == FALSE ) return()
		isolate({
			dat <- heatmap_data()
			complex_heatmap (data_mat=dat$counts, metadata=dat$sub_pdata, property=input$heat_property,
				n_genes=input$cluster_ngenes, clustering_method=input$clustering_method, type="Fit in Screen")
		})
	}
		,height= eventReactive (input$draw_heat,{
			m <- sub_meta_dat()
			df <- m[,which(colnames(m) %in% input$heat_property),drop=FALSE]
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
		})
	)
	
	output$HeatStaticscr = renderPlot({
		if (vheatscr$doPlot == FALSE ) return()
		isolate({
			dat <- heatmap_data()
			complex_heatmap (data_mat=dat$counts, metadata=dat$sub_pdata, property=input$heat_property, 
				input$cluster_ngenes, clustering_method=input$clustering_method, type="Scrollable")			
		})
	}
		,height= eventReactive (input$draw_heatscr,{
			if(input$cluster_ngenes>=100) {
				value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24))
			} else {
				value <- (length(input$heat_property)*25) + 400+(input$cluster_ngenes*20)
			}
			return(value)
		})
	)
	
	output$download <- downloadHandler(
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
				png(file, width = 15, height = 12, units = "in",res=300)
				dat <- heatmap_data()
				complex_heatmap (data_mat=dat$counts, metadata=dat$sub_pdata, property=input$heat_property,
					n_genes=input$cluster_ngenes, clustering_method=input$clustering_method, type="Fit in Screen")
				dev.off()
			} else if(input$download_type_fit == "csv"){
				dat2 <- data_dl()
				write.csv(dat2, file)
			} else {
				return()
			}
		}
	)
	
	output$downloadscr <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scr=="png"){
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
			if(input$download_type_scr == "png") {
				height = reactive({
					if(input$cluster_ngenes>=1000) {
						value <- max((length(input$heat_property)*25)+(4*96),(length(input$heat_property)*25)+(input$cluster_ngenes*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property)*25)+(10*96),(length(input$heat_property)*25)+(input$cluster_ngenes*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				dat <- heatmap_data()
					complex_heatmap (data_mat=dat$counts, metadata=dat$sub_pdata, property=input$heat_property,
					n_genes=input$cluster_ngenes, clustering_method=input$clustering_method, type="Scrollable")
				dev.off()
			} else{
				dat2 <- data_dl()
				write.csv(dat2, file)
			}			
		}
	)

	### Interactive
	
	height_in = eventReactive(input$draw_heatint,{
		if(!is.null(input$heat_property)) {
			df <- sub_meta_dat()[,input$heat_property,drop=FALSE]
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

	output$iheatmapui = renderUI({
		if (vheatint$doPlot == FALSE) return()
		isolate({
			iheatmaprOutput("iheatmap", width="auto", height = height_in()) %>% withSpinner(type=5)
		})
	})
	
	output$iheatmap = renderIheatmap({
		if (vheatint$doPlot == FALSE) return()
		isolate({
			dat <- heatmap_data()
			print(interactive_heatmap (data_mat=dat$counts, metadata=dat$sub_pdata, property=input$heat_property, 
				n_genes=input$cluster_ngenes, clustering_method=input$clustering_method, signature=FALSE))			
		})
	})

	######
	
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="static" && input$fir_or_scr=="fitinscreen"){
			if (vheat$doPlot == FALSE ) return()
			isolate({
				shinyjs::show("download")
				shinyjs::show("download_type_fit")
			})
		} else {
			shinyjs::hide("download")
			shinyjs::hide("download_type_fit")
		}
	})
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="static" && input$fir_or_scr=="scrollable"){
			shinyjs::show("downloadscr")
			shinyjs::show("download_type_scr")
		} else {
			shinyjs::hide("downloadscr")
			shinyjs::hide("download_type_scr")
		}
	})
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="static" && input$fir_or_scr=="fitinscreen"){
			shinyjs::show("draw_heat")
		} else {
			shinyjs::hide("draw_heat")
		}
	})
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="static" && input$fir_or_scr=="scrollable"){
			shinyjs::show("draw_heatscr")
		} else {
			shinyjs::hide("draw_heatscr")
		}
	})
	observe({
		if (input$tab2=="plots" && input$preplots=="heatmap" && input$heatvis=="inter"){
			shinyjs::show("draw_heatint")
		} else {
			shinyjs::hide("draw_heatint")
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

	############# PCA ############

	load("data/allColors.rda")
	
	output$pca_prop <- renderUI({
		coln <- colnames(sub_meta_dat())
		selectInput(inputId="pca_property",label=h4('Select a property'),
			choices=coln, selected=coln[length(coln)]
		)
	})
	vpca <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_pca, {
		vpca$doPlot <- input$draw_pca
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vpca$doPlot <- FALSE
		}
	}) 

	vpca3 <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_pca3, {
		vpca3$doPlot <- input$draw_pca3
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vpca3$doPlot <- FALSE
		}
	}) 

	output$pca_warn <- renderText({
		validate(
			need(metac()>300,""),
			need(input$draw_pca==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the PCA plot.")
	})
	observeEvent(input$draw_pca,{
		shinyjs::hide("pca_warn")
	})
	output$pca3_warn <- renderText({
		validate(
			need(metac()>300,""),
			need(input$draw_pca3==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the PCA plot.")
	})
	observeEvent(input$draw_pca3,{
		shinyjs::hide("pca3_warn")
	})

	output$legendplot = renderPlot({
		if (vpca$doPlot == FALSE ) return()
		isolate({
			m <- meta_dat()$sub_meta
			allColors = unique(allColors)[-c(6,7,11,13,15,16,19,20,22,24,25,28,35)]
			groups <- m[,which(colnames(m)==input$pca_property)]		
			groups[is.na(groups)] <- "NA"
			groups <- relevel(groups, unique(as.character(groups)[1]))
			if (length(unique(groups))>22) {
				colcols <- colorRampPalette(brewer.pal(9,"Set1"))(length(levels(groups)))
			} else {
				colcols <- allColors[seq_along(levels(groups))]
			}

			if (length(unique(groups))>=50) {
				ncolm=2
			} else {
				ncolm=1
			} 

			plot(1, type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i")
			legend("top",title = input$pca_property
			,bty = "n",cex = 1, ncol=ncolm,	fill=colcols, 
			legend= sapply(levels(groups), function(x) paste(strsplit(x, "(?<=.{20})", perl = TRUE)[[1]], collapse="\n")))
		})
	}	
		,height = eventReactive(input$draw_pca,{
			m <- meta_dat()$sub_meta
			groups <- m[,which(colnames(m)==input$pca_property)]
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
		})
	)
	output$legendplot3 = renderPlot({
		if (vpca3$doPlot == FALSE ) return()
		isolate({
			m <- meta_dat()$sub_meta
			allColors = unique(allColors)[-c(6,7,11,13,15,16,19,20,22,24,25,28,35)]
			groups <- m[,which(colnames(m)==input$pca_property)]		
			groups[is.na(groups)] <- "NA"
			groups <- relevel(groups, unique(as.character(groups)[1]))
			if (length(unique(groups))>22) {
				colcols <- colorRampPalette(brewer.pal(9,"Set1"))(length(levels(groups)))
			} else {
				colcols <- allColors[seq_along(levels(groups))]
			}

			if (length(unique(groups))>=50) {
				ncolm=2
			} else {
				ncolm=1
			} 

			plot(1, type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i")
			legend("top",title = input$pca_property
			,bty = "n",cex = 1, ncol=ncolm,	fill=colcols, 
			legend= sapply(levels(groups), function(x) paste(strsplit(x, "(?<=.{20})", perl = TRUE)[[1]], collapse="\n")))
		})
	}	
		,height = eventReactive(input$draw_pca3,{
			m <- meta_dat()$sub_meta
			groups <- m[,which(colnames(m)==input$pca_property)]
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
		})
	)
	
	valuespca <- reactive ({
		countsTable <- collapse_data()$counts
		rownames(countsTable) <- countsTable[,1]
		countsTable <- countsTable[,-c(1,2)]
		data_mat <- cpm(countsTable[,which(colSums(countsTable)>0)], normalized.lib.sizes=FALSE, log=TRUE)
		exps_pca <- data_mat[complete.cases(data_mat),]
		pseudoCount2 <- exps_pca[abs(rowSums(exps_pca))>0,]
		pseudoCount <- pseudoCount2[apply(pseudoCount2, 1, var, na.rm=TRUE) != 0, ]
		b.PCA <- prcomp(t(pseudoCount),retx=TRUE,center=TRUE, scale=TRUE)
		dat <- data.frame(rowname=rownames(b.PCA$x), b.PCA$x, stringsAsFactors=FALSE)
		rownames(dat) <- dat[,1]
		dat <- data.matrix(dat[,-1])
		return(dat)
	})

	output$p2D_ui <- renderUI({
		pairsD3Output("p2D",width = "auto",height=800) %>% withSpinner(type=5)
	})
	output$pca3Dplot_ui <- renderUI({
		rglwidgetOutput("pca3Dplot", width = "auto", height = "700px") %>% withSpinner(type=5)
	})

	output$p2D <- renderPairsD3({
		if (vpca$doPlot == FALSE) return()
		isolate({
			print(pca(valuespca(), sub_meta_dat(), input$pca_property, "pca2D"))
		})
	})
	output$pca3Dplot <- renderRglwidget({
		if (vpca3$doPlot == FALSE) return()
		isolate({
			try(rgl.close())
			pca(valuespca(), sub_meta_dat(), input$pca_property, "pca3D")
			rglwidget()
		})
	})
	output$downloadPCA2D = downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			return(paste(datasetname,"_",input$pca_property,"_pca_2D.html",sep=""))
		},
		content = function(file) {
			savePairs(pca(valuespca(), sub_meta_dat(), input$pca_property, "pca2D"), file)
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
	
	#### t-SNE
	
	vtsne3 <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_tsne3, {
		vtsne3$doPlot <- input$draw_tsne3
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vtsne3$doPlot <- FALSE
		}
	}) 

	vtsne2 <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_tsne2, {
		vtsne2$doPlot <- input$draw_tsne2
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vtsne2$doPlot <- FALSE
		}
	}) 

	values_tsne <- reactive ({
		countsTable <- counts_dat()[,-1]
		exps_pca <- countsTable[complete.cases(countsTable),]
		pseudoCount2 <- exps_pca[abs(rowSums(exps_pca))>0,]
		pseudoCount <- pseudoCount2[apply(pseudoCount2, 1, var, na.rm=TRUE) != 0, ]
		x <- dist(t(as.matrix(pseudoCount)))	
		perplexity <- if (ncol(pseudoCount)>100){
							30
						} else if(ncol(pseudoCount)>3){
							(ncol(pseudoCount)/3)-1
						} else {
							0
						}
		tsne_dat <- Rtsne(x, is_distance=TRUE, perplexity=perplexity, dims=3)
		dat <- tsne_dat$Y
		return(dat)
	})

	output$tsne3Dplot_ui <- renderUI({
		plotlyOutput("tsne3D",width=900,height=600) %>% withSpinner(type=5)
	})
	output$tsne3D <- renderPlotly({
		if (vtsne3$doPlot == FALSE) return()
		isolate({
			print(pca(values_tsne(), sub_meta_dat(), input$pca_property, "tsne3D"))
		})
	})

	output$t2Dplot_ui <- renderUI({
		plotlyOutput("tsne2D",width=900,height=600) %>% withSpinner(type=5)
	})

	output$tsne2D <- renderPlotly({
		if (vtsne2$doPlot == FALSE) return()
		isolate({
			pca(values_tsne(), sub_meta_dat(), input$pca_property, "tsne2D")
		})
	})

	################################################ Analysis nav page ################################################
	
	observe({
		shinyjs::hide(selector = "#grin li a[data-value=analyze]")
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

	observeEvent(input$analysisbtn,{
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
		dat <- collapse_data2()
		group <- dat$sub_pdata[,1]
		if(any(as.numeric(table(group))<2)){
			y <- DGEList(counts=dat$counts, genes= dat$fdata[,2])
			keep <- rowSums(cpm(y)>1) >= min(table(group))
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
			y <- DGEList(counts=dat$counts, genes= dat$fdata[,2], group=group)
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
		numericInput("sig_lev", label=h4('Level of significance (alpha)') ,min = 0, max = 0.1, value = 0.05, step= 0.01)
	})
	output$effect <- renderUI({
		sliderInput("fc", label=h4('Effect size (fold change)') ,min = 1, max = 5, value = 2, step= 0.1)
	})
	
	output$samples <- renderUI({
		val <- nrow(sub_meta_dat())
		numericInput("n_samples", label=h4('No. of samples') ,min = 2, max= val+5, step= 2, value = val)
	})

	output$depth_ui <- renderUI({
		dat <- data.matrix(counts_dat()[,-1])
		val <- round(unique(pow_dat()$depth),2)
		numericInput("depth", label=h4('Sequencing depth (in million)') ,min = 2, max= val+5, step= 0.5, value = val, width="230px")
	})
	observe({
		if(input$sigtab=="power_curve"){
			shinyjs::show("depth")
		} else {
			shinyjs::hide("depth")
		}
	})
	
	output$gen_pow_btn_ui <- renderUI({
		validate(
			need(input$sigtab=="power_curve","")
		)
		actionButton("gen_pow_btn", "Draw power curve",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
	})	
	vpow <- reactiveValues(doPlot = FALSE)
	observeEvent(input$gen_pow_btn, {
		vpow$doPlot <- input$gen_pow_btn
	})
	observeEvent(input$analysisbtn, {
		vpow$doPlot <- FALSE
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vpow$doPlot <- FALSE
		}
	}) 
	
	output$power1_warn <- renderText({
		validate(
			need(metac()>200,""),
			need(input$gen_pow_btn==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the Power curve.")
	})
	observeEvent(input$gen_pow_btn,{
		shinyjs::hide("power1_warn")
	})

	df.power <- eventReactive(vpow$doPlot,{
		if (vpow$doPlot == FALSE) return()
		isolate({
			power <- rnapower(n=c(1:input$n_samples), depth=input$depth, cv= unique(pow_dat()$cv), alpha=input$sig_lev, effect=input$fc)
			df.power <- as.data.frame(cbind(samples=1:input$n_samples, power))
			return(df.power)
		})
	})
	
	output$powerplot <- renderPlotly({	
		if (vpow$doPlot == FALSE) return()
		isolate({
			if(input$n_samples>10) {
				p.line <- ggplot(df.power()) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ xlim(1,input$n_samples) +labs(x="No. of samples", y="Power") #+scale_x_discrete(limit = c(1:input$n_samples))
			} else {
				p.line <- ggplot(df.power()) + geom_line(size=1, color="blue", aes(x=samples, y=power))+ labs(x="No. of samples", y="Power") +scale_x_discrete(limit = c(1:input$n_samples))		
			}
			ggplotly(p.line)
		})
	})	
	
	observeEvent(input$signature,{
		output$pow_method <- renderUI({
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
			need(input$sigtab=="detec_gene","")
		)
		val = 0.8
		numericInput("power_level", label=h4('Power') ,min = 0.1, max= 1, step= 0.1, value = val)
	})

	output$gen_detect_btn_ui <- renderUI({
		validate(
			need(input$sigtab=="detec_gene","")
		)
		actionButton("gen_detect_btn", "Generate plot",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
	})	

	v2 <- reactiveValues(doPlot2 = FALSE)
	observeEvent(input$gen_detect_btn, {
		v2$doPlot2 <- input$gen_detect_btn
	})
	observeEvent(input$analysisbtn, {
		v2$doPlot2 <- FALSE
	})  
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			v2$doPlot2 <- FALSE
		}
	}) 

	output$power2_warn <- renderText({
		validate(
			need(metac()>200,""),
			need(input$gen_detect_btn==0,"")
		)
		paste0("This dataset has ",metac()," sample runs and might take longer to create the detectability of genes plot.")
	})
	observeEvent(input$gen_detect_btn,{
		shinyjs::hide("power2_warn")
	})

	genes <- reactive ({
		f <- collapse_data()$fdata
		genes <- f[complete.cases(f), 1:2, drop=F]
		return(genes)
	})

	output$gene_warn <- renderText({
		validate(
			need(input$search_gene!="","")
		)
		g <- unlist(strsplit(input$search_gene,","))
		gvalid <- g[which(toupper(g) %in% toupper(genes()[,2]))]
		if(length(gvalid)>=1){
			gnull <- gvalid[which(!toupper(gvalid) %in% toupper(to_vary()$genes))]
			if(length(gnull)>1){
				paste("Gene:",paste0(gnull, collapse=","),"are filtered out because of low counts")
			} else if(length(gnull)==1){
				paste("Gene:",paste0(gnull, collapse=","),"is filtered out because of low counts")
			} else {
				return()
			}
		} else {
			return()
		}
	})
	output$gene_sym_warn <- renderText({
		validate(
			need(input$search_gene!="","")
		)
		g <- unlist(strsplit(input$search_gene,","))
		if(any(!toupper(g) %in% toupper(genes()[,2]))){
			gnull <- g[which(!toupper(g) %in% toupper(genes()[,2]))]
			if(length(gnull)>1){
				paste("Genes:",paste0(gnull, collapse=","),"are invalid gene symbols")
			} else if(length(gnull)==1){
				paste("Gene:",paste0(gnull, collapse=","),"is an invalid gene symbol")
			} else {
				return()
			}
		} else {
			return()
		}
	})
	output$gene_sym_green <- renderText({
		validate(
			need(input$search_gene!="","")
		)
		g <- unlist(strsplit(input$search_gene,","))
		if(any(toupper(g) %in% toupper(to_vary()$genes))){
			gnull <- g[which(toupper(g) %in% toupper(to_vary()$genes))]
			if(length(gnull)>1){
				paste("Click 'Generate plot' button to search for genes:",paste0(gnull, collapse=","))
			} else if(length(gnull)==1){
				paste("Click 'Generate plot' button to search for gene:",paste0(gnull, collapse=","))
			} else {
				return()
			}
		} else {
			return()
		}
	})
	
	observe({
		if(input$sigtab=="detec_gene"){
			shinyjs::show("search_gene")
			shinyjs::show("gene_warn")
			shinyjs::show("gene_sym_warn")
		} else {
			shinyjs::hide("search_gene")
			shinyjs::hide("gene_warn")
			shinyjs::hide("gene_sym_warn")
		}
	})

	
	to_vary <- reactive ({
		logCPM_range <- if(min(pow_dat()$AveLogCPM)>0) {
							0:(max(pow_dat()$AveLogCPM)*100) /100
						} else {
							(min(pow_dat()$AveLogCPM)*100):(max(pow_dat()$AveLogCPM)*100) /100
						}
		depth_range <- 2^logCPM_range * unique(pow_dat()$depth)
		
		effect_vary <- NULL
		effect_vary <- c(effect_vary,input$fc)
		alpha_vary <- NULL
		alpha_vary <- c(alpha_vary,input$sig_lev)
		samples_vary <- NULL
		samples_vary <- c(samples_vary,input$n_samples)
		power_vary <- NULL
		power_vary <- c(power_vary,input$power_level)
		
		bcov_range <- as.matrix(rnapower(depth=depth_range, n=samples_vary, alpha=alpha_vary, power=power_vary, effect=effect_vary))
		colnames(bcov_range) <- "Line of detectability"
		
		df_logCPM_bcov_rnapower <- melt(data.frame(logCPM=logCPM_range, bcov=bcov_range), id="logCPM")
		df_logCPM_bcov <- data.frame(logCPM=pow_dat()$AveLogCPM, bcov=pow_dat()$bcov_tagwise, genes=pow_dat()$genes)
		df_logCPM_bcov_rnapower_overlap <- data.frame(logCPM=df_logCPM_bcov$logCPM, variable="Genes", value=df_logCPM_bcov$bcov, genes=df_logCPM_bcov$genes)
		df_logCPM_bcov_rnapower$variable <- gsub(c("bcov.Line.of.detectability|Line.of.detectability"), paste0("Line of detectability"), df_logCPM_bcov_rnapower$variable)
		df_logCPM_bcov_rnapower$genes <- paste0(rep(" ", nrow(df_logCPM_bcov_rnapower)))
		df_logCPM_bcov_rnapower_overlap2 <- rbind(df_logCPM_bcov_rnapower_overlap, df_logCPM_bcov_rnapower)
		return(df_logCPM_bcov_rnapower_overlap2)
	})
			
	output$detect_gene <- renderPlotly({
		if (v2$doPlot2 == FALSE) return()
		isolate({
			dataset <- to_vary()
			dat <- dataset[which(toupper(dataset$genes) %in% toupper(unlist(strsplit(input$search_gene,",")))),]

			p <- ggplot(dataset, aes(x=logCPM, y=value, colour=variable)) + 
				geom_point(size=0.7, alpha=0.6, aes(text=paste0("Gene symbol: ",dataset$genes))) +
				labs(x="Average log2(counts per million)", y="Biological coefficient of variation") + 
				guides(colour=guide_legend(override.aes = list(size=10, alpha=1), title.position="top", title="")) + 
				theme(legend.text=element_text(size = 11), axis.text=element_text(size=10), axis.title=element_text(size=12))
			gene_det <- p + geom_point(data=dat,color="darkgreen",size=1)
			ggplotly(gene_det, source="gene_det_source", tooltip=c("text","x","y"))
		})
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

	################### Signature ##################
		
	observeEvent(input$geo_acc2, {
		output$metadata_analysis <- renderUI ({
			col_names <- colnames(sub_meta_dat())
			selectInput(inputId="analysis_property",label=h4('Variable of interest'),
				choices=c('',col_names), selected=NULL
			)
		})
	})
	
	metadata_small <- reactive({
		metadata <- sub_meta_dat()
		metadata <- metadata[, input$analysis_property,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- input$analysis_property
		return(metadata)
	})

	output$ana_type <- renderUI({
		validate(
			need(input$analysis_property!='',"")
		)
		metadata <- sub_meta_dat()
		selectInput(inputId="analysis_type",label=h4('Type of comparison'),
			choices= 
				if(length(levels(metadata_small()[,1]))==2 && ncol(metadata)>1){
					c('Two group without covariate','Two group with covariate')
				} else if(length(levels(metadata_small()[,1]))==2 && ncol(metadata)==1 && any(as.numeric(table(metadata_small()[,1]))<2)){
					'Two group without covariate'
				} else if(length(levels(metadata_small()[,1]))>2 && ncol(metadata)==1 && all(as.numeric(table(metadata_small()[,1]))>=2)){
					c('Two group without covariate','Multi group')
				} else if(length(levels(metadata_small()[,1]))>2 && ncol(metadata)>1 && all(as.numeric(table(metadata_small()[,1]))>=2)){
					c('Two group without covariate','Two group with covariate', 'Multi group')
				} else {
					'Two group without covariate'
				},
			selected= 'Two group without covariate'
		)
	})
	output$group1 <- renderUI ({
		validate(
			need(any(c(input$analysis_type=="Two group without covariate",input$analysis_type=="Two group with covariate")),""),
			need(input$analysis_property!='',"")
		)		
		selectInput(inputId="input_group1", label=h4('Experimental group'),
		choices=c('',levels(metadata_small()[,1])), selected = NULL)
	})
	output$group2 <- renderUI ({
		validate(
			need(any(c(input$analysis_type=="Two group without covariate",input$analysis_type=="Two group with covariate")),""),
			need(input$analysis_property!='',"")
		)		
		available <- levels(metadata_small()[,1])[which(levels(metadata_small()[,1])!=input$input_group1)]
		selectInput( inputId="input_group2", label=h4('Control group'),
		choices=c('',available), selected = NULL)
	})	
	output$cov <- renderUI({
		validate(
			need(input$analysis_property!='',""),
			need(input$analysis_type=="Two group with covariate","")
		)
		metadata <- sub_meta_dat()
		selectInput(inputId="covar",label=h4('Select covariates'), 
			choices= colnames(metadata)[which(colnames(metadata)!=input$analysis_property)]
			, multiple = TRUE, selected= NULL
		)
	})
	output$gen_reg_sig_ui <- renderUI({
		validate(
			need(input$analysis_property!='',""),
			need(all(design_levels()>=2), "")
		)
		if(input$analysis_type=="Two group without covariate" && input$input_group2!=''){
			actionButton("gen_reg_sig", "Generate signature",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
		} else if(input$analysis_type=="Two group with covariate" & !is.null(input$covar)){
			if(is.null(design_rank())){
				actionButton("gen_reg_sig", "Generate signature",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
			} else {
				return()
			}
		} else if(input$analysis_type=="Multi group"){
			actionButton("gen_reg_sig", "Generate signature",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
		} else {
			return()
		}
	})	
	output$sig1_warn <- renderText({
		validate(
			need(nrow(design_data())>100,""),
			need(input$gen_reg_sig==0,""),
			need(input$input_group2!="","")
		)
		paste0("This dataset has ",nrow(design_data())," sample runs and might take longer to generate signature.")
	})
	observeEvent(input$gen_reg_sig,{
		shinyjs::hide("sig1_warn")
	})

	metadata_exp_cov <- reactive({
		metadata <- sub_meta_dat()
		metadata <- metadata[, c(input$analysis_property,input$covar),drop=FALSE]
		return(metadata)
	})
	design_data <- reactive({
		if(input$analysis_type=="Two group without covariate"){
			dat1 <- metadata_small()[metadata_small()[,1] %in% input$input_group1,,drop=FALSE]
			dat1$"Selected groups"[dat1[,1]==input$input_group1] <- "experimental"
			dat2 <- metadata_small()[metadata_small()[,1] %in% input$input_group2,,drop=FALSE]
			dat2$"Selected groups"[dat2[,1]==input$input_group2] <- "control"			
			dat <- rbind(dat1,dat2)
			dat <- dat[, c("Selected groups", input$analysis_property)]
			return(dat)
		} else if(input$analysis_type=="Two group with covariate"){
			dat1 <- metadata_exp_cov()[metadata_exp_cov()[,1] %in% input$input_group1,,drop=FALSE]
			dat1$"Selected groups"[dat1[,1]==input$input_group1] <- "experimental"
			dat2 <- metadata_exp_cov()[metadata_exp_cov()[,1] %in% input$input_group2,,drop=FALSE]
			dat2$"Selected groups"[dat2[,1]==input$input_group2] <- "control"
			dat <- rbind(dat1,dat2)
			dat <- dat[ ,c("Selected groups", input$analysis_property, input$covar)]
			return(dat)
		} else if(input$analysis_type=="Multi group"){
			dat <- metadata_small()
			return(dat)
		} else {
			return()
		}
	})

	design_levels <- reactive({
		if(input$analysis_type=="Two group without covariate"){
			metad <- design_data()
			design_data <- metad[,-1,drop=F]
			nlev <- as.numeric(apply(design_data, 2, function(x)length(unique(x))))
			return(nlev)
		} else if(input$analysis_type=="Two group with covariate"){
			metad <- design_data()
			design_data <- metad[,-1,drop=F]
			nlev <- as.numeric(apply(design_data, 2, function(x)length(unique(x))))
			return(nlev)
		} else if(input$analysis_type=="Multi group"){
			design_data <- sub_meta_dat()
			nlev <- as.numeric(apply(design_data, 2, function(x)length(unique(x))))
			return(nlev)			
		} else {
			return()
		}
	})

	design_rank <- reactive({
		if(all(design_levels()>=2)){
			if(input$analysis_type!="Multi group"){
				metad <- design_data()
				analysis_metadata <- metad[which(metad[,"Selected groups"]=="experimental" | metad[,"Selected groups"]=="control"),c(input$analysis_property,input$covar),drop=FALSE]
				group <- analysis_metadata[,input$analysis_property, drop=F]
				covar <- analysis_metadata[,input$covar, drop=F]
				design_data <- cbind(covar, group)
				design_data[] <- lapply(design_data, factor)
				design <- model.matrix(~., data=design_data)
				ne <- nonEstimable(design)	
				return(ne)
			} else {
				design <- model.matrix(~., data=design_data())
				ne <- nonEstimable(design)	
				return(ne)
			}
		} else {
			return()
		}
	})
	
	output$rank_warn <- renderText({
		validate(
			need(any(input$analysis_type=='Two group with covariate',input$analysis_type=='Multi group'),"")
			,need(all(design_levels()>=2), "")
			,need(!is.null(design_rank()), "")			
		)
		paste("Design matrix is not of full rank.  The following coefficients are not estimable:\n", paste(design_rank(), collapse = " ; "))
	})
	output$level_warn <- renderText({
		validate(
			need(any(design_levels()<2), "")
			,need(nrow(design_data())>0,"")
		)	
		paste("All the variables need to have 2 or more levels")
	})

	sigdata <- reactive({
		filepath <- paste('data/',toupper(input$geo_acc2),'/',toupper(input$geo_acc2),'_signatureData.txt', sep='')
		analysis_genes <- data.frame(Ensembl_ID=rownames(counts_dat()), Gene_symbol=counts_dat()[,1])
		metadata <- sub_meta_dat()
		if(input$analysis_type=='Two group without covariate' && input$input_group2!=''){
			signature_data= analysis_2g_wocov(counts=counts_dat()[,-1], metadata=metadata, genes=analysis_genes, property=input$analysis_property,
				group1=input$input_group1, group2=input$input_group2)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,"Log_FoldChange"], Significance_pvalue=signature_data$signaturesData_final[,"Adjusted_pvalue"], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$analysis_type=='Two group with covariate' && input$covar!=''){
			signature_data= analysis_2g_withcov(counts=counts_dat()[,-1], metadata=metadata, genes=analysis_genes, property=input$analysis_property,
				covariate=input$covar, group1=input$input_group1, group2=input$input_group2)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,"Log_FoldChange"], Significance_pvalue=signature_data$signaturesData_final[,"Adjusted_pvalue"], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$analysis_type=='Multi group') {
			signature_data= analysis_mg_wocov(counts=counts_dat()[,-1], metadata=metadata, genes=analysis_genes, property=input$analysis_property)
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,ncol(signature_data$signaturesData_final)-3], Significance_pvalue=signature_data$signaturesData_final[,ncol(signature_data$signaturesData_final)], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else {
			return()
		}
		return(signature_data)
	})
	
	output$sig_meta = DT::renderDataTable({
		validate(
			need(input$cr_reg_sig=="meta_reg","")
		)
		datatable( design_data()
			,selection = 'none'
			,extensions = 'Buttons'
			,options = list(dom = 'Bfrtip',buttons = 'csv', pageLength = if(nrow(design_data())>=15) {15} else {nrow(design_data())}, scrollX = TRUE)
			,rownames= TRUE
		)
	}, server = FALSE) 


	output$sig_data <- DT::renderDataTable({
		validate(
			need(input$gen_reg_sig,"")
		)
		datatable(
			sigdata()$signaturesData_final[,which(!colnames(sigdata()$signaturesData_final) %in% c("logCPM", "PValue"))]
			,selection = 'none'
			,rownames = FALSE
			,filter = 'top' 
			,options = list(pageLength = 15, scrollX = TRUE,columnDefs = list(list(className = 'dt-center', targets ="_all"))
			)
		)
	})	
	
	observeEvent(input$gen_reg_sig, {
		if (input$gen_reg_sig) {
			updateTabsetPanel(session, "cr_reg_sig", selected = "sig_reg")
		}
    })

	output$downloadsigdata <- downloadHandler(
		filename = function() { 
			paste(toupper(input$geo_acc2), '_signatureData.csv', sep='') 
		},
		content = function(file) {
			write.csv(sigdata()$signaturesData_final[,-1], file, row.names=FALSE )
		}
	)

	observe({
		shinyjs::hide("gen_reg_sig")
	})
	observeEvent(input$covar,{
		if(!is.null(input$covar) && !is.null(design_rank())){
			shinyjs::hide("gen_reg_sig")
		} else {
			shinyjs::show("gen_reg_sig")
		}
	})
	
	observe({
		hide(selector = "#cr_reg_sig li a[data-value=sig_reg]")
	})
	observeEvent(input$gen_reg_sig,{
		if (!is.null(input$gen_reg_sig)) {
			shinyjs::show(selector = "#cr_reg_sig li a[data-value=sig_reg]")
		} else {
			shinyjs::hide(selector = "#cr_reg_sig li a[data-value=sig_reg]")
		}
	})
	observeEvent(input$covar, {
		if(input$covar!="" & !is.null(design_rank())){
			shinyjs::hide("downloadsigdata")
		} else {
			shinyjs::show("downloadsigdata")
		}
	})
	observeEvent(input$cr_reg_sig,{
		if(input$cr_reg_sig=="sig_reg"){
			shinyjs::show("downloadsigdata")
		}	
	})
	observeEvent(input$signature,{
		output$sig_method <- renderUI({
			if (input$signature=='sig_data') {
				includeMarkdown("www/signature_method.Rmd")
			}
		})
	})
	
	
	##### show to ilincs
	
	output$dynamiclink <- renderUI({
		if(input$analysis_type=='Two group without covariate'){
			if (input$input_group2!='' && input$analysis_type=='Two group without covariate') {		
				actionButton(inputId='ul2ilincs', label=div("Upload signature to iLINCS",img(src="images/ilincs_logo.png", height = "20px", width="20px")), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			} 
		} else if(input$analysis_type=='Multi group'){
			actionButton(inputId='ul2ilincs', label=div("Upload signature to iLINCS",img(src="images/ilincs_logo.png", height = "20px", width="20px")), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
		} else if(input$analysis_type=='Two group with covariate'){
			if (is.null(design_rank())) {
				actionButton(inputId='ul2ilincs', label=div("Upload signature to iLINCS",img(src="images/ilincs_logo.png", height = "20px", width="20px")), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
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
		if (input$analysis_type=='Multi group') {
			shinyjs::hide("dynamiclink")
		}
	})

	####### modal Heatmap in regular signature tab #######

	vheatfit_regsig <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatfit_regsig, {
		vheatfit_regsig$doPlot <- input$draw_heatfit_regsig
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatfit_regsig$doPlot <- FALSE
		}
	}) 
	observeEvent(input$cr_reg_sig, {
		vheatfit_regsig$doPlot <- FALSE
	})
	
	vheatscr_regsig <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatscr_regsig, {
		vheatscr_regsig$doPlot <- input$draw_heatscr_regsig
	})	
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatscr_regsig$doPlot <- FALSE
		}
	}) 
	observeEvent(input$cr_reg_sig, {
		vheatscr_regsig$doPlot <- FALSE
	})
		
	vheatint_regsig <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatint_regsig, {
		vheatint_regsig$doPlot <- input$draw_heatint_regsig
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatint_regsig$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$cr_reg_sig, {
		vheatint_regsig$doPlot <- FALSE
	})

	output$heat_typeui_mod_regsig <- renderUI({	
		radioButtons("heat_type_mod_regsig", h4("Type of heatmap"), c("Heatmap across all the samples"="all","Heatmap across the comparison samples"="notall"),
			selected="notall")
	})
	output$heat_prop_mod_regsig <- renderUI({	
		col_names <- colnames(sub_meta_dat())
		checkboxGroupInput("heat_property_mod_regsig", h4("Select properties"), choices  = col_names, 
			selected=col_names[which(col_names %in% input$analysis_property)])
	})
	output$cluster_meth_mod_regsig <- renderUI({	
		selectInput("clustering_method_mod_regsig", label = h4("Grouping genes and samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected= "Pearson correlation"
		)
	})	
	output$cluster_genes_mod_regsig <- renderUI({
		numericInput("cluster_ngenes_mod_regsig", label=h4('Number of top DE genes')
			,value= 100
			,min = 1, max = 10000
		)
	})

	plot_counts_mod_regsig <- reactive({
		counts <- if(input$heat_type_mod_regsig=="all"){
					sigdata()$counts_all
				} else {
					sigdata()$counts_sig
				}
		rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
		return(counts)
	})
	data_dl_mod_regsig <- reactive({
		dge <- if(input$heat_type_mod_regsig=="all"){
					DGEList(counts=sigdata()$counts_all)
				} else {
					DGEList(counts=sigdata()$counts_sig)
				}
		dge_cpm <- cpm(dge, log=FALSE)
		expr_data<- data.matrix(dge_cpm)
		expr_data_in <- expr_data[1:input$cluster_ngenes_mod_regsig,]
		gene_sym <- counts_dat()[match(rownames(expr_data_in), counts_dat()[,1]),1,drop=F]
		dat2 <- data.frame(GeneSymbol=gene_sym[,1], expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	})

	output$fitscreen_mod_regsig = renderPlot({
		if (vheatfit_regsig$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod_regsig=="all"){
					sigdata()$counts_all
				} else {
					sigdata()$counts_sig
				}
			rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
			complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[which(rownames(sub_meta_dat()) %in% colnames(counts)),,drop=F]
			, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
			, type="Fit in Screen")
		})
	}
		,height = eventReactive (input$draw_heatfit_regsig,{
			df <- sub_meta_dat()[,which(colnames(sub_meta_dat()) %in% input$heat_property_mod_regsig),drop=FALSE]
			df[is.na(df)] <- "NA"
			df[] <- lapply( df, factor)
			colnames(df) <- input$heat_property_mod_regsig
			x=sapply(df, function(x) length(unique(x)))
			nrow = max(x)/2
			if(input$cluster_ngenes_mod_regsig<=50)
				return((length(input$heat_property_mod_regsig)*25) + (160+ (input$cluster_ngenes_mod_regsig*5.50))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
			else {
				return((length(input$heat_property_mod_regsig)*25) + 500 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
			}
		})
	)
	
	output$downloadfit_mod_regsig <- downloadHandler(
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
				counts <- if(input$heat_type_mod_regsig=="all"){
						sigdata()$counts_all
					} else {
						sigdata()$counts_sig
					}
				rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
				complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
				, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
				, type="Fit in Screen")
				dev.off()
			} else if(input$download_type_fit_mod_regsig == "csv"){
				dat2 <- data_dl_mod_regsig()
				write.csv(dat2, file)
			} else {
				return()
			}
		}
	)

	#### Scrollable heatmap

	output$scroll_mod_regsig = renderPlot({
		if (vheatscr_regsig$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod_regsig=="all"){
					sigdata()$counts_all
				} else {
					sigdata()$counts_sig
				}
			rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
			complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
			, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
			, type="Scrollable")
		})
	}
		,height = eventReactive (input$draw_heatscr_regsig,{
			if(input$cluster_ngenes_mod_regsig>=100) {
				value <- max((length(input$heat_property_mod_regsig)*25)+(10*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*24))
			} else {
				value <- (length(input$heat_property_mod_regsig)*25) + 400+(input$cluster_ngenes_mod_regsig*20)
			}
			return(value)				
		})
	)
	
	output$downloadscr_mod_regsig <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scr_mod_regsig=="png"){
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
			if(input$download_type_scr_mod_regsig == "png") {
				height = reactive({
					if(input$cluster_ngenes_mod_regsig>=1000) {
						value <- max((length(input$heat_property_mod_regsig)*25)+(4*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property_mod_regsig)*25)+(10*96),(length(input$heat_property_mod_regsig)*25)+(input$cluster_ngenes_mod_regsig*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				counts <- if(input$heat_type_mod_regsig=="all"){
					sigdata()$counts_all
				} else {
					sigdata()$counts_sig
				}
				rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
				complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
				, property=input$heat_property_mod_regsig, input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig
				, type="Scrollable")
				dev.off()
			} else{
				dat2 <- data_dl_mod_regsig()
				write.csv(dat2, file)
			}	
		}
	)

	### Interactive heatmap
	
	height_in_mod_regsig = eventReactive (input$draw_heatint_regsig,{
		df <- sub_meta_dat()[,which(colnames(sub_meta_dat()) %in% input$heat_property_mod_regsig),drop=FALSE]
		df[is.na(df)] <- "NA"
		df[] <- lapply( df, factor)
		colnames(df) <- input$heat_property_mod_regsig
		x=sapply(df, function(x) length(unique(x)))
		nrow = max(x)/2
		if(input$cluster_ngenes_mod_regsig<=50)
			return((length(input$heat_property_mod_regsig)*25) + (160+ (input$cluster_ngenes_mod_regsig*5.5))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
		else {
			return((length(input$heat_property_mod_regsig)*25) + 650 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
		}
	})
	
	output$iheatmapui_mod_regsig = renderUI({
		if (vheatint_regsig$doPlot == FALSE ) return()
		isolate({
			iheatmaprOutput("iheatmap_mod_regsig", width="auto", height = height_in_mod_regsig()) %>% withSpinner(type=5)
		})
	})
	
	output$iheatmap_mod_regsig = renderIheatmap({
		if (vheatint_regsig$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod_regsig=="all"){
				sigdata()$counts_all
			} else {
				sigdata()$counts_sig
			}
			rownames(counts) <- paste0(as.character(sigdata()$signaturesData_final[,1]), " : ", as.character(sigdata()$signaturesData_final[,2]))
			interactive_heatmap (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F], 
			property=input$heat_property_mod_regsig, n_genes=input$cluster_ngenes_mod_regsig, clustering_method=input$clustering_method_mod_regsig, signature=TRUE)
		})
	})

	
	### MA plot
	
	vma_regsig_regsig <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_ma_regsig, {
		vma_regsig_regsig$doPlot <- input$draw_ma_regsig
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vma_regsig_regsig$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$cr_reg_sig, {
		vma_regsig_regsig$doPlot <- FALSE
	}) 	

	output$fdrui_mod_regsig <- renderUI({
		numericInput("fdr_mod_regsig", label=h4('False discovery rate')
			,value= 0.05, min = 0, max = 0.5 , step= 0.01
		)
	})
	output$maplotui_mod_regsig = renderUI({
		plotOutput("maplot_mod_regsig", click = "plot_click",brush = brushOpts(id = "plot_brush")) %>% withSpinner(type=5)
	})
	
	
	output$maplot_mod_regsig <- renderPlot({	
		if (vma_regsig_regsig$doPlot == FALSE) return()
		isolate({
			res <- sigdata()$signaturesData_final
			res$Adjusted_pvalue <- res$Adjusted_pvalue #round(as.numeric(res$Adjusted_pvalue), 6)
			plot(res[, "logCPM"], res[, "Log_FoldChange"]
				,cex=0.75, pch = "*", col = ifelse(res[, "Adjusted_pvalue"] < input$fdr_mod_regsig, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
			legend("topright", legend=paste0("FDR < ",input$fdr_mod_regsig), text.col = 2, bty = "n")
			abline(h=c(-1,1), col="blue")
		})
	})	
	
	observeEvent(input$plot_click, {
		res1 <- sigdata()$signaturesData_final
		res1$Adjusted_pvalue <- res1$Adjusted_pvalue #round(as.numeric(res1$Adjusted_pvalue), 6)
		if(!is.null(input$plot_click)){
			res <- nearPoints(sigdata()$signaturesData_final, input$plot_click, "logCPM", "Log_FoldChange")
		} 
		if (nrow(res) == 0)
			return()
		output$maplot_mod_regsig <- renderPlot({
			if (vma_regsig_regsig$doPlot == FALSE) return()
			isolate({
				plot(res1[, "logCPM"], res1[, "Log_FoldChange"]
					,cex=0.75, pch = "*", col = ifelse(res1[, "Adjusted_pvalue"] < input$fdr_mod_regsig, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
				points(res1[which(res1[,1] %in% res[,1]),"logCPM"], res1[which(res1[,1] %in% res[,1]),"Log_FoldChange"]
					,cex=0.75, pch = "*",col="black")
				legend("topright", legend=paste0("FDR < ",input$fdr_mod_regsig), text.col = 2, bty = "n")
				abline(h=c(-1,1), col="blue")
			})
		}) 
		
		output$maplot_clickedpoints <- DT::renderDataTable({		
			
			datatable(
				res[,which(!colnames(res) %in% c("PValue","logCPM"))]
				,selection = 'none'
				,rownames = FALSE
				,extensions = 'Buttons'
				,options = list(dom = 'Bfrtip', pageLength = 10, scrollX = TRUE, buttons = c('csv')
					,columnDefs = list(list(className = 'dt-center', targets ="_all"))
				)			 
			)
		})
    })
	
	observeEvent(input$plot_brush, {
		res1 <- sigdata()$signaturesData_final
		res1$Adjusted_pvalue <- res1$Adjusted_pvalue #round(as.numeric(res1$Adjusted_pvalue), 6)
		if(!is.null(input$plot_brush)){
			res <- brushedPoints(sigdata()$signaturesData_final, input$plot_brush, "logCPM", "Log_FoldChange")
		} 
		if (nrow(res) == 0)
			return()
		output$maplot_mod_regsig <- renderPlot({
			if (vma_regsig_regsig$doPlot == FALSE) return()
			isolate({
				plot(res1[, "logCPM"], res1[, "Log_FoldChange"]
					,cex=0.75, pch = "*", col = ifelse(res1[, "Adjusted_pvalue"] < input$fdr_mod_regsig, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
				points(res1[which(res1[,1] %in% res[,1]),"logCPM"], res1[which(res1[,1] %in% res[,1]),"Log_FoldChange"]
					,cex=0.75, pch = "*",col="black")
				legend("topright", legend=paste0("FDR < ",input$fdr_mod_regsig), text.col = 2, bty = "n")
				abline(h=c(-1,1), col="blue")
			})
		}) 
		output$maplot_clickedpoints <- DT::renderDataTable({
			
			datatable(
				res[,which(!colnames(res) %in% c("PValue","logCPM"))]
				,selection = 'none'
				,rownames = FALSE
				,extensions = 'Buttons'
				,options = list(dom = 'Bfrtip', pageLength = 10, scrollX = TRUE, buttons = c('csv')
					,columnDefs = list(list(className = 'dt-center', targets ="_all"))
				)			 
			)
		})		
	
	})

	output$downloadMA = downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			return(paste0(datasetname,"_MAplot.png"))
		},
		content = function(file) {
			res <- sigdata()$signaturesData_final
			res$Adjusted_pvalue <- res$Adjusted_pvalue #round(as.numeric(res$Adjusted_pvalue), 6)
			res$threshold <- as.factor(res$Adjusted_pvalue < input$fdr_mod_regsig)

			png(file)
			plot(res[, "logCPM"], res[, "Log_FoldChange"]
				,cex=0.75, pch = "*", col = res$threshold, xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
			abline(h=c(-1,1), col="blue")
			dev.off()
		}	
	)

	observe({
		if (input$heatmodal_regsig=="static_mod_regsig" && input$fit_or_scr_mod_regsig=="fitinscreen_mod_regsig"){
			shinyjs::show("downloadfit_mod_regsig")
			shinyjs::show("download_type_fit_mod_regsig")
		} else {
			shinyjs::hide("downloadfit_mod_regsig")
			shinyjs::hide("download_type_fit_mod_regsig")
		}
	})

	############################ User made experimental design ############################
	
	metadata_user_exp <- reactive({
		#metadata <- meta_dat()$sub_meta
		metadata <- sub_meta_dat()
		metadata <- metadata[, input$user_analysis_prop,drop=FALSE]
		metadata[is.na(metadata)] <- "NA"
		metadata[] <- lapply( metadata, factor)
		colnames(metadata) <- input$user_analysis_prop
		rownames(metadata) <- paste0(rownames(metadata), " : ",metadata[,input$user_analysis_prop])
		return(metadata)
	})

	observeEvent(input$geo_acc2, {
		output$user_meta_prop <- renderUI ({
			metadata <- sub_meta_dat()
			col_names <- colnames(metadata)
			selectInput(inputId="user_analysis_prop",label=h4('Variable of interest'),
				choices=c('',col_names), selected=NULL
			)
		})
	})
	output$user_g1_ui <- renderUI ({
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
			#need(input$signature=="exp_desg",""),
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
			need(input$user_analysis_prop!='',"")
		)
		selectInput(inputId="user_analysis_type",label=h4('Type of comparison'),
			choices= 
				if(ncol(sub_meta_dat())>1){
					c('Two group without covariate','Two group with covariate')
				} else {
					'Two group without covariate'
				},
			selected= 'Two group without covariate'
		)
	})
	output$user_cov_ui <- renderUI({
		validate(
			need(input$user_analysis_prop!='',"")
			,need(input$user_analysis_type=="Two group with covariate","")
			,need(input$user_g2!="","")
		)
		dat <- sub_meta_dat()
		selectInput(inputId="user_cov",label=h4('Select covariates'), 
			choices= colnames(dat)[which(colnames(dat)!=input$user_analysis_prop)]
			, multiple = TRUE, selected= NULL
		)
	})
	output$gen_sig_btn_ui <- renderUI({
		validate(
			need(input$user_analysis_prop!='',"")
			,need(input$user_g1!="","")
			,need(input$user_g2!="","")
			,need(all(design_user_levels()>=2), "")
		)
		if(input$user_analysis_type=="Two group without covariate"){
			actionButton("gen_sig_btn", "Generate signature",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
		} else if(input$user_analysis_type=="Two group with covariate" & !is.null(input$user_cov)){
			if(is.null(design_user_rank())){
				actionButton("gen_sig_btn", "Generate signature",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
			} else {
				return()
			}
		} else {
			return()
		}
	})	
	
	output$sig2_warn <- renderText({
		validate(
			need(nrow(design_tab())>100,""),
			need(input$gen_sig_btn==0,""),
			need(input$user_g2!="","")
		)
		paste0("This dataset has ",nrow(design_tab())," sample runs and might take longer to generate signature.")
	})
	observeEvent(input$gen_sig_btn,{
		shinyjs::hide("sig2_warn")
	})

	metadata_user_exp_cov <- reactive({
		#metadata <- meta_dat()$sub_meta
		metadata <- sub_meta_dat()
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
		} else {
			dat1 <- metadata_user_exp_cov()[rownames(metadata_user_exp_cov()) %in% input$user_g1,,drop=FALSE]
			dat1$"Selected groups"[rownames(dat1)==input$user_g1] <- "experimental"
			dat2 <- metadata_user_exp_cov()[rownames(metadata_user_exp_cov()) %in% input$user_g2,,drop=FALSE]
			dat2$"Selected groups"[rownames(dat2)==input$user_g2] <- "control"
			dat <- rbind(dat1,dat2)
			dat <- dat[ ,c("Selected groups", input$user_analysis_prop, input$user_cov),drop=F]
			return(dat)
		} 
	})
	design_user_levels <- reactive({
		if(input$user_analysis_type=="Two group without covariate"){
			metad <- design_tab()
			design_data <- metad[,-1,drop=F]
			nlev <- as.numeric(apply(design_data, 2, function(x)length(unique(x))))
			return(nlev)
		} else if(input$user_analysis_type=="Two group with covariate"){
			metad <- design_tab()
			design_data <- metad[,-1,drop=F]
			nlev <- as.numeric(apply(design_data, 2, function(x)length(unique(x))))
			return(nlev)
		} else {
			return()
		}
	})
	design_user_rank <- reactive({
		if(all(design_user_levels()>=2)){
			metad <- design_tab()
			analysis_metadata <- metad[which(metad[,"Selected groups"]=="experimental" | metad[,"Selected groups"]=="control"),c("Selected groups",input$user_cov),drop=FALSE]
			group <- analysis_metadata[,"Selected groups", drop=F]
			covar <- analysis_metadata[,input$user_cov, drop=F]
			design_data <- cbind(covar, group)
			design_data[] <- lapply(design_data, factor)
		
			design <- model.matrix(~., data=design_data)
			ne <- nonEstimable(design)	
			return(ne)
		} else {
			return()
		}
	})
	
	output$user_rank_warn <- renderText({
		validate(
			need(input$user_analysis_type=='Two group with covariate',"")
			#,need(input$user_cov!="","")
			,need(all(design_user_levels()>=2), "")
			,need(!is.null(design_user_rank()), "")
		)	
		paste("Design matrix is not of full rank. The following coefficients are not estimable:\n", paste(design_user_rank(), collapse = " "))
	})
	
	output$user_level_warn <- renderText({
		validate(
			need(any(design_user_levels()<2), "")
			,need(nrow(design_tab())>0,"")
		)	
		paste("All the variables need to have 2 or more levels")
	})

	output$designtable = DT::renderDataTable({
		datatable( design_tab()
			,selection = 'none'
			,options = list(pageLength = if(nrow(design_tab())>=10) {10} else {nrow(design_tab())}, scrollX = TRUE)
			,rownames= TRUE
		)
	}, server = FALSE) 
	
	sigdata_user <- reactive({
		
		filepath <- paste('data/',toupper(input$geo_acc2),'/',toupper(input$geo_acc2),'_signatureData.txt', sep='')
		analysis_genes <- data.frame(Ensembl_ID=rownames(counts_dat()), Gene_symbol=counts_dat()[,1])
		metadata <- design_tab()
		rownames(metadata) <- sapply(strsplit(rownames(design_tab())," : "), `[`, 1)
		if(input$user_analysis_type=='Two group without covariate'){
			signature_data= analysis_2g_wocov(counts=counts_dat()[,-1], metadata=metadata, genes=analysis_genes, property="Selected groups",
				group1="experimental", group2="control")
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,"Log_FoldChange"], Significance_pvalue=signature_data$signaturesData_final[,"Adjusted_pvalue"], stringsAsFactors=F)
			write.table(ul_sigdata, filepath, row.names=FALSE, quote=F, append=F, sep="\t" )
		} else if(input$user_analysis_type=='Two group with covariate'){
			signature_data= analysis_2g_withcov(counts=counts_dat()[,-1], metadata=metadata, genes=analysis_genes, property="Selected groups",
				covariate=input$user_cov, group1="experimental", group2="control")
			ul_sigdata <- data.frame(Name_GeneSymbol=toupper(signature_data$signaturesData_final[,2]), Value_LogDiffExp= signature_data$signaturesData_final[,"Log_FoldChange"], Significance_pvalue=signature_data$signaturesData_final[,"Adjusted_pvalue"], stringsAsFactors=F)			
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
			sigdata_user()$signaturesData_final[,which(!colnames(sigdata_user()$signaturesData_final) %in% c("logCPM", "PValue"))]
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
				actionButton(inputId='user_ul2ilincs', label=div("Upload signature to iLINCS",img(src="images/ilincs_logo.png", height = "20px", width="20px")), style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
					onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
				)
			} else {
				return()
			}
		} else if(input$user_analysis_type=='Two group with covariate'){
			actionButton(inputId='user_ul2ilincs', label=div("Upload signature to iLINCS",img(src="images/ilincs_logo.png", height = "20px", width="20px")), style="color: #fff; background-color: #337ab7; border-color: #2e6da4", #icon = icon("paper-plane")
				onclick =paste0("window.open('http://www.ilincs.org/ilincs/signatures/main/upload?source=GRIN&fileName=",toupper(input$geo_acc2),"_signatureData.txt","', '_blank')")
			)
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

	vheatfit_mod <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatfit, {
		vheatfit_mod$doPlot <- input$draw_heatfit
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatfit_mod$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$expr_design, {
		vheatfit_mod$doPlot <- FALSE
	}) 	

	vheatscr_mod <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatscr_mod, {
		vheatscr_mod$doPlot <- input$draw_heatscr_mod

	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatscr_mod$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$expr_design, {
		vheatscr_mod$doPlot <- FALSE
	}) 	

	vheatint_mod <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_heatint_mod, {
		vheatint_mod$doPlot <- input$draw_heatint_mod
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vheatint_mod$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$expr_design, {
		vheatint_mod$doPlot <- FALSE
	}) 	

	output$heat_typeui_mod <- renderUI({	
		radioButtons("heat_type_mod", h4("Type of heatmap"), c("Heatmap across all the samples"="all", "Heatmap across the comparison samples"="notall"),
			selected="notall")
	})
	output$heat_prop_mod <- renderUI({	
		col_names <- colnames(sub_meta_dat())
		checkboxGroupInput("heat_property_mod", h4("Select properties"), choices  = col_names, 
			selected=col_names[which(col_names %in% input$user_analysis_prop)])
	})

	output$cluster_meth_mod <- renderUI({	
		selectInput("clustering_method_mod", label = h4("Grouping genes and samples"),
			choices = list("Group by properties" , "Pearson correlation", "Euclidean distance"), 
			selected="Pearson correlation"
		)
	})	
	output$cluster_genes_mod <- renderUI({
		numericInput("cluster_ngenes_mod", label=h4('Number of top DE genes')
			,value= 100
			,min = 1, max = 10000
		)
	})

	plot_counts_mod <- reactive({
		counts <- if(input$heat_type_mod=="all"){
			sigdata_user()$counts_all
		} else {
			sigdata_user()$counts_sig
		}
		rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
		return(counts)
	})

	data_dl_mod <- reactive({
		dge <- if(input$heat_type_mod=="all"){
				DGEList(counts=sigdata_user()$counts_all)
			} else {
				DGEList(counts=sigdata_user()$counts_sig)
			}
		dge_cpm <- cpm(dge, log=FALSE)
		expr_data<- data.matrix(dge_cpm)
		expr_data_in <- expr_data[1:input$cluster_ngenes_mod,]
		gene_sym <- counts_dat()[match(rownames(expr_data_in), counts_dat()[,1]),1,drop=F]
		dat2 <- data.frame(GeneSymbol=gene_sym[,1], expr_data_in, stringsAsFactors=F, check.names=F)
		return(dat2)
	})

	output$fitscreen_mod = renderPlot({
		if (vheatfit_mod$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod=="all"){
				sigdata_user()$counts_all
			} else {
				sigdata_user()$counts_sig
			}
			rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
			complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[which(rownames(sub_meta_dat()) %in% colnames(counts)),,drop=F]
			, property=input$heat_property_mod, input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod
			, type="Fit in Screen")

		})
	}
		,height = eventReactive (input$draw_heatfit,{
			if(!is.null(input$heat_property_mod)) {
				df <- sub_meta_dat()[,input$heat_property_mod,drop=FALSE]
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

	output$downloadfit_mod <- downloadHandler(
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
				counts <- if(input$heat_type_mod=="all"){
						sigdata_user()$counts_all
					} else {
						sigdata_user()$counts_sig
					}
				rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
				complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
				, property=input$heat_property_mod, input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod
				, type="Fit in Screen")
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

	output$scroll_mod = renderPlot({
		if (vheatscr_mod$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod=="all"){
					sigdata_user()$counts_all
				} else {
					sigdata_user()$counts_sig
				}
			rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
			complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
			, property=input$heat_property_mod, input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod
			, type="Scrollable")
		})
	}
		,height = eventReactive (input$draw_heatscr_mod,{
			if(input$cluster_ngenes_mod>=100) {
				value <- max((length(input$heat_property_mod)*25)+(10*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*24))
			} else {
				value <- (length(input$heat_property_mod)*25) + 400+(input$cluster_ngenes_mod*20)
			}
			return(value)				
		})
	)
	
	output$downloadscr_mod <- downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			if(input$download_type_scr_mod=="png"){
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
			if(input$download_type_scr_mod == "png") {
				height = reactive({
					if(input$cluster_ngenes_mod>=1000) {
						value <- max((length(input$heat_property_mod)*25)+(4*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*24)) #1 inch=96 pixel
					} else {
						value <- max((length(input$heat_property_mod)*25)+(10*96),(length(input$heat_property_mod)*25)+(input$cluster_ngenes_mod*30))
					}
					return(value)
				})
				png(file, width = 3000, height = height(), res=175, units = "px")
				counts <- if(input$heat_type_mod=="all"){
					sigdata_user()$counts_all
				} else {
					sigdata_user()$counts_sig
				}
				rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
				complex_heatmap_sig (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F]
				, property=input$heat_property_mod, input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod
				, type="Scrollable")
				dev.off()
			} else{
				dat2 <- data_dl_mod()
				write.csv(dat2, file)
			}	
		}
	)

	### Interactive
	
	height_in_mod = eventReactive (input$draw_heatint_mod,{
		df <- sub_meta_dat()[,which(colnames(sub_meta_dat()) %in% input$heat_property_mod),drop=FALSE]
		df[is.na(df)] <- "NA"
		df[] <- lapply( df, factor)
		colnames(df) <- input$heat_property_mod
		x=sapply(df, function(x) length(unique(x)))
		nrow = max(x)/2
		if(input$cluster_ngenes_mod<=50)
			return((length(input$heat_property_mod)*25) + (160+ (input$cluster_ngenes_mod*5.5))+(nrow*10)+ (max(nchar(rownames(df)))*6.8))
		else {
			return((length(input$heat_property_mod)*25) + 650 +(nrow*10)+ (max(nchar(rownames(df)))*6.8))
		}
	})

	output$iheatmapui_mod = renderUI({
		if (vheatint_mod$doPlot == FALSE ) return()
		isolate({
			iheatmaprOutput("iheatmap_mod", width="auto", height = height_in_mod()) %>% withSpinner(type=5)
		})
	})
	output$iheatmap_mod = renderIheatmap({
		if (vheatint_mod$doPlot == FALSE ) return()
		isolate({
			counts <- if(input$heat_type_mod=="all"){
				sigdata_user()$counts_all
			} else {
				sigdata_user()$counts_sig
			}
			rownames(counts) <- paste0(as.character(sigdata_user()$signaturesData_final[,1]), " : ", as.character(sigdata_user()$signaturesData_final[,2]))
			interactive_heatmap (data_mat=counts, metadata=sub_meta_dat()[colnames(counts),,drop=F], 
			property=input$heat_property_mod, n_genes=input$cluster_ngenes_mod, clustering_method=input$clustering_method_mod, signature=TRUE)
		})
	})
	
	### MA plot
	
	vma <- reactiveValues(doPlot = FALSE)
	observeEvent(input$draw_ma, {
		vma$doPlot <- input$draw_ma
	})
	observeEvent(input$grin, {
		if(input$grin=="datasets") {
			vma$doPlot <- FALSE
		}
	}) 	
	observeEvent(input$expr_design, {
		vma$doPlot <- FALSE
	}) 	

	output$fdrui_mod <- renderUI({
		numericInput("fdr_mod", label=h4('False discovery rate')
			,value= 0.05, min = 0, max = 0.5 , step= 0.01
		)
	})
	output$maplotui_mod = renderUI({
		plotOutput("maplot_mod", click = "plot_click_mod",brush = brushOpts(id = "plot_brush_mod")) %>% withSpinner(type=5)
	})
	
	output$maplot_mod <- renderPlot({	
		if (vma$doPlot == FALSE) return()
		isolate({
			res <- sigdata_user()$signaturesData_final
			#res$Adjusted_pvalue <- round(as.numeric(res$Adjusted_pvalue), 6)
			plot(res[, "logCPM"], res[, "Log_FoldChange"]
				,cex=0.75, pch = "*", col = ifelse(res[, "Adjusted_pvalue"] < input$fdr_mod, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
			legend("topright", legend=paste0("FDR < ",input$fdr_mod), text.col = 2, bty = "n")
			abline(h=c(-1,1), col="blue")
		})
	})	

	observeEvent(input$plot_click_mod, {
		res1 <- sigdata_user()$signaturesData_final
		#res1$Adjusted_pvalue <- round(as.numeric(res1$Adjusted_pvalue), 6)
		if(!is.null(input$plot_click_mod)){
			res <- nearPoints(sigdata_user()$signaturesData_final, input$plot_click_mod, "logCPM", "Log_FoldChange")
		} 
		if (nrow(res) == 0)
			return()
		output$maplot_mod <- renderPlot({
			if (vma$doPlot == FALSE) return()
			isolate({
				plot(res1[, "logCPM"], res1[, "Log_FoldChange"]
					,cex=0.75, pch = "*", col = ifelse(res1[, "Adjusted_pvalue"] < input$fdr_mod, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
				points(res1[which(res1[,1] %in% res[,1]),"logCPM"], res1[which(res1[,1] %in% res[,1]),"Log_FoldChange"]
					,cex=0.75, pch = "*",col="black")
				legend("topright", legend=paste0("FDR < ",input$fdr_mod), text.col = 2, bty = "n")
				abline(h=c(-1,1), col="blue")
			})
		}) 
		
		output$maplot_clickedpoints_mod <- DT::renderDataTable({		
			datatable(
				res[,which(!colnames(res) %in% c("PValue","logCPM"))]
				,selection = 'none'
				,rownames = FALSE
				,extensions = 'Buttons'
				,options = list(dom = 'Bfrtip', pageLength = 10, scrollX = TRUE, buttons = c('csv')
					,columnDefs = list(list(className = 'dt-center', targets ="_all"))
				)			 
			)
		})
    })

	observeEvent(input$plot_brush_mod, {
		res1 <- sigdata_user()$signaturesData_final
		#res1$Adjusted_pvalue <- round(as.numeric(res1$Adjusted_pvalue), 6)
		if(!is.null(input$plot_brush_mod)){
			res <- brushedPoints(sigdata_user()$signaturesData_final, input$plot_brush_mod, "logCPM", "Log_FoldChange")
		} 
		if (nrow(res) == 0)
			return()
		output$maplot_mod <- renderPlot({
			if (vma$doPlot == FALSE) return()
			isolate({
				plot(res1[, "logCPM"], res1[, "Log_FoldChange"]
					,cex=0.75, pch = "*", col = ifelse(res1[, "Adjusted_pvalue"] < input$fdr_mod, 2, 8), xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
				points(res1[which(res1[,1] %in% res[,1]),"logCPM"], res1[which(res1[,1] %in% res[,1]),"Log_FoldChange"]
					,cex=0.75, pch = "*",col="black")
				legend("topright", legend=paste0("FDR < ",input$fdr_mod), text.col = 2, bty = "n")
				abline(h=c(-1,1), col="blue")
			})
		}) 
		output$maplot_clickedpoints_mod <- DT::renderDataTable({
		
			datatable(
				res[,which(!colnames(res) %in% c("PValue","logCPM"))]
				,selection = 'none'
				,rownames = FALSE
				,extensions = 'Buttons'
				,options = list(dom = 'Bfrtip', pageLength = 10, scrollX = TRUE, buttons = c('csv')
					,columnDefs = list(list(className = 'dt-center', targets ="_all"))
				)			 
			)
		})		
	
	})
	
	output$downloadMA_mod = downloadHandler(
		filename = function() {
			datasetname <- toupper(input$geo_acc)
			return(paste0(datasetname,"_MAplot.png"))
		},
		content = function(file) {
			res <- sigdata_user()$signaturesData_final
			#res$Adjusted_pvalue <- round(as.numeric(res$Adjusted_pvalue), 6)
			res$threshold <- as.factor(res$Adjusted_pvalue < input$fdr_mod)

			png(file)
			plot(res[, "logCPM"], res[, "Log_FoldChange"]
				,cex=0.75, pch = "*", col = res$threshold, xlab = "Average log(Counts per million)", ylab = "Log2(Fold change)")
			abline(h=c(-1,1), col="blue")
			dev.off()
		}	
	)

	observe({
		if (input$heatmodal=="static_mod" && input$fit_or_scr_mod=="fitinscreen_mod"){
			shinyjs::show("downloadfit_mod")
			shinyjs::show("download_type_fit_mod")
		} else {
			shinyjs::hide("downloadfit_mod")
			shinyjs::hide("download_type_fit_mod")
		}
	})


	################################################ Processing users data nav page ################################################
	
	observe({
		shinyjs::hide(selector = "#grin li a[data-value=out_cons]")
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

	vstart_process <- reactiveValues(doPlot = FALSE)
	observeEvent(input$start_process, {
		vstart_process$doPlot <- input$start_process
	})

	observeEvent(input$start_process,{
		#if (input$start_process) {
		if(vstart_process$doPlot==FALSE) return()
		isolate ({
			if(input$user_geo!=''){
				updateTabsetPanel(session, "grin", selected = "out_cons")
				user_geo = trimws(input$user_geo, which ="both")
				if(!file.exists(paste0("data/user_geo_request/", toupper(user_geo), ".txt")) & 
					input$user_geo!='' & !dir.exists(paste0("data/user_geo_request/", toupper(user_geo))))
				{
					system(paste0("touch data/user_geo_request/", toupper(user_geo), ".txt"))
				}
				
				filenames <- reactive({
					details <- file.info(list.files(path="data/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
					files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
					x=gsub("^.*/", "", sub(".txt","",files))
					return(x)
				})
				create_dir <- reactive({
					system(paste0("mkdir data/user_geo_request/", filenames()[1]))
					system(paste0("chmod -R 777 data/user_geo_request/", filenames()[1]))
					dir_name <- list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE)
					return(dir_name)
				})

				updateTextInput(session, 'geo_current', 
					value = ifelse(length(list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE))==0, create_dir(), 
						list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE))
				)
				updateSelectInput(session, 'geo_next', label="Waiting to process", choices=filenames(), selected=filenames()[1])
			}
		})
    })	
	
	output$steps <- renderText({
		paste0("Data is being processed by our back-end pipeline.")
	})

	vprocess_log <- reactiveValues(doPlot = FALSE)
	observeEvent(input$process_log, {
		vprocess_log$doPlot <- input$process_log
	})

	observeEvent(input$process_log,{
		if(vprocess_log$doPlot==FALSE) return()
		isolate ({
			updateTabsetPanel(session, "grin", selected = "out_cons")

			filenames <- reactive({
				details <- file.info(list.files(path="data/user_geo_request",pattern = ".txt$", recursive = F,full.names=T))
				files <- rownames(details[with(details, order(as.POSIXct(mtime))), ])
				x=gsub("^.*/", "", sub(".txt","",files))
				return(x)
			})
			dir_created <- reactive({
				dir_name <- list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE)
				return(dir_name)
			})

			updateTextInput(session, 'geo_current', 
				value = ifelse(length(list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE))==0, dir_created(), 
					list.dirs(path = "data/user_geo_request/", full.names = F, recursive = FALSE))
			)
			updateSelectInput(session, 'geo_next', label="Waiting to process", choices=filenames(), selected=filenames()[1])
		})
	})
	
	observeEvent(input$show_log,{
		if(input$show_log){
			fileReaderData <- reactiveFileReader(500, session, paste0("data/user_geo_request/",input$geo_current,"/log.txt"), readLines)			
			output$process_log_text <- renderText({
				validate(
					need(input$geo_current!="","")
					,need(input$show_log!=0,"")
				)
				text <- fileReaderData()
				text[is.na(text)] <- ""
				paste(text, collapse = '\n')
			})
		}
    })	
	
})