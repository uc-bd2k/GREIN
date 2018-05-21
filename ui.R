options(rgl.useNULL = TRUE)
library(shiny)	
library(Biobase)
library("RNASeqPower")
library("ComplexHeatmap")
library("circlize")
library("edgeR")
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library("rgl")
library("rsconnect")
library("shinyRGL")
library(Rtsne)
library("rglwidget")
library("RColorBrewer")
library("pairsD3")
library("data.table")
library("feather")
library("plotly")
library("ggplot2")
library("reshape2")
library(DT)
library("iheatmapr")


source("R/dropdownButton.R", local = TRUE)

shinyUI(
	navbarPage(
		inverse=TRUE, windowTitle = 'GREIN' 
		,title= tags$a(href='https://shiny.ilincs.org/grein',target="_blank", tags$img(src='images/logo1c.png', height="40px",width="100px", style='margin: -20px 0px;'))
		,id='grin', selected="datasets"

		,tabPanel(tags$b('Help'), value='help'
			,tags$style(
				'.navbar-nav {width: 90%;}
				.navbar-nav :first-child{float:right}'
			)
			,column(1				
			)
			,column(10
				,div(style="color:#337ab7; text-align:center;", titlePanel(h2("Step-by-step guide of GREIN"))), br()
				,htmlOutput('help_pdf')
			)
			,column(1
				
			)			
		)
		,tabPanel(tags$b('About'), value='about',			
			tags$head(
				tags$style(type="text/css", "
				   #loadmessage {
					 position: fixed;
					 top: 80px;
					 left: 0px;
					 width: 100%;
					 padding: 5px 0px 5px 0px;
					 text-align: center;
					 font-weight: bold;
					 font-size: 100%;
					 color: #000000;
					 background-color: #DDDDDD;
					 z-index: 105;
				   }
					  .js-irs-3 .irs-grid-text{
					  display:none;
					  }
					/* basic positioning */
					.legend { list-style: none; }
					.legend li { float: left; margin-right: 10px; }
					.legend span { border: 1px solid #ccc; float: left; width: 12px; height: 12px; margin: 2px; }
				")
				,tags$script('window.onbeforeunload = function() { return "Please use the button on the webpage"; };')
				,tags$link(rel = "icon", type = "image/ico", href = "images/grein_favicon.ico")
				,tags$style(type="text/css",
					".shiny-output-error { visibility: hidden; }",
					".shiny-output-error:before { visibility: hidden; }")
				,tags$link(href="css/chardinjs.css",rel="stylesheet")
				,tags$link(href="css/ilincs.css",rel="stylesheet")
			)
			,useShinyjs()
			,column(1
				
			)
			,column(10,
				wellPanel(
					includeMarkdown("www/description1.Rmd")
				)
			)
			,column(1
				
			)
	 		,includeHTML("www/html/footer.html")
		),
		tabPanel(tags$b('GEO datasets'), value='datasets'
			,includeMarkdown("www/title.Rmd")
			,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
			,column(4, offset=0, style='padding:0px 20px 0px 0px;'
				,tags$br()
				,plotlyOutput("samplesize_den", width = "auto", height = "300px") %>% withSpinner(type=5), br()
			)
			,column(4, offset=0, style='padding:10px 0px 0px 10px;'
				,br()
				, plotlyOutput("study_stat", width = "auto", height = "300px") %>% withSpinner(type=5), br()
			)
			,column(4, offset=0, style='padding:10px 20px 0px 0px;'
				,tags$br()
				,plotlyOutput("sample_stat", width = "auto", height = "300px") %>% withSpinner(type=5), br()
			)
			,column(12, offset=0, style='padding:0px;'
				,div(style="height:2px;border-width:0;color:gray;background-color:gray", hr())
				,div(style="color:#337ab7;", titlePanel(""))
			)
			,column(4, offset=0, style='padding:0px;'
				,tags$br()
				,div(style="display: inline-block;vertical-align:right; ",h4('Search for GEO series (GSE) accession'))
				,div(style="display: inline-block;vertical-align:right;", bsButton("user_geo_bs", label = "", icon = icon("info-circle"),style="", size = "medium"))
				,bsPopover(id = "user_geo_bs", title = "",
					content = paste0("Please type your GEO accession id (e.g., GSE46597) in this box to see if the dataset is already processed. If not, you can initialize the processing of the data. The dataset will be processed and uploaded soon."),
					placement = "right", trigger = "hover", options = list(container = "body")
				)
				,tags$style(type='text/css', '#user_geo_bs {background-color: white; border:none;}')
				,div(style="display: inline-block;vertical-align:right;line-height:0px;",textInput('user_geo', '', width='340px'))
				
				,div(style="width:400px;", verbatimTextOutput("warn"))
				,div(style="width:400px;",verbatimTextOutput("warn2"))
				,div(style="width:400px;",verbatimTextOutput("warn3"))
				,div(style="width:400px;",verbatimTextOutput("warn4"))
                ,div(style="width:400px;",verbatimTextOutput("warn5"))
				,tags$style(type='text/css', '#warn {color: red; font-size:15px;}')
				,tags$style(type='text/css', '#warn2 {color: red; font-size:15px;}')
				,tags$style(type='text/css', '#warn3 {color: red; font-size:15px;}')
				,tags$style(type='text/css', '#warn4 {color: red; font-size:15px;}')
                ,tags$style(type='text/css', '#warn5 {color: red; font-size:15px;}')
				
				,div(style="display: ;vertical-align:bottom;", shinyjs::hidden(actionButton("start_process","Start processing", icon("refresh"),
						style="color: #fff; background-color: #337ab7; border-color: #2e6da4")))
				,tags$style(type='text/css', "label {font-size: 18px;}")
			)
			,column(4, offset=0, style='padding:0px;'
				,tags$br()
				,div(style="display: inline-block;vertical-align:right; ",h4('Search by ontologies', tags$a(href='http://metasra.biostat.wisc.edu/',target="_blank","(MetaSRA)", style="color: #337ab7")))
				,bsPopover(id = "sample_ont_ser_bs", title = "",
					content = paste0("Search by ontological annotations of samples and datasets provided  by MetaSRA project."),
					placement = "right", trigger = "hover", options = list(container = "body")
				)
				,tags$style(type='text/css', '#sample_ont_ser_bs {background-color: white; border:none;}')
				,div(style="display: inline-block;vertical-align:right;", bsButton("sample_ont_ser_bs", label = "", icon = icon("info-circle"),style="", size = "medium"))
				,div(style="display: inline-block;vertical-align:right;line-height:0px;",textInput('sample_ont_ser', '', width="325px"))

				,div(style="width:320px;",verbatimTextOutput("warn_ont"))
				,div(style="width:320px;",verbatimTextOutput("warn_ont2"))
				,tags$style(type='text/css', '#warn_ont {color: red; font-size:15px;}') 
				,tags$style(type='text/css', '#warn_ont2 {color: red; font-size:15px;}') 

			)	
			,column(4, offset=0, style='padding:0px;'
				,div(style="line-height:20px", br())
				#User requested datasets in the processing queue
				,div(style="display: inline-block;vertical-align:right; ",h4('User requested datasets'))
				,div(style="display: inline-block;vertical-align:right;", bsButton("process_log_bs", label = "", icon = icon("info-circle"),style="", size = "medium"))
				,bsPopover(id = "process_log_bs", title = "",
					content = paste0("Click `Processing console` button to see status of the user requested datasets in the processing queue."),
					placement = "right", trigger = "hover", options = list(container = "body")
				)
				,tags$style(type='text/css', '#process_log_bs {background-color: white; border:none;}')			
				
				,div(style="line-height:0px;", selectInput("geo_next_lp", " ", choices="", selected=""))
				,div(style="margin:-10px;", br()), actionButton("process_log","Processing console", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
			)
			,column(12, offset=0, style='padding:0px;'
				,tags$hr()
				,div(style="color:#337ab7; text-align:center;", titlePanel(h1("Select a dataset")))
			)
			,DT::dataTableOutput('datatable'), br()
			,includeHTML("www/html/footer.html")
		)
		,tabPanel(tags$b('Explore dataset'), value='explore'
			,column(3, offset=0
				,wellPanel(
					textInput('geo_acc', h3('Selected study'))
					,conditionalPanel("input.geo_acc != ''"
						,actionButton("analysisbtn", "Analyze", icon("th"), style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #337ab7; border-color: #2e6da4")
						,conditionalPanel("input.tab2 == 'metadata'"
							,br(),tags$div(class="metadata_column_class", uiOutput("metadata_column"))
							,tags$head(tags$style(HTML(".metadata_column_class label{font-size: 16px;}")))
							,tags$div(class="full_metacol_class", uiOutput("full_metacol"))
							,tags$head(tags$style(HTML(".full_metacol_class label{font-size: 16px;}")))
							,br(), shinyjs::hidden(downloadButton('downloadmeta', label = "Download metadata", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
						)
						,conditionalPanel("input.tab2 == 'countsTable'"
							
							,br(),uiOutput("counts_sample_ui")
							,actionButton(inputId = "load_counts", label = "Show counts table",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
							,conditionalPanel("input.load_counts != ''"
								,br(),br(),downloadButton('downloadcounts', label = "Download data", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
								,tags$div(class="dl_counts_type_class", radioButtons('download_type_counts', label = "", choices = c("Gene level", "Transcript level"), inline=T))
								,tags$head(tags$style(HTML(".dl_counts_type_class label{font-size: 16px;}")))
							)
						)
						,conditionalPanel("input.preplots == 'corr' | input.preplots == 'dens'"
							,br(),br(),uiOutput("cor_sample_ui")
							,uiOutput("draw_cor_ui")
							,uiOutput("draw_dens_ui")
						)
						,conditionalPanel("input.preplots == 'heatmap'"
							,br(),tags$div(class="heat_prop_class", uiOutput("heat_prop"))
							,tags$head(tags$style(HTML(".heat_prop_class label{font-size: 16px;}")))
							,uiOutput("cluster_meth")
							,uiOutput("cluster_genes")
							,conditionalPanel("input.heatvis=='static' && input.fir_or_scr=='fitinscreen'"	
								,actionButton("draw_heat", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
								,conditionalPanel("input.draw_heat != ''"
									,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('download', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
									,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit', label = "", choices = c("png", "csv"), inline = T))								
								)
							)
							,conditionalPanel("input.heatvis=='static' && input.fir_or_scr=='scrollable'"
								,actionButton("draw_heatscr", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
								,conditionalPanel("input.draw_heatscr != ''"
									,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('downloadscr', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
									,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scr', label = "", choices = c("png", "csv"), inline = T))							
								)
							)
							,conditionalPanel("input.heatvis=='inter'"
								,actionButton("draw_heatint", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
							)
						)
						,conditionalPanel("input.tab2 == 'plots' && input.preplots == 'pca'"						
							,br(),uiOutput("pca_prop")							
							,conditionalPanel("input.pcap == 'pca2'"						
								,actionButton("draw_pca", "Plot PCA",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
								,br(), plotOutput("legendplot", height='auto')
							)
							,conditionalPanel("input.pcap == 'pca3'"													
								,actionButton("draw_pca3", "Plot PCA",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
								,br(), plotOutput("legendplot3", height='auto')
							)
							,conditionalPanel("input.pcap == 'tsne2'"													
								,actionButton("draw_tsne2", "Plot t-SNE",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
							)							
							,conditionalPanel("input.pcap == 'tsne3'"													
								,actionButton("draw_tsne3", "Plot t-SNE",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
							)							
						)
					)
				)
			)
			,column(9, offset=0
				,wellPanel(
					conditionalPanel("input.geo_acc != ''"
						,tabsetPanel( id="tab2"
							,tabPanel("Description", value="summary", br()
								,DT::dataTableOutput('geo_summary') %>% withSpinner(type=5)
							)
							,tabPanel("Metadata",value="metadata", br()
								,DT::dataTableOutput('metadata') %>% withSpinner(type=5)
								,DT::dataTableOutput('metadata_full', width = 900)
								,DT::dataTableOutput('ontology', width = 900)
								)
							,tabPanel("Counts table",value="countsTable", br()
								,verbatimTextOutput("counts_warn")
								,tags$style(type='text/css', '#counts_warn {color: red; font-weight: bold; font-size:18px;}')
								,DT::dataTableOutput('counts', width=900) %>% withSpinner(type=5)
							)
							,tabPanel("QC report",value="qc", br()
								,downloadButton('downloadqc', label = "Download QC report", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
								,htmlOutput('multiqc') %>% withSpinner(type=5)
							)
							,tabPanel("Visualization",value="plots", br()
								,tabsetPanel(id="preplots"
									,tabPanel("Correlation plot",value="corr", br()
										,verbatimTextOutput("corr_warn")
										,tags$style(type='text/css', '#corr_warn {color: red; font-weight: bold; font-size:18px;}')
										,plotlyOutput("corplot") %>% withSpinner(type=5)
									)
									,tabPanel("Density plot",value="dens", br()
										,verbatimTextOutput("dens_warn")
										,tags$style(type='text/css', '#dens_warn {color: red; font-weight: bold; font-size:18px;}')
										,plotlyOutput("densityp", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
									,tabPanel("Heatmap",value="heatmap"
										,br() 
										,verbatimTextOutput("heat_warn")
										,tags$style(type='text/css', '#heat_warn {color: red; font-weight: bold; font-size:18px;}')
										,tabsetPanel(id="heatvis"
											,tabPanel(div(style="display: inline-block;vertical-align:right;","Static heatmap", bsButton("heat_stat_bs", label = "", icon = icon("info-circle"),style="", size = "small")), value="static"
												,bsPopover(id = "heat_stat_bs", title = "",
													content = paste0("Highly variable genes are selected based on median absolute deviation (MAD) values of the genes."),
													placement = "bottom", trigger = "hover", options = list(container = "body")
												)	
												,tags$style(type='text/css', '#heat_stat_bs {background-color: transparent; border:none;}')
												,br()
												,tabsetPanel(id="fir_or_scr"
													,tabPanel("Fit in screen",value="fitinscreen", br()
														,plotOutput("HeatStatic", height= "auto") %>% withSpinner(type=5)
													)
													,tabPanel("Scrollable",value="scrollable", br()
														,plotOutput("HeatStaticscr", height= "auto") %>% withSpinner(type=5)
													)
												)	
											)
											,tabPanel(div(style="display: inline-block;vertical-align:right;","Interactive heatmap", bsButton("heat_int_bs", label = "", icon = icon("info-circle"),style="", size = "small")), value="inter"
												,bsPopover(id = "heat_int_bs", title = "",
													content = paste0("You can select any area in the following heatmap to zoom in and double-click to zoom out."),
													placement = "bottom", trigger = "hover", options = list(container = "body")
												)	
												,tags$style(type='text/css', '#heat_int_bs {background-color: transparent; border:none;}')												
												,br()
												,uiOutput("iheatmapui", width = "auto", height = "auto"), br(), br()
											)
										)
									),
									tabPanel("PCA & t-SNE",value="pca", br()
										,tabsetPanel(id="pcap"
											,tabPanel("2D PCA plot",value="pca2"	
												,br()
												,conditionalPanel("input.draw_pca != ''"
													,downloadButton('downloadPCA2D', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
												)
												,br(),uiOutput("p2D_ui")
												,br()
											)
											,tabPanel("3D PCA plot",value="pca3"
												,br()
												,conditionalPanel("input.draw_pca3 != ''"
													,downloadButton('downloadPCA3D', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
												)
												,br(), uiOutput("pca3Dplot_ui")
												,br()
											)
											,tabPanel("2D t-SNE plot",value="tsne2"
												,br()
												,uiOutput("t2Dplot_ui")
											)											
											,tabPanel("3D t-SNE plot",value="tsne3"
												,br()
												,uiOutput("tsne3Dplot_ui")
											)
										)
									)
								)
							)
						)
					)
				)
			)
			,includeHTML("www/html/footer.html")
		)
		,tabPanel(tags$b('Analyze dataset'), value='analyze'
			,column(3, offset=0.5, 
				wellPanel(
					textInput('geo_acc2', h3('Selected study'))
					,conditionalPanel("input.geo_acc != ''"
						,conditionalPanel("input.signature == 'power'"
							,uiOutput("power_ui")
							,uiOutput("samples")
							,uiOutput("alpha")
							,uiOutput("effect")	
							,conditionalPanel("input.sigtab== 'power_curve'"
								,uiOutput("depth_ui")
							)
							,br(), uiOutput("gen_pow_btn_ui")
							,uiOutput("gen_detect_btn_ui")
						)
						,conditionalPanel("input.signature == 'sig_data'"
							,conditionalPanel("input.cr_reg_sig == 'meta_reg'"
								,uiOutput("metadata_analysis")
								,uiOutput("ana_type")
								,uiOutput("group1")
								,uiOutput("group2")
								,uiOutput("cov")
								,br(), uiOutput("gen_reg_sig_ui")
							)
							,conditionalPanel("input.signature == 'sig_data' && input.cr_reg_sig == 'sig_reg'"
								,br(),actionButton("heatmap_mod_regsig", "Show visualization", style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
								,bsModal("modalheatmap_regsig", "", "heatmap_mod_regsig", size = "large"
									,column(4, offset=-0.5
										,wellPanel(
											conditionalPanel("input.heatmodal_regsig=='static_mod_regsig' || input.heatmodal_regsig=='inter_mod_regsig'"
												,tags$div(class="heat_typeui_mod_regsig_class",uiOutput("heat_typeui_mod_regsig"))
												,tags$head(tags$style(HTML(".heat_typeui_mod_regsig_class label{font-size: 15px;}")))
												,tags$div(class="heat_prop_mod_regsig_class",uiOutput("heat_prop_mod_regsig"))
												,tags$head(tags$style(HTML(".heat_prop_mod_regsig_class label{font-size: 15px;}")))
												,uiOutput("cluster_meth_mod_regsig")
												,uiOutput("cluster_genes_mod_regsig"), br()
												
												,conditionalPanel("input.heatmodal_regsig=='static_mod_regsig' && input.fit_or_scr_mod_regsig=='fitinscreen_mod_regsig'"	
													,actionButton("draw_heatfit_regsig", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
													,conditionalPanel("input.draw_heatfit_regsig != ''"
														,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('downloadfit_mod_regsig', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
														,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit_mod_regsig', label = "", choices = c("png", "csv"), inline = T))
													)
												)
												,conditionalPanel("input.heatmodal_regsig=='static_mod_regsig' && input.fit_or_scr_mod_regsig=='scrollable_mod_regsig'"
													,actionButton("draw_heatscr_regsig", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
													,conditionalPanel("input.draw_heatscr_regsig != ''"
														,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('downloadscr_mod_regsig', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
														,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scr_mod_regsig', label = "", choices = c("png", "csv"), inline = T))							
													)
												)
												,conditionalPanel("input.heatmodal_regsig=='inter_mod_regsig'"
													,actionButton("draw_heatint_regsig", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
												)
											)
											,conditionalPanel("input.heatmodal_regsig=='ma_mod_regsig'"
												,uiOutput("fdrui_mod_regsig"), br()
												,actionButton("draw_ma_regsig", "Draw MA plot",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
											)
										)
									)
									,column(8, offset=-0.5
										,wellPanel(
											tabsetPanel(id="heatmodal_regsig",
												tabPanel("Static heatmap",value="static_mod_regsig", br()
													,conditionalPanel("input.heatmodal_regsig == 'static_mod_regsig'"
														,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("Heatmap of top DE genes (ranked by adjusted p-values)")))
													)
													,tabsetPanel(id="fit_or_scr_mod_regsig"
														,tabPanel("Fit in screen",value="fitinscreen_mod_regsig", br()
															,plotOutput("fitscreen_mod_regsig", height="auto") %>% withSpinner(type=5)
														)
														,tabPanel("Scrollable",value="scrollable_mod_regsig", br()
															,plotOutput("scroll_mod_regsig", height="auto") %>% withSpinner(type=5)
														)
													)
												),
												tabPanel("Interactive heatmap",value="inter_mod_regsig", br() 
													,conditionalPanel("input.heatmodal_regsig == 'inter_mod_regsig'"
														,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("Heatmap of top DE genes (ranked by adjusted p-values)")))
													)	
													,conditionalPanel("input.draw_heatint_regsig != ''"
														,uiOutput("iheatmapui_mod_regsig", width = "auto", height = "auto")
													)
												),
												tabPanel("MA plot",value="ma_mod_regsig", br() 
													,conditionalPanel("input.heatmodal_regsig == 'ma_mod_regsig'"
														,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("log-fold change (M) vs. log-average expression (A) plot")))
													)																								
													,downloadButton('downloadMA', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), br()
													,br(), uiOutput("maplotui_mod_regsig", width = "auto", height = "auto"), br()
													,br(), DT::dataTableOutput("maplot_clickedpoints", width=500)
												)
											)
										)
									)
								)				
								,br(), br(), downloadButton('downloadsigdata', label = "Download signature", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;")
								,br(), br(),uiOutput("dynamiclink")
							)
						)
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='sig_meta'"	
							,uiOutput("user_meta_prop")
							,uiOutput("user_ana_type")
							,uiOutput("user_cov_ui")
							,br(), uiOutput("gen_sig_btn_ui"), br()
						)
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='user_exp_sig'"
							,br(), actionButton("heatmap_mod", "Show visualization", style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
							,bsModal("modalheatmap", "", "heatmap_mod", size = "large"
								,column(4, offset=-0.5
									,wellPanel(
										conditionalPanel("input.heatmodal=='static_mod' || input.heatmodal=='inter_mod'"
											,tags$div(class="heat_type_mod_class", uiOutput("heat_typeui_mod"))
											,tags$head(tags$style(HTML(".heat_type_mod_class label{font-size: 15px;}")))
											,tags$div(class="heat_prop_mod_class",uiOutput("heat_prop_mod"))
											,tags$head(tags$style(HTML(".heat_prop_mod_class label{font-size: 15px;}")))
											,uiOutput("cluster_meth_mod")
											,uiOutput("cluster_genes_mod"), br()
											
											,conditionalPanel("input.heatmodal=='static_mod' && input.fit_or_scr_mod=='fitinscreen_mod'"
												,actionButton("draw_heatfit", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
												,conditionalPanel("input.draw_heatfit != ''"
													,br(),div(style="display: inline-block;vertical-align:right;", downloadButton('downloadfit_mod', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
													,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit_mod', label = "", choices = c("png", "csv"), inline = T))
												)
											)
											,conditionalPanel("input.heatmodal=='static_mod' && input.fit_or_scr_mod=='scrollable_mod'"
												,actionButton("draw_heatscr_mod", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
												,conditionalPanel("input.draw_heatscr_mod != ''"
													,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('downloadscr_mod', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
													,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scr_mod', label = "", choices = c("png", "csv"), inline = T))							
												)
											)
											,conditionalPanel("input.heatmodal=='inter_mod'"
												,actionButton("draw_heatint_mod", "Draw heatmap",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000")
											)
										)
										,conditionalPanel("input.heatmodal=='ma_mod'"
											,uiOutput("fdrui_mod"), br()
											,actionButton("draw_ma", "Draw MA plot",style="color: #fff; box-shadow: 1px 1px 1px 1px #888; background-color: #800000; border-color: #800000"), br()
										)
									)
								)
								,column(8, offset=-0.5
									,wellPanel(
										tabsetPanel(id="heatmodal",
											tabPanel("Static heatmap",value="static_mod", br()
												,conditionalPanel("input.heatmodal == 'static_mod'"
													,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("Heatmap of top DE genes (ranked by adjusted p-values)")))
												)
												,tabsetPanel(id="fit_or_scr_mod"
													,tabPanel("Fit in screen",value="fitinscreen_mod", br()
														,plotOutput("fitscreen_mod", height="auto") %>% withSpinner(type=5)
													)
													,tabPanel("Scrollable",value="scrollable_mod", br()
														,plotOutput("scroll_mod", height="auto") %>% withSpinner(type=5)
													)
												)								
											),
											tabPanel("Interactive heatmap",value="inter_mod", br() 
												,conditionalPanel("input.heatmodal == 'inter_mod'"
													,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("Heatmap of top DE genes (ranked by adjusted p-values)")))
												)	
												,conditionalPanel("input.draw_heatint_mod != ''"
													,uiOutput("iheatmapui_mod", width = "auto", height = "auto")
												)
											),
											tabPanel("MA plot",value="ma_mod", br() 
												,conditionalPanel("input.heatmodal == 'ma_mod'"
													,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h4("log-fold change (M) vs. log-average expression (A) plot")))
												)																								
												,downloadButton('downloadMA_mod', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), br()
												,br(), uiOutput("maplotui_mod", width = "auto", height = "auto"), br()
												,br(), DT::dataTableOutput("maplot_clickedpoints_mod", width=500)
											)									
										)
									)
								)
							)
						)
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='user_exp_sig'"
							,br(), shinyjs::hidden(downloadButton('download_user_sigdata', label = "Download signature", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;")), br()
							,br(), uiOutput("user_dynamiclink")
						)						
					)
				)
			)
			,column(9, offset=0.5,
				wellPanel(
					conditionalPanel("input.geo_acc !='' && input.geo_acc2 !=''"
						,tabsetPanel( id="signature"
							,tabPanel("Power analysis",value="power", br()
								,tabsetPanel( id="sigtab"
									,tabPanel(div(style="display: inline-block;vertical-align:right;","Power curve", bsButton("power_bs", label = "", icon = icon("info-circle"),style="", size = "small")),value="power_curve", br()
										,bsPopover(id = "power_bs", title = "",
											content = paste0("We use Bioconductor package <strong>RNASeqPower</strong> to estimate power."),
											placement = "bottom", trigger = "hover", options = list(container = "body")
										)	
										,tags$style(type='text/css', '#power_bs {background-color: transparent; border:none;}')										
										,plotlyOutput("powerplot", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
									,tabPanel(div(style="display: inline-block;vertical-align:right;","Detectability of genes", bsButton("detect_bs", label = "", icon = icon("info-circle"),style="", size = "small")),value="detec_gene", br()
										,bsPopover(id = "detect_bs", title = "",
											content = paste0("The line of detectability is the estimated biological coefficient of variation values for different sequencing depths. Genes below this line have higher power to be detected as differentially expressed."),
											placement = "bottom", trigger = "hover", options = list(container = "body")
										)	
										,tags$style(type='text/css', '#detect_bs {background-color: transparent; border:none;}')										
										
										,div(style="display: inline-block;vertical-align:right; ",h4('Search for genes of interest'))
										,div(style="display: inline-block;vertical-align:right;", bsButton("search_gene_bs", label = "", icon = icon("info-circle"),style="", size = "small"))
										,bsPopover(id = "search_gene_bs", title = "",
											content = paste0("You can type one or more gene symbols seperated by commas (e.g., TNMD,FGR) and click `Generate plot` button."),
											placement = "right", trigger = "hover", options = list(container = "body")
										)
										,tags$style(type='text/css', '#search_gene_bs {background-color: transparent; border:none;}')										
										,div(style="line-height:0px;",textInput('search_gene', '',width='300px' ))

										,verbatimTextOutput("gene_warn")
										,tags$style(type='text/css', '#gene_warn {color: red; font-weight: bold; font-size:18px;}')
										,verbatimTextOutput("gene_sym_warn")
										,tags$style(type='text/css', '#gene_sym_warn {color: red; font-weight: bold; font-size:18px;}')
										,verbatimTextOutput("gene_sym_green")
										,tags$style(type='text/css', '#gene_sym_green {color: green; font-weight: bold; font-size:18px;}')
										,plotlyOutput("detect_gene", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
								)
							)
							,tabPanel(div(style="display: inline-block;vertical-align:right;","Create a signature",bsButton("reg_sig_bs", label = "", icon = icon("info-circle"),style="", size = "small")),value="sig_data", br()
								,conditionalPanel("input.signature == 'sig_data'"
									,bsPopover(id = "reg_sig_bs", title = "",
										content = paste0("We use Bioconductor package <B>edgeR</B> to find differentially expressed genes between groups. Control group is considered as reference in computing log of fold change and p-values are adjusted for multiple comparison."),
										placement = "bottom", trigger = "hover", options = list(container = "body")
									)	
									,tags$style(type='text/css', '#reg_sig_bs {background-color: transparent; border:none;}')										
								), br()
								,shinyjs::hidden(verbatimTextOutput("sig1_warn"))
								,tags$style(type='text/css', '#sig1_warn {color: red; font-weight: bold; font-size:18px;}')
								,tabsetPanel(id="cr_reg_sig"
									,tabPanel("Metadata",value="meta_reg", br()
										,verbatimTextOutput("rank_warn")
										,tags$style(type='text/css', '#rank_warn {color: red; font-weight: bold; font-size:18px;}')
										,verbatimTextOutput("level_warn")
										,tags$style(type='text/css', '#level_warn {color: red; font-weight: bold; font-size:18px;}')
										,DT::dataTableOutput('sig_meta', width = 900)
									)
									,tabPanel("Signature",value="sig_reg", br()
										,DT::dataTableOutput('sig_data', width = 900) %>% withSpinner(type=5)
									)
								)
							)
							,tabPanel(div(style="display: inline-block;vertical-align:right;","User specific design", bsButton("user_sig_bs", label = "", icon = icon("info-circle"),style="", size = "small")),value="exp_desg", br()
								,conditionalPanel("input.signature == 'exp_desg'"
									,bsPopover(id = "user_sig_bs", title = "",
										content = paste0("You can modify the experimental design by placing samples to your group of interest."),
										placement = "bottom", trigger = "hover", options = list(container = "body")
									)	
									,tags$style(type='text/css', '#user_sig_bs {background-color: transparent; border:none;}')										
								), br()
								,shinyjs::hidden(verbatimTextOutput("sig2_warn"))
								,tags$style(type='text/css', '#sig2_warn {color: red; font-weight: bold; font-size:18px;}')
								,tabsetPanel(id="expr_design"
									,tabPanel("Metadata",value="sig_meta", br()
										,verbatimTextOutput("user_design_title")
										,tags$style(type='text/css', '#user_design_title {color: navyblue; font-size:18px; font-weight: bold;}')

										,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='sig_meta'"
											,tags$style(".container { border:1px solid steelblue; margin: 0 auto; font-size:7px; width: 100%; height: 270px; overflow-y: scroll; }")
											,column(6
												,wellPanel(
													uiOutput("user_g1_ui")
												)
											)
											,column(6
												,wellPanel(
													uiOutput("user_g2_ui")
												)
											)
										), tags$hr()
										,verbatimTextOutput("user_rank_warn")
										,tags$style(type='text/css', '#user_rank_warn {color: red; font-weight: bold; font-size:18px;}')
										,verbatimTextOutput("user_level_warn")
										,tags$style(type='text/css', '#user_level_warn {color: red; font-weight: bold; font-size:18px;}')
										,DT::dataTableOutput('designtable', width = 900)										
									)
									,tabPanel("Signature",value="user_exp_sig", br()
										,DT::dataTableOutput('sigtable', width = 900) %>% withSpinner(type=5)
									)
								)
							)
						)
					)
				)
			)
			,includeHTML("www/html/footer.html")
		)
		,tabPanel(tags$b('Processing console'), value='out_cons'
			,column(3, offset=0.5, 
				textInput('geo_current', h4('Currently processing')), br(),
				selectInput("geo_next", "Waiting to process",choices="", selected="")
			)
			,column(9, offset=0.5,
				tabsetPanel( id="process",
					tabPanel("Processing log",value="log", br()
						,verbatimTextOutput("steps")
						,actionButton("show_log","Show log file", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"), br()
						,br(), verbatimTextOutput('process_log_text')
					)
				)		
			)
			,includeHTML("www/html/footer.html")
		)
	)
)
