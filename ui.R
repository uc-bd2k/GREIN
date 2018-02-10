options(rgl.useNULL = TRUE)
library(shiny)	
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(Biobase)
library(DT)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(RMySQL)
require(rgl)
require(shinyRGL)
library(rsconnect)
require(pairsD3)
library(rglwidget)
library(edgeR)
library(locfit)
library(plotly)
library(ggplot2)
library(reshape2)
library(RNASeqPower)
library(httr)
library(DESeq2)
library(iheatmapr)
#options(shiny.sanitize.error=FALSE)
#options(shiny.trace=TRUE)
source("R/dropdownButton.R", local = TRUE)

shinyUI(
	navbarPage(
		inverse=TRUE, windowTitle = 'GRIN' 
		#title= "GRIN",
		,title= tags$a(href='https://www.ilincs.org',target="_blank", tags$img(src='images/iLINCSnewLogo.png', height="45px",width="55px", style='padding:0px 0px 13px 10px;'))
		#tags$img(src='images/grin_logo1.png', height="60px",width="40px", style='padding:0px 0px 24px 0px;')),
		#,title=img(src='images/iLINCSnewLogo.png', style="float:right; padding-right:25px", height="45px",width="55px")
		,id='grin', selected="datasets"
	
		,tabPanel('About', value='about',
			
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
				,tags$style(type="text/css",
					".shiny-output-error { visibility: hidden; }",
					".shiny-output-error:before { visibility: hidden; }")
				,tags$link(href="css/chardinjs.css",rel="stylesheet")
				,tags$link(href="css/ilincs.css",rel="stylesheet")
				#,includeHTML("www/html/nav.html")
			)
			,useShinyjs()
			,column(1
				
			)
			,column(10,
				wellPanel(
					includeMarkdown("www/description1.Rmd")
					,tags$img(src='images/About_steps.png',height="500px",width="800px", style='display: block; margin-left: auto; margin-right: auto;')
					,includeMarkdown("www/description2.Rmd")
				)
			)
			,column(1
				
			)
	 		,includeHTML("www/html/footer.html")
		),
		tabPanel('GEO datasets', value='datasets'
			
			,includeMarkdown("www/title.Rmd")
			,tags$hr()
			,column(4, offset=0, style='padding:0px;'
				, div(style="display: inline-block;vertical-align:right;", textInput('user_geo', h4('Search for GEO series (GSE) accession'), width='auto'))
				, div(style="display: inline-block;vertical-align:right;", bsButton("q1", label = "", icon = icon("question"),style="info", size = "extra-small"))
				, bsPopover(id = "q1", title = "",
						content = paste0("Please type your GEO accession id (e.g., GSE46597) in this box to see if the dataset is already processed. If not, you can initialize the processing of the data. The dataset will be processed and uploaded soon."),
						placement = "right", trigger = "focus", options = list(container = "body")
					)

				, verbatimTextOutput("warn")
				, verbatimTextOutput("warn2")
				, shinyjs::hidden(verbatimTextOutput("warn3"))
				, verbatimTextOutput("warn4")
				, tags$style(type='text/css', '#warn {color: red; font-size:15px;}') #font-weight: bold;
				, tags$style(type='text/css', '#warn2 {color: red; font-size:15px;}')
				, tags$style(type='text/css', '#warn3 {color: red; font-size:15px;}')
				, tags$style(type='text/css', '#warn4 {color: red; font-size:15px;}')
				
				, div(style="display: ;vertical-align:bottom;", shinyjs::hidden(actionButton("start_process","Start processing", icon("refresh"),
						style="color: #fff; background-color: #337ab7; border-color: #2e6da4")))
				, tags$style(type='text/css', "label {font-size: 18px;}")
				, br(), selectInput("geo_next_lp", "User requested datasets in the processing queue", choices="", selected="")
				, actionButton("process_log","Processing console", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
			)
			,column(4, offset=0, style='padding:0px;',
				plotlyOutput("study_stat", width = "auto", height = "300px"), br()
			)
			,column(4,
				plotlyOutput("sample_stat", width = "auto", height = "300px"), br()
			)
			,br(), div(style="color:#337ab7;", titlePanel("Select a dataset"))
			#,h3('Select a dataset'), br()
			, br(), DT::dataTableOutput('datatable'), br()
			,includeHTML("www/html/footer.html")
		),
		tabPanel('Explore dataset', value='explore',
			column(3, offset=0.5, 
				wellPanel(
					#conditionalPanel("input.geo_acc != ''"
						#,uiOutput("geo_acc_ui")
						textInput('geo_acc', h3('Selected study'))
					,conditionalPanel("input.geo_acc != ''"
						,actionButton("analysisbtn", "Analyze", icon("th"), 
							style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
						,br(), br(), tags$div(class="metadata_column_class", uiOutput("metadata_column"))
						,tags$head(tags$style(HTML(".metadata_column_class label{font-size: 16px;}")))
						,br(), tags$div(class="full_metacol_class", uiOutput("full_metacol"))
						,tags$head(tags$style(HTML(".full_metacol_class label{font-size: 16px;}")))
						,br(), shinyjs::hidden(downloadButton('downloadmeta', label = "Download metadata", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
						#,shinyjs::hidden(downloadButton('downloadcounts', label = "Download gene level counts data", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
						,downloadButton('downloadcounts', label = "Download counts table", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
						,br(), br(), div(style="font-size: 2px;display: inline-block;vertical-align:right;", radioButtons('download_type_counts', label = "", choices = c("Gene level", "Transcript level"), inline = T))
						,tags$div(class="heat_prop_class", uiOutput("heat_prop"))
						,tags$head(tags$style(HTML(".heat_prop_class label{font-size: 16px;}")))
						,uiOutput("cluster_meth")
						,uiOutput("cluster_genes")
						,br(), div(style="display: inline-block;vertical-align:right;", downloadButton('downloadFit', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
						,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit', label = "", choices = c("png", "csv"), inline = T))
						,div(style="display: inline-block;vertical-align:right;", downloadButton('downloadScroll', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
						,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scroll', label = "", choices = c("png", "csv"), inline = T))
						,uiOutput("pca_prop")
						#,tags$style(".container_pca { border:1px solid steelblue; margin: 0 auto; font-size:7px; width=200%; height: 200x; overflow-y: scroll; overflow-x: scroll; }")
						,uiOutput("legend")	
						#,br(), uiOutput("heatmaptext")
						#,br(), uiOutput("metadatatext")
					)
				)
			),
			column(9, offset=0.5,
				wellPanel(
					conditionalPanel("input.geo_acc != ''"
						,tabsetPanel( id="tab2"
							,tabPanel("Description", value="summary", br()
								,htmlOutput("geo_summary") %>% withSpinner(type=5)
							)
							,tabPanel("Metadata",value="metadata", br()
								,wellPanel( uiOutput("metadatatext")), br()
								,DT::dataTableOutput('metadata') %>% withSpinner(type=5)
								,DT::dataTableOutput('metadata_full', width = 900)
							)
							,tabPanel("Counts table",value="countsTable", br()
								,DT::dataTableOutput('counts', width = 900) %>% withSpinner(type=5)
							)
							,tabPanel("QC report",value="qc", br()
								,downloadButton('downloadqc', label = "Download QC report", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
								,htmlOutput('multiqc') %>% withSpinner(type=5)
							)
							,tabPanel("Visualization",value="plots"
								,tabsetPanel(id="preplots"
									,tabPanel("Correlation plot",value="corr", br()
										,plotlyOutput("corplot") %>% withSpinner(type=5)
									)
									,tabPanel("Density plot",value="dens", br()
										,plotlyOutput("densityp", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
									,tabPanel("Heatmap",value="heatmap"
										,br(), wellPanel(uiOutput("heatmaptext"))
										,tabsetPanel(id="heatvis"
											,tabPanel("Fit in screen",value="fitscrn", br()
												#,uiOutput("heatmaptext"), br()
												,uiOutput("fitscreen", width="auto", height="auto"),br(), br()
											)
											,tabPanel("Scrollable",value="scrolb", br() 
												,uiOutput("scroll", width="auto", height="auto"), br(), br()
											)
											,tabPanel("Interactive heatmap",value="iheat", br() 
												,uiOutput("iheatmapui", width = "auto", height = "auto")
											)
										)
									),
									tabPanel("Principal component analysis",value="pca"
										,tabsetPanel(
											tabPanel("Intercative 2-D PCA plot",br() 											
												,downloadButton('downloadPCA2D', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
												,uiOutput("pca2Dplot")
												,br()
											)
											,tabPanel("Intercative 3-D PCA plot", br() 
												,downloadButton('downloadPCA3D', label = "Download image", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
												,rglwidgetOutput("pca3Dplot", width = "auto", height = "700px") %>% withSpinner(type=5)
												,br()
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
		,tabPanel('Analyze dataset', value='analyze'
			,conditionalPanel("input.signature == 'sig_data'"
				,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h3("Construct differential expression signature for the selected dataset")))
			)
			,conditionalPanel("input.signature == 'exp_desg'"
				,div(style="margin-top: -20px;margin-bottom: 20px; color:#337ab7;", titlePanel(h3("Construct differential expression signature for user modified design")))
			)
			,column(3, offset=0.5, 
				wellPanel(
						textInput('geo_acc2', h3('Selected study'))
						#,uiOutput("geo_acc2_ui")
					,conditionalPanel("input.geo_acc != ''"
						,uiOutput("user_meta_prop")
						,uiOutput("user_ana_type")
						,uiOutput("user_cov_ui")
						
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='sig_meta'"
							,br(), uiOutput("gen_sig_btn_ui"), br()
						)
						# ,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='user_exp_sig'"
							# ,br(), uiOutput("user_dynamiclink")
						# )
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='user_exp_sig'"
							,br(), actionButton("heatmap_mod", "Show heatmap", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
							,bsModal("modalheatmap", "Heatmap of top DE genes (ranked by adjusted p-values)", "heatmap_mod", size = "large"
								,column(4, offset=-0.5
									,wellPanel(
										tags$div(class="heat_type_mod_class", uiOutput("heat_typeui_mod"))
										,tags$head(tags$style(HTML(".heat_type_mod_class label{font-size: 15px;}")))
										,tags$div(class="heat_prop_mod_class",uiOutput("heat_prop_mod"))
										,tags$head(tags$style(HTML(".heat_prop_mod_class label{font-size: 15px;}")))
										,uiOutput("cluster_meth_mod")
										,uiOutput("cluster_genes_mod"), br()
										,div(style="display: inline-block;vertical-align:right;", downloadButton('downloadFit_mod', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
										,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit_mod', label = "", choices = c("png", "csv"), inline = T))
										,div(style="display: inline-block;vertical-align:right;", downloadButton('downloadScroll_mod', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
										,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scroll_mod', label = "", choices = c("png", "csv"), inline = T))
									)
								)
								,column(8, offset=-0.5
									,wellPanel(
										tabsetPanel(id="heatmodal",
											tabPanel("Fit in screen",value="fitscrn_mod", br()
												,uiOutput("fitscreen_mod", width="auto", height="auto") 
											),
											tabPanel("Scrollable",value="scrolb_mod", br() 
												,uiOutput("scroll_mod", width="auto", height="auto")
											)
										)
									)
								)
							)
						)
						,uiOutput("metadata_analysis")
						,uiOutput("ana_type")
						,uiOutput("group1")
						,uiOutput("group2")
						,uiOutput("cov")
						#,br(), uiOutput("dynamiclink")
						,conditionalPanel("input.signature == 'sig_data'"
							,br(),actionButton("heatmap_mod_regsig", "Show heatmap", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
							,bsModal("modalheatmap_regsig", "Heatmap of top DE genes (ranked by adjusted p-values)", "heatmap_mod_regsig", size = "large"
								,column(4, offset=-0.5
									,wellPanel(
										tags$div(class="heat_typeui_mod_regsig_class",uiOutput("heat_typeui_mod_regsig"))
										,tags$head(tags$style(HTML(".heat_typeui_mod_regsig_class label{font-size: 15px;}")))
										,tags$div(class="heat_prop_mod_regsig_class",uiOutput("heat_prop_mod_regsig"))
										,tags$head(tags$style(HTML(".heat_prop_mod_regsig_class label{font-size: 15px;}")))
										,uiOutput("cluster_meth_mod_regsig")
										,uiOutput("cluster_genes_mod_regsig"), br()
										,div(style="display: inline-block;vertical-align:right;", downloadButton('downloadFit_mod_regsig', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
										,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_fit_mod_regsig', label = "", choices = c("png", "csv"), inline = T))
										,div(style="display: inline-block;vertical-align:right;", downloadButton('downloadScroll_mod_regsig', label = "Download", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
										,div(style="display: inline-block;vertical-align:right;", radioButtons('download_type_scroll_mod_regsig', label = "", choices = c("png", "csv"), inline = T))
									)
								)
								,column(8, offset=-0.5
									,wellPanel(
										tabsetPanel(id="heatmodal_regsig",
											tabPanel("Fit in screen",value="fitscrn_mod_regsig", br()
												,uiOutput("fitscreen_mod_regsig", width="auto", height="auto") 
											),
											tabPanel("Scrollable",value="scrolb_mod_regsig", br() 
												,uiOutput("scroll_mod_regsig", width="auto", height="auto")
											)
										)
									)
								)
							)
						)
						,br(), shinyjs::hidden(downloadButton('downloadsigdata', label = "Download signatures", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;"))
						,br(), br(), uiOutput("dynamiclink")
						#,br(), uiOutput("sig_method")
						,conditionalPanel("input.signature == 'exp_desg' && input.expr_design=='user_exp_sig'"
							,shinyjs::hidden(downloadButton('download_user_sigdata', label = "Download signatures", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;")), br()
							,br(), br(), uiOutput("user_dynamiclink")
						)
						
						,uiOutput("power_ui")
						,uiOutput("samples")
						,uiOutput("alpha")
						,uiOutput("effect")						
						,uiOutput("depth_ui")
						,uiOutput("search_gene_ui")
						#,br(), uiOutput("gene_detect_ui")						
						#,br(), uiOutput("pow_method")
					)
				)
			)
			,column(9, offset=0.5,
				wellPanel(
					conditionalPanel("input.geo_acc !='' && input.geo_acc2 !=''"
						,tabsetPanel( id="signature"
							,tabPanel("Power analysis",value="power"
								,tabsetPanel( id="sigtab"
									,tabPanel("Power curve",value="power_curve", br()
										,wellPanel(uiOutput("pow_method"))
										,plotlyOutput("powerplot", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
									,tabPanel("Detectability of genes",value="detec_gene", br()
										,wellPanel(uiOutput("gene_detect_ui"))
										,plotlyOutput("detect_gene", width = "auto", height = "auto") %>% withSpinner(type=5)
									)
								)
							)
							,tabPanel("Create a signature",value="sig_data", br()
								,wellPanel(uiOutput("sig_method"))
								,verbatimTextOutput("rank_warn") 
								,tags$style(type='text/css', '#rank_warn {color: red; font-weight: bold; font-size:18px;}')
								,DT::dataTableOutput('sig_data', width = 900) %>% withSpinner(type=5)
								#,br(), shinyjs::hidden(downloadButton('downloadsigdata', label = "Download signatures", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;")), br()
							)
							,tabPanel("User specific design",value="exp_desg"
								,tabsetPanel(id="expr_design"
									,tabPanel("Metadata",value="sig_meta", br()
										,wellPanel(uiOutput("user_method"))
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
										# ,conditionalPanel(
											# condition="$('html').hasClass('shiny-busy')",
											# HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
											# <span id=\"bar-text\">Generating signatures...</span></div></div>")
										# )
										,DT::dataTableOutput('sigtable', width = 900) %>% withSpinner(type=5)
										#,br(), shinyjs::hidden(downloadButton('download_user_sigdata', label = "Download signatures", style="color: #fff; background-color: #337ab7; border-color: #2e6da4;"))
									)
								)
							)
						)
					)
				)
			)
			,includeHTML("www/html/footer.html")
		)
		,tabPanel('Processing console', value='out_cons'
			,column(3, offset=0.5, 
				textInput('geo_current', h4('Currently processing')), br(),
				selectInput("geo_next", "Waiting to process",choices="", selected="")
			)
			,column(9, offset=0.5,
				tabsetPanel( id="process",
					tabPanel("Processing log",value="log", br()
						,verbatimTextOutput("steps")
						,verbatimTextOutput('process_log')
					)
				)		
			)
			,includeHTML("www/html/footer.html")
		)
	)
)
