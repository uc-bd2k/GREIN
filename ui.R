options(rgl.useNULL = TRUE)
library(shiny)	
library(shinyjs)
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
options(shiny.sanitize.error=FALSE)

shinyUI(navbarPage(
	title = 'GRSN', id='explorgeo',
	
	tabPanel('GEO Datasets', value='datasets',
		
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
			,includeHTML("www/html/nav.html")
		),
		useShinyjs(),
		includeMarkdown("www/description.Rmd"),
 
			tags$hr(),
			h3('Select a dataset'), br(),
			conditionalPanel(
				condition="$('html').hasClass('shiny-busy')",
				HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
				<span id=\"bar-text\">Loading datasets...</span></div></div>")
			),
			DT::dataTableOutput('datatable'), br()
		,includeHTML("www/html/footer.html")
	),
	
	tabPanel('Explore dataset', value='explore',
		column(3, offset=0.5, 
			wellPanel(
				conditionalPanel("input.geo_acc != ''",
					textInput('geo_acc', h3('Selected study'))
					,actionButton("analysisbtn", "Analyze", icon("paper-plane"), 
						style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
					,uiOutput("metadata_column"), br()
					,uiOutput("full_metacol")
					,uiOutput("heat_prop")
					,selectInput("clustering_method", label = h4("Grouping samples"),
						choices = list("Group by properties" , "Pearson Correlation", "Euclidean distance"), selected="Pearson Correlation")
					,uiOutput("cluster_genes")
					,uiOutput("pca_prop")
					,uiOutput("legend")	
					, br(),uiOutput("countstext")
					, br(), uiOutput("heatmaptext")
					, br(), uiOutput("metadatatext")
				)
			)
		),
	
		column(9, offset=0.5,
			wellPanel(
				conditionalPanel("input.geo_acc != ''",
					tabsetPanel( id="tab2",
						tabPanel("Description", value="summary", br(), 
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Loading dataset...</span></div></div>")
							),
							htmlOutput("geo_summary")
						),
						tabPanel("Metadata",value="metadata", br(),
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Loading metadata...</span></div></div>")
							),
							DT::dataTableOutput('metadata', width = 900),br(),
							downloadButton('downloadmeta', label = "Download metadata")
						),
						tabPanel("Counts table",value="countsTable", br(),
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Loading expression data...</span></div></div>")
							),
							DT::dataTableOutput('counts'),br(),
							downloadButton('downloadcounts', label = "Download full counts table")
						),
						tabPanel("QC report",value="qc", br(), 
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Loading QC report...</span></div></div>")
							),
							downloadButton('downloadqc', label = "Download QC report"),
							htmlOutput('multiqc')
						),
						tabPanel("Precomputed plots",value="plots",
							tabsetPanel(id="preplots",
								tabPanel("Correlation plot",value="corr", br(),
									conditionalPanel(
										condition="$('html').hasClass('shiny-busy')",
										HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
										<span id=\"bar-text\">Loading correlation plot...</span></div></div>")
									),
									plotlyOutput("corplot")
								),
								tabPanel("Density plot",value="dens", br(),
									conditionalPanel(
										condition="$('html').hasClass('shiny-busy')",
										HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
										<span id=\"bar-text\">Loading density plot...</span></div></div>")
									),
									plotlyOutput("densityp", width = "auto", height = "auto")
								),
								tabPanel("Heatmap",value="heatmap",
									tabsetPanel(
										tabPanel("Fit in Screen",br(), 
											downloadButton('downloadFit', label = "Download image"), 
											conditionalPanel(
												condition="$('html').hasClass('shiny-busy')",
												HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
												<span id=\"bar-text\">Loading heatmap...</span></div></div>")),
											uiOutput("fitscreen", width="auto", height="auto"),br(), br()
										),
										tabPanel("Scrollable", br(), 
											downloadButton('downloadScroll', label = "Download image"),
											conditionalPanel(
												condition="$('html').hasClass('shiny-busy')",
												HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
												<span id=\"bar-text\">Loading heatmap...</span></div></div>")),
											plotOutput("scroll", width="auto", height="auto"), br(), br()
										)
									)
								),
								tabPanel("Principal component analysis",value="pca",
									tabsetPanel(
										tabPanel("Intercative 2-D PCA plot",br(), 
											conditionalPanel(
												condition="$('html').hasClass('shiny-busy')",
												HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
												<span id=\"bar-text\">Loading PCA plot...</span></div></div>")
											),
										downloadButton('downloadPCA2D', label = "Download image"),
										uiOutput("pca2Dplot"), 
										br()),
										tabPanel("Intercative 3-D PCA plot", br(), 
											conditionalPanel(
												condition="$('html').hasClass('shiny-busy')",
												HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
												<span id=\"bar-text\">Loading PCA plot...</span></div></div>")
											),
										downloadButton('downloadPCA3D', label = "Download image"),
										rglwidgetOutput("pca3Dplot", width = "auto", height = "700px"), 
										br())
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
	,tabPanel('Analyze dataset', value='analyze',

		column(3, offset=0.5, 
			wellPanel(
				conditionalPanel("input.geo_acc != ''"
					,textInput('geo_acc2', h3('Selected study')), br()
					,uiOutput("metadata_analysis")
					,uiOutput("group1")
					,uiOutput("group2")
					, br(), uiOutput("sig_method")
					, uiOutput("alpha")
					, uiOutput("effect")
					, uiOutput("samples")
				)
			)
		),
	
		column(9, offset=0.5,
			wellPanel(
				conditionalPanel("input.input_group2 !='' && input.geo_acc2 !='' ",
					tabsetPanel( id="signature",
						tabPanel("Signatures",value="sig_data", br(),
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Calculating signatures...</span></div></div>")),
							DT::dataTableOutput('sig_data'),br(),
							downloadButton('downloadsigdata', label = "Download signatures")
						),
						tabPanel("Power analysis",value="power", br(),
							conditionalPanel(
								condition="$('html').hasClass('shiny-busy')",
								HTML("<div class=\"progress\" style=\"height:25px !important\"><div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"40\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"width:100%\">
								<span id=\"bar-text\">Loading power curve...</span></div></div>")),
							plotlyOutput("powerplot", width = "auto", height = "auto")
						)
					)
				)
			)
		)
		,includeHTML("www/html/footer.html")
	)
			
)
)