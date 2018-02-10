FROM shiny.ilincs.org:5000/shiny:3.3.1-1.5.1
#FROM shiny.ilincs.org:5000/shiny
COPY . /srv/shiny-server

# RUN apt-get update
# RUN apt-get install libmariadb-client-lgpl-dev
# RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biobase'); biocLite('ComplexHeatmap'); biocLite('circlize'); install.packages('shinyjs'); install.packages('RMySQL'); install.packages('reshape2'); "
RUN apt-get update
RUN apt-get install -y libmariadb-client-lgpl-dev

RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biobase'); biocLite('ComplexHeatmap'); biocLite('RNASeqPower'); "
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('circlize'); biocLite('edgeR'); "
RUN R -e "install.packages(c('shinyjs', 'RMySQL'), repos = 'http://cran.us.r-project.org');"

RUN apt-get install -y libxml2-dev libx11-dev
RUN apt-get install -y libglu1-mesa-dev freeglut3-dev mesa-common-dev
RUN R -e "install.packages(c('devtools', 'rsconnect'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('shinyRGL', 'rgl', 'rglwidget'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('RColorBrewer', 'pairsD3'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('locfit', 'DT'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('ggplot2', 'plotly', 'reshape2'), repos = 'http://cran.us.r-project.org');"

