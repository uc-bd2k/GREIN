FROM rocker/shiny
COPY . /srv/shiny-server

# RUN apt-get update
# RUN apt-get install libmariadb-client-lgpl-dev
# RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biobase'); biocLite('ComplexHeatmap'); biocLite('circlize'); install.packages('shinyjs'); install.packages('shinyBS'); install.packages('RMySQL'); install.packages('reshape2'); "
RUN apt-get update
RUN apt-get install -y libssl-dev libssh2-1-dev
RUN apt-get install -y libmariadb-client-lgpl-dev
RUN apt-get install -y libxml2-dev libx11-dev
RUN apt-get install -y libglu1-mesa-dev freeglut3-dev mesa-common-dev

RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biobase'); biocLite('ComplexHeatmap'); biocLite('RNASeqPower'); biocLite('rols'); biocLite('rhdf5');"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('circlize'); biocLite('edgeR'); biocLite('DESeq2'); biocLite('limma');"
RUN R -e "install.packages(c('shinyjs', 'shinyBS', 'RMySQL', 'shinycssloaders'), repos = 'http://cran.us.r-project.org');"

RUN R -e "install.packages(c('devtools', 'rsconnect', 'httr', 'dplyr'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('shinyRGL', 'rgl', 'rglwidget', 'Rtsne'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('RColorBrewer', 'pairsD3'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('locfit', 'feather', 'data.table'), repos = 'http://cran.us.r-project.org');"
RUN R -e "install.packages(c('ggplot2', 'plotly', 'reshape2'), repos = 'http://cran.us.r-project.org');"
RUN R -e "devtools::install_github('ropensci/iheatmapr');"
RUN R -e "devtools::install_github('rstudio/DT');"

