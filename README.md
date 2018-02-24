# GRIN : GEO RNA-seq interactive navigator for Processing, Exploring, and Analyzing GEO RNA-seq Datasets

The Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) is a public repository of gene expression data 
that hosts more than 40,000 RNA-Seq samples and this number is increasing. Most of these samples are deposited in raw sequencing 
format which needs to be downloaded and processed. With an aim to transform all these datasets in an analysis-ready format, 
we are currently processing the available RNA-Seq samples of human, mouse, and rat from GEO using an R-based automated pipeline. 
This pipeline simultaneously downloads and processes RNA-Seq raw sequencing data available in GEO. We demonstrate the results in a web 
application, GRIN (GEO RNA-Seq Interactive Navigator) using the shiny framework. This interactive and intuitive application allows a user 
with little or no computational programming background to explore and analyze processed GEO RNA-Seq datasets in a point-and-click manner. 
GRIN provides the flexibility to analyze and create ilincs compliant signatures that can be uploaded to ilincs for further in-depth analysis. 
In addition, GRIN produces publication-ready graphs and let the user to download all the analysis results. By accumulating the processed 
GEO RNA-Seq datasets in a common platform, we present GRIN as a prominent choice to the practitioner for analyzing GEO RNA-Seq datasets.

![GEO RNA-seq pipeline workflow](../master/www/images/About_steps2.png)