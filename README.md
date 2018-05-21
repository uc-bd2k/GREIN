![](../master/www/images/logo1c_icon.png)
---
# GREIN : GEO RNA-seq experiments interactive navigator for re-analyzing GEO RNA-seq data

---
GREIN is accessible at:   
https://shiny.ilincs.org/grein
---
The Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) is a public repository of gene expression data 
that hosts more than 6,000 RNA-Seq datasets and this number is increasing. Most of these samples are deposited in raw sequencing 
format which needs to be downloaded and processed. With an aim to transform all these datasets in an analysis-ready format, 
we are currently processing the available RNA-Seq samples of human, mouse, and rat from GEO using an R-based automated pipeline [GREP2](https://github.com/uc-bd2k/GREP2). 
This pipeline simultaneously downloads and processes RNA-Seq raw sequencing data available in GEO. We demonstrate the results in a web 
application, GREIN (GEO RNA-Seq Experiments Interactive Navigator) using the shiny framework. This interactive and intuitive application allows a user 
with little or no computational programming background to explore and analyze processed GEO RNA-Seq datasets in a point-and-click manner. 
GREIN provides the flexibility to analyze and create ilincs compliant signatures that can be uploaded to ilincs for further in-depth analysis. 
In addition, GREIN produces publication-ready graphs and let the user to download all the analysis results. By accumulating the processed 
GEO RNA-Seq datasets in a common platform, we present GREIN as a prominent choice to the practitioner for analyzing GEO RNA-Seq datasets.

![GEO RNA-seq pipeline workflow](../master/www/images/About_steps2.png)

---
If you use [GREIN](https://shiny.ilincs.org/grein), please cite:

* Al Mahi,N. *et al* (2018) GREIN: An interactive web platform for re-analyzing GEO RNA-seq data. bioRxiv 326223; doi: https://doi.org/10.1101/326223

---
## Issues and bug reports

Please submit any bug reports, issues, or comments here: https://github.com/uc-bd2k/GREIN/issues