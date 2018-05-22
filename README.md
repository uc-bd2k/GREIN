![](../master/www/images/logo1c_icon.png)
---
# GREIN : GEO RNA-seq experiments interactive navigator for re-analyzing GEO RNA-seq data

---
GREIN is accessible at:   
https://shiny.ilincs.org/grein
---
The Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) is a public repository of gene expression data 
that hosts more than 6,000 RNA-seq datasets and this number is increasing. Most of these samples are deposited in raw sequencing 
format which needs to be downloaded and processed. With an aim to transform all these datasets in an analysis-ready format, 
we are currently processing the available RNA-seq samples of human, mouse, and rat from GEO using an R-based automated pipeline [GREP2](https://github.com/uc-bd2k/GREP2). 
This pipeline simultaneously downloads and processes RNA-seq raw sequencing data available in GEO. We demonstrate the results in a web 
application, GREIN (GEO RNA-seq Experiments Interactive Navigator) using the shiny framework. This interactive and intuitive application allows a user 
with little or no computational programming background to explore and analyze processed GEO RNA-seq datasets in a point-and-click manner. 
GREIN provides the flexibility to analyze and create ilincs compliant signatures that can be uploaded to ilincs for further in-depth analysis. 
In addition, GREIN produces publication-ready graphs and let the user to download all the analysis results. By accumulating the processed 
GEO RNA-seq datasets in a common platform, we present GREIN as a prominent choice to the practitioner for analyzing GEO RNA-seq datasets.

![GEO RNA-seq pipeline workflow](../master/www/images/About_steps2.png)

---
## If you use [GREIN](https://shiny.ilincs.org/grein), please cite:

* Al Mahi,N. *et al* (2018) GREIN: An interactive web platform for re-analyzing GEO RNA-seq data. bioRxiv 326223; doi: https://doi.org/10.1101/326223

---
## Docker instructions

### Installation of Docker

Ubuntu: follow [the instructions](https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/) to get Docker CE for Ubuntu.

Mac: follow [the instructions](https://store.docker.com/editions/community/docker-ce-desktop-mac) to install [`the Stable verion of Docker CE`](https://download.docker.com/mac/stable/Docker.dmg) on Mac.

Windows: follow [the instructions](https://docs.docker.com/toolbox/toolbox_install_windows/) to install [`Docker Toolbox`](https://download.docker.com/win/stable/DockerToolbox.exe) on Windows.

---
### Download and run the docker container
To obtain the latest docker image, run the following in your command line:
```
[sudo] docker pull ucbd2k/grein
```
Linux users may need to use `sudo` to run Docker.

To run the container execute the following command:

```
[sudo] docker run -d -p <an available port>:3838 ucbd2k/grein
```
Typically one can use port 3838 if not already used by another application. In that case the commad is

```
[sudo] docker run -d -p 3838:3838 ucbd2k/grein
```

First make sure that port 3838 is free to use. If not free, you can stop and kill any othe docker containers on this port by

```
[sudo] docker stop <container ID> && docker rm <container ID>
```
To know the container ID run this command:
```
docker ps -a
```

Before starting GREIN, please download the github repository first. We have put an example dataset (GSE22666) in the `data` folder.
You should keep all the processed datasets with corresponding files (eset.RData, multiqc_report.html, and transcripts.RData) in this folder. You can process GEO RNA-seq daatsets using our 
[GREP2](https://github.com/uc-bd2k/GREP2) pipeline. If you want to process data on the fly using GREIN, then you will have to run the GREP2 pipeline
in the backend (run `GEO_data_processing.R` for any new dataset you want to process and change the parameters accordingly) which will grab the GEO accession ID from
the [user_geo_request](https://github.com/uc-bd2k/GREIN/tree/master/data/user_geo_request) folder. GREIN will look for the log file within
this directory. After you process the data, you will have to run `after_processing.R` which will update the datatable used by GREIN.

To start GREIN, open a browser and type in the address bar ``<Host URL>:<available port as specified>``. For example `http://localhost:3838` on Mac or Linux systems when 3838 port is used.

Host URL on Ubuntu and Mac is `localhost`, if accessed locally. On Windows, the IP is shown when Docker is launched by double-clicking the Docker Quickstart Terminal icon on desktop, or it can be obtained from the output of `docker-machine ls` in the interactive shell window.

---
## Issues and bug reports

Please submit any bug reports, issues, or comments here: https://github.com/uc-bd2k/GREIN/issues