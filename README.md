# ShinY-GL - Shiny Your Gene List

An R-Shiny application for the downstream analysis of gene list of interest

#### Developed by [Coralie Capron](https://github.com/Browco) and [Costas Bouyioukos](https://github.com/cbouyio)

_________

A simple to use R-Shiny application for the downstream gene enrichment and pathway analysis any gene list of interest.

The application takes a simple **un-ordered** gene list, or a gene list **ordered** by any arbritary criterion of interest (expression log-Fold-Change, enrichment in peaks, as so on) and performs the following analyses:

1. Gene Ontologies enrichment
2. Gene Set enrichments (in case an ordered list is provided)
3. REACTOME pathway enrichment
4. KEGG pathway enrichment

The application creates the relative visulasation of each analysis interms of dotplots and/or barplots.

### Requirements
Rstudio is absolutely needed to launch the app
Some R librairies are used in this application and need to be installed (it might take a few minutes to install them):
- shinyWidgets
- shinydashboard
- shiny
- clusterProfiler
- ReactomePA
- msigdbr
- shinycssloaders
- ggplot2
- shinyjs
- some Bioconductor databases (see below)
```
install.packages(c("shinyWidgets","shinydashboard","shiny","clusterProfiler","ReactomePA", "msigdbr", "shinycssloaders", "ggplot2", "shinyjs"))
```

You will need BiocManager to install Bioconductor databases 
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("org.At.tair.db", "org.Ce.eg.db","org.Dm.eg.db","org.Dr.eg.db","org.EcK12.eg.db","org.Hs.eg.db","org.Mm.eg.db","org.Sc.sgd.db", "org.Xl.eg.db"))
```

### How to launch the app ? 

Everything you need to achieve the app demonstration is on the shiny_app directory. 
You must open the file shiny_app/app.R in an Rstudio session. Then, click on the "Run app" button (upper right side of the code window).
Some tests files are provided : 
- test_gene_ensembl  : This is a 2 genes ID unordered file 
- test_gene_ensembl2col : This is a 5 genes ID ordered file used for GSEA  

### Versions of the different librairies used for the implementation of this project 
- R : V.3.6.1
- rstudio :  V.1.2.5019
- clusterProfiler : V.3.14.0
- ggplot2 : V.3.2.1
- ReactomePA : V.1.30.0
- msigdbr : V.7.0.1
- shiny : V.1.4.0
- shinycssloaders : V.0.2.0
- shinydashboard : V.0.7.1
- shinyjs : V.1.0
- shinyWidgets : V.0.5.0  

- all org.db are in V.3.10.0

### Support or Contact

Contact Costas Bouyioukos at: costas.bouyioukos@u-paris.fr
