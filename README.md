# ShinY-GL - Shiny Your Gene List

An R-Shiny application for the downstream enrichment analyses of gene lists of interest

#### Developed by [Coralie Capron](https://github.com/Browco) and [Costas Bouyioukos](https://github.com/cbouyio)

_________
A simple to use R-Shiny application for the downstream enrichment analysis, including ontologies, pathways and gene sets, of any gene list of interest.

The application takes a simple **un-ordered** gene list, or a gene list **ordered** by any arbitrary criterion of interest (e.g. expression log-Fold-Change, enrichment in peaks, as so on) and performs the following analyses:

1. Gene Ontologies enrichment
2. Gene Set enrichments (in case an ordered list is provided)
3. REACTOME pathway enrichment
4. KEGG pathway enrichment

The application creates the relative visualisation of each analysis in terms of dotplots and/or barplots.

## Download
Current release [https://github.com/parisepigenetics/shiny-your-gene-list/releases/tag/firsLabRelease](https://github.com/parisepigenetics/shiny-your-gene-list/releases/tag/firsLabRelease)

## Requirements
The Rstudio suite is needed to launch the app. (please install a relatively new version of [RStudio](https://rstudio.com/products/rstudio/download/) as this will make easier the use and run of the app)
The following R libraries are required by this application and need to be installed before you try to run it:
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
You can install all requirements by running the following R command:
```
install.packages(c("shinyWidgets","shinydashboard","shiny","clusterProfiler","ReactomePA", "msigdbr", "shinycssloaders", "ggplot2", "shinyjs"))
```

To make sure you will have the databases for annotations of several model organisms install them with the following command:
(you will need BiocManager to install Bioconductor databases)
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("org.At.tair.db", "org.Ce.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db", "org.Xl.eg.db"))
```

## How to launch the app ?
Everything you need to make an app demonstration is on the shiny_app directory.
Open the file shiny_app/app.R in an Rstudio session. Then, click on the "Run app" button (upper right side of the code window) and then if everything is in place a window with the GUI should appear.

There are some test data provided to show the application's capabilites:
- test_gene_ensembl: This is a 2 genes ID unordered file
- test_gene_ensembl2col: This is a 5 genes ID ordered file used for GSEA  

### (for developers) versions of the different libraries used for the implementation of this project
- R: V.3.6.1
- rstudio:  V.1.2.5019
- clusterProfiler: V.3.14.0
- ggplot2: V.3.2.1
- ReactomePA: V.1.30.0
- msigdbr: V.7.0.1
- shiny: V.1.4.0
- shinycssloaders: V.0.2.0
- shinydashboard: V.0.7.1
- shinyjs: V.1.0
- shinyWidgets: V.0.5.0  

- all org.db are in V.3.10.0

## Support or Contact
Contact Costas Bouyioukos at: costas.bouyioukos@u-paris.fr
