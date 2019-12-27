library(shiny)
library(org.At.tair.db) #ataliana
library(org.Ce.eg.db) #celegans
library(org.Dm.eg.db) #dmelanogaster
library(org.Dr.eg.db) #drerio
library(org.EcK12.eg.db) #ecoli straink12
library(org.Hs.eg.db) #hsapiens
library(org.Mm.eg.db) #mmusculus
library(org.Sc.sgd.db) #scerevisiae
library(org.Xl.eg.db) #xenopus
library(clusterProfiler)
library(ReactomePA) #enrichpathway
library(shinycssloaders) #spinner while plot is loading
library(shinyWidgets) #switch button dotplot or barplot
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Xl.eg.db")
# BiocManager::install("org.At.tair.db")
# BiocManager::install("org.Ce.eg.db")
# a
# BiocManager::install("org.Dm.eg.db")
# a
# BiocManager::install("org.Dr.eg.db")
# a
# BiocManager::install("org.EcK12.eg.db")
# a
# BiocManager::install("org.Hs.eg.db")
# a
# BiocManager::install("org.Mm.eg.db")
# a
# BiocManager::install("org.Sc.sgd.db")
# a

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Uploading Files"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      helpText("Load a genes list and extract some informations (...) "),
      # Input: Select a file ----
      fileInput("file", "Choose a CSV gene list File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      textAreaInput("manualEntry", "...Or write your list of genes", width = "250px"),
      helpText("Write a list of genes like this : geneID,geneID2,geneID3"),
      verbatimTextOutput("value"), #Will allow the attribution of manualEntry in var output$value

      # Only show the file options if there is a file uploaded 
      conditionalPanel(
        condition = "output.file_listener",
      # Horizontal line ----
      tags$hr(),
      h4("Uploaded file organization"),
      # Input: Checkbox if file has header ----
      checkboxInput("header", "File has header ?", FALSE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t",
                               WhiteSpace = " "),
                   selected = ","),
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"')

      ),
      # Horizontal line ----
      tags$hr(),
      selectInput("organism", "Select organism", list(
        "Vertebrate"  = c("Human (Homo sapiens)"="org.Hs.eg.db","Mouse (Mus musculus)"="org.Mm.eg.db",
                          "Western clawed frog (Xenopus tropicalis)"="org.Xl.eg.db","Zebrafish (Danio rerio)"="org.Dr.eg.db"
                          ),
        "Eukaryotic" = c("Fruit fly (Drosophila melanogaster)"="org.Dm.eg.db",
                          "Nematode worm (Caenorhabditis elegans)"="org.Ce.eg.db","Yeast (Saccharomyces cerevisiae)"="org.Sc.sgd.db"),
        "Prokaryotic" = c("Bacteria (Escherichia coli)"="org.EcK12.eg.db"),
        "Plantae" = c("Angiosperma (Arabidopsis thaliana)"="org.At.tair.db")), selected = NULL, multiple = FALSE),
 
      radioButtons("idType", "Gene identifier type",
                   choices = c(ENTREZ = "ENTREZID",
                               ENSEMBL = "ENSEMBL",
                               SYMBOL = "SYMBOL",
                               "COMMON NAME" = "COMMON"),
                   selected = "ENSEMBL"),
      
      checkboxGroupInput("ontology", "Ontologies",
                         c("Molecular Function" = "MF",
                           "Cellular Component" = "CC",
                           "Biological Process" = "BP")),
    # submit button
      actionButton("submit", "Submit"), width = 2),
    # Main panel for displaying outputs ----
      mainPanel(
      h2("This is the main title of the page"),
      tabsetPanel(id="tabs",
        tabPanel(value= "tab1","GeneOntology enrichment", fluidRow(
          #Deal the number of cellWidths with the number of checkbox, divide 100/number and replace
          splitLayout(cellWidths = c("33%", "33%", "33%"), plotOutput("barplot_BP"),switchInput(
            inputId = "switchBP",
            onLabel = "barplot",
            offLabel = "dotplot"
          ),
                      withSpinner(plotOutput("barplot_MF"),type = 5,color="#0dc5c1"), plotOutput("barplot_CC")))),
        tabPanel(value= "tab2","KEGG enrichment",
                    withSpinner(plotOutput("barplot_ekegDEGs"), type = 5,color="#0dc5c1"),
                    plotOutput("barplot_ekegMDGEs"),
                    
                 # )
        # )
      ),        
           
        
        tabPanel(value= "tab3","Msigdbr"),
        tabPanel(value= "tab4","Reactome",
                 withSpinner(plotOutput("barplot_ekePDEGs"),type = 5,color="#0dc5c1"),
                 plotOutput("dotplot_ekePDEGs"))
      ))
  )
)

####SERVER####################

# Define server logic to read selected file ----
server <- function(input, output) {
  output$file_listener <- reactive({
    return(!is.null(input$file))
  })
  
  #Listener is always active 
  #  trim <- function (x) gsub("^\\s+|\\s+$", "", x) #remove trailing and leading whitespace
  
  #Wait for the user to finish filling the informations and update all the variables
  observe({text_reactive()})
  text_reactive <- eventReactive( input$submit, {

   
 ### Retrieve list of genes ----
    geneList = reactive({
      #Check if the user wants to upload a file or paste a list of genes
      if(!is.null(input$file)){
        geneList=read.csv(input$file$datapath,
                          header = input$header,
                          sep = input$sep,
                          strip.white=TRUE,
                          quote = input$quote
        )
       
        
        geneList = as.character(geneList[,1])
        geneList = sort(geneList, decreasing = TRUE)
        # return(list('data'=geneList))
        return(geneList)
      }
      else{
        genes = c(strsplit(input$manualEntry, "[,\\t;  ]")) 
        #will recognize any character used to write the list
        return(genes)
      }
    })
    
  ### Isolate the reactive values to store them and manipulate them ----
  isolate({
    ###Format of the IDs ?
    req(geneList())
    all_types = c(ENTREZ = "ENTREZID",
                  ENSEMBL = "ENSEMBL",
                  SYMBOL = "SYMBOL")

    #Starting to observe if tabs are clicked
    observeEvent(input$tabs, {
      
    # x= geneList()$data
    gene_list= geneList()
    # observe(print(class(x)))
    # vecgenes= x[[1]]
    # names(vecgenes)=x[[1]]
    # gene_list<-na.omit(vecgenes)
    # gene_list = sort(gene_list, decreasing = TRUE)
    geneConv.df = bitr(gene_list, fromType = input$idType,
                   toType = c(all_types[!all_types %in% input$idType]),
                   OrgDb = input$organism)
    
    #/!\ input$idType might be different than ENSEMBL !!
    rownames(geneConv.df) <- geneConv.df[["ENSEMBL"]]
    observe(print(geneConv.df))
    #remove ensembl column
    geneConv.df$ENSEMBL <- NULL
    #add only symbol elements
    geneListSYMB <- geneConv.df[["SYMBOL"]]
    observe(print(geneListSYMB))
    observe(print(geneConv.df))
    #observe plotbutton or dotplot
    if(input$tabs == "tab1"){
    ### GO enrichments ------------
    # GO group enrichment
    observe(print(class(input$ontology)))
    for(i in 1: length(input$ontology)){
      if(input$ontology[i] == "MF"){
        ggoDEGs_MF <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
        output$barplot_MF <-renderPlot({
          barplot(ggoDEGs_MF, showCategory = 30,  title = "GroupGO DEGs MF")
        })

      }else if (input$ontology[i] == "BP"){
        ggoDEGs_BP <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
        output$barplot_BP <-renderPlot({
          barplot(ggoDEGs_BP, showCategory = 30,  title = "GroupGO DEGs BP")
        })
      }else if (input$ontology[i] == "CC"){
        ggoDEGs_CC <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
        output$barplot_CC <-renderPlot({
          barplot(ggoDEGs_CC, showCategory = 30,  title = "GroupGO DEGs CC")
        })
      }
      
    }
    }

    # # Create the "universe" gene set (this is a set of all the "detected" genes)
    # row_NON_zero <- apply(tpmPPglu, 1, function(row) all(row != 0))
    # tpmPPgluClean <- tpmPPglu[row_NON_zero,]
    # univPPglu <- rownames(tpmPPgluClean)
    # geneConv2 <- bitr(univPPglu, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
    # univPPglu2 <- geneConv2[["SYMBOL"]]
    # 
    # ## MSIGDB enrichment
    # m_df <- msigdbr(species = "Homo sapiens")
    # m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!
    # 
    # # Gene enrichment
    # esigDEGs <- enricher(genesSYMB, universe = univPPglu2, TERM2GENE = m_df, minGSSize = 50, maxGSSize = 1000)  # One can play arounf with the set sizes to obtain something meaningfull.
    # barplot(esigDEGs, showCategory = 50, title = "BarPlot DEGs MsiGDB enrichment")
    # dotplot(esigDEGs, showCategory = 50, title = "DotPlot DEGs MsiGDB enrichment")
    # 
    # # Gene set enrichment
    # esigsDEGs <- GSEA(geneListSYMB, minGSSize = 20, TERM2GENE = m_df)
    # dotplot(esigsDEGs, showCategory = 50, title = "DotPlot MsiGDB DEGS GSEA")
    
 ## KEGG pathway enrichment -----
    # Convert ENSEMBL IDs to ENTREZ <===== WHY NOT USE BITR ?
    # genesENTREZ <- as.character(mapIds(input$organism, gene_list, to ="ENTREZID", from= "ENSEMBL"))

    if(input$tabs == "tab2"){
      
    genesENTREZ= geneConv.df[["ENTREZID"]]
    allOrganisms_KEGG <- c("hsa","mmu","xtr","dre","dme","cel","sce","ecok","ath")
    names(allOrganisms_KEGG)= c("org.Hs.eg.db","org.Mm.eg.db","org.Xl.eg.db","org.Dr.eg.db","org.Dm.eg.db","org.Ce.eg.db","org.Sc.sgd.db",
    "org.EcK12.eg.db","org.At.tair.db")

    # Enrich KEGG pathways
    ekegDEGs <- enrichKEGG(gene = genesENTREZ, organism = allOrganisms_KEGG[[input$organism]], pvalueCutoff = 0.05)
    output$barplot_ekegDEGs <-renderPlot({
      barplot(ekegDEGs, title = "DEGs KEGG enrichment")  # Only one "ribosome"
    })
    # Enrich KEGG modules
    ekegMDGEs <- enrichMKEGG(gene = genesENTREZ, organism = allOrganisms_KEGG[[input$organism]], pvalueCutoff = 0.05)
    output$barplot_ekegMDGEs <-renderPlot({
      barplot(ekegMDGEs, title = "DEGs KEGG modules enrichment")  # Only one "ribosome"
    })
   
    }
    if(input$tabs == "tab4"){
    # Enrich REACTOME Pathways
    ###/!\ A.thaliana and E.coli are not in the list, here human and yeast were written but will give a result error
    allOrganisms_enrichKEGG = c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly","yeast","human")
    names(allOrganisms_enrichKEGG)= c("org.Hs.eg.db","org.Mm.eg.db","org.Xl.eg.db","org.Dr.eg.db","org.Dm.eg.db","org.Ce.eg.db","org.Sc.sgd.db", "org.EcK12.eg.db","org.At.tair.db")
 
    ekePDEGs <- enrichPathway(gene = genesENTREZ, organism = allOrganisms_enrichKEGG[[input$organism]], pvalueCutoff = 0.05)
    output$barplot_ekePDEGs <-renderPlot({
      barplot(ekePDEGs, showCategory = 30, title = "DEGs REACTOME Pathways enrichment")
    })
    output$dotplot_ekePDEGs <-renderPlot({
      dotplot(ekePDEGs, showCategory = 20, title = "DEGs REACTOME Pathways enrichment")
    })
    }
    })
  })
  })
}  

# Run the app ----
shinyApp(ui, server)