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
library(shinycssloaders) #spinner while plot is loading

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
      tabsetPanel(
        tabPanel("Plot", fluidRow(
          splitLayout(cellWidths = c("33%", "33%", "33%"),  withSpinner(plotOutput("barplot_MF"), type = 5,color="#0dc5c1"),
                      withSpinner(plotOutput("barplot_CC"),type = 5,color="#0dc5c1"), withSpinner(plotOutput("barplot_BP"),type = 5,color="#0dc5c1")))),
        tabPanel("Summary", verbatimTextOutput("summary")), 
        tabPanel("Table", tableOutput("table"))
      ),
       textOutput("orga"),
      
      )
  )
)

####SERVER####################

# Define server logic to read selected file ----
server <- function(input, output) {
  output$file_listener <- reactive({
    return(!is.null(input$file))
  })
  
  #Listener is always active 
  outputOptions(output, "file_listener", suspendWhenHidden = FALSE)
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

    # ggoDEGs_BP2 <- groupGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", level = 2, keyType = "ENSEMBL", readable = TRUE)
    # barplot(ggoDEGs_BP2, showCategory = 30,  title = "GroupGO DEGs BP_2")
    #
    
  })
  })
}  
  
  

  
  # })



  
# Run the app ----
shinyApp(ui, server)