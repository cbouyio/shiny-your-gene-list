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
      textAreaInput("manualEntry", "...Or write your list of genes", width = "200px"),
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
      # Horizontal line ----
      tags$hr(),
      radioButtons("idType", "Gene identifier type",
                   choices = c(ENTREZ = "ENTREZID",
                               ENSEMBL = "ENSEMBL",
                               SYMBOL = "SYMBOL",
                               "COMMON NAME" = "COMMON"),
                   selected = "ENSEMBL"),
      # submit button
      actionButton("submit", "Submit")),
    # Main panel for displaying outputs ----
      mainPanel(
      h3("This is the main title of the page"),
       textOutput("orga")
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
  
  observe({text_reactive()})
  text_reactive <- eventReactive( input$submit, {

   
    ####Retrieve list of genes
    geneList = reactive({
      #Check if the user wants to upload a file or paste a list of genes
      if(!is.null(input$file)){
        geneList=read.csv(input$file$datapath,
                          header = input$header,
                          sep = input$sep,
                          strip.white=TRUE,
                          quote = input$quote
        )
        # y = data.frame(x = names(read.csv(input$file$datapath,
        #                                        header = input$header,
        #                                        sep = input$sep, 
        #                                        strip.white=TRUE,
        #                                        quote = input$quote
        # )))
        # y[] <- lapply(y, as.character)
        # return(list('genes'=csv))
        
        # geneList = as.character(geneList[,1])
        # geneList = sort(geneList, decreasing = TRUE)
        return(list('data'=geneList))
      }
      else{
        genes = c(strsplit(input$manualEntry, "[,\\t;  ]")) 
        #will recognize any character used to write the list
        return(list('data'=genes))
      }
    })
    
  # Isolate the reactive values to store them and manipulate them
  isolate({
    ###Format of the IDs ?
    req(geneList())
    all_types = c(ENTREZ = "ENTREZID",
                  ENSEMBL = "ENSEMBL",
                  SYMBOL = "SYMBOL")


    observe(print("GeneList"))
    x= geneList()$data
    vecgenes= x[[1]]
    names(vecgenes)=x[[1]]
    gene_list<-na.omit(vecgenes)
    gene_list = sort(gene_list, decreasing = TRUE)
    observe(print(gene_list))
    observe(print(names(vecgenes)))

    gene.df = bitr(names(vecgenes), fromType = input$idType,
                   toType = c(all_types[!all_types %in% input$idType]),
                   OrgDb = input$organism)

    observe(print(gene.df))
    observe(print(summary(gene.df)))
  })
  })
}  
  
  
  ### GO enrichments ------------
  # GO group enrichment
  # ggoDEGs_MF2 <- groupGO(gene = genesENS, OrgDb = database(), ont = "MF", level = 2, keyType = "ENSEMBL", readable = TRUE)
  # barplot(ggoDEGs_MF2, showCategory = 30,  title = "GroupGO DEGs MF_2")
  # 
  # 
  # ggoDEGs_BP2 <- groupGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "BP", level = 2, keyType = "ENSEMBL", readable = TRUE)
  # barplot(ggoDEGs_BP2, showCategory = 30,  title = "GroupGO DEGs BP_2")
  #   
  
  
  # })



  
# Run the app ----
shinyApp(ui, server)