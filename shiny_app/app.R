
# IMPORT ---------
library(shiny)
library(shinydashboard)
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
library(ggplot2) #used for cnet plot
library(shinyjs) #show/hide barplot/dotplot in GOenrichment
# UI ---------------------------------------------------------------- 
  ## header ----
header <- dashboardHeader(
  title = "Shiny your gene list"
)
  ## sidebar ---- 
sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs",
  menuItem(text = "Importing file", tabName = "datafile",icon = icon("file")
  ),
  menuItem("Enrichment", tabName = "enrichoose", icon = icon("dna")),
  menuItem("Analysis", tabName = "analysis", icon = icon("clipboard-check"))
)
)
  ## body ----
body <- dashboardBody(
  fluidPage(
    useShinyjs(), #set up shinyjs, useful to show or hide barplot/dotplot
    tabItems(
      tabItem(tabName = "datafile",
              fluidRow(
                box(width = 12,
                    fileInput(inputId = "file",
                              label = "Choose a file",
                              accept = c(".csv")
                    ),
                    tableOutput(outputId = "Contents"),
                    verbatimTextOutput(outputId = "Data")
                ),
                box(width = 12,
                    textAreaInput(inputId="manualEntry",
                                  label = "...Or write your list of genes",
                                  value = NULL ,
                                  width = "250px",
                                  resize = "both", placeholder ="geneID1\ngeneID2 ",),
                   
                    helpText("Write an ordered or unordered list of genes lines by lines with the same separator"),
                    verbatimTextOutput("Manual"), #Will allow the attribution of manualEntry in var output$value)
              ))),
      tabItem(tabName = "enrichoose",
              fluidRow(
                box(width = 12,
                    selectInput("organism", "Select organism", 
                        list("Vertebrate"  = c("Human (Homo sapiens)"="org.Hs.eg.db","Mouse (Mus musculus)"="org.Mm.eg.db",
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
                    actionButton("submit", "Submit"))
                )),
      tabItem(tabName = "analysis",
              fluidRow(
                tabBox(width=12,title = "Enrichment analysis",
                  # The id lets us use input$tabset1 on the server to find the current tab
                  id = "tabset",
                  tabPanel(value= "tab1","GeneOntology enrichment", fluidRow( 
                    
                    box(width = 12, title = "Biological Process", collapsible = TRUE, status = "primary",
                        plotOutput("barplot_BP"), plotOutput("dotplot_BP"),switchInput(
                    inputId = "switchBP",
                    value = TRUE, #true is barplot
                    onLabel = "barplot",
                    offLabel = "dotplot"),plotOutput("cnetBP")),
                    
                    box(width = 12, title ="Molecular Function", collapsible = TRUE, status = "primary",
                       plotOutput("barplot_MF"), plotOutput("dotplot_MF"), switchInput(
                         inputId = "switchMF",
                         value = TRUE, #true is barplot
                         onLabel = "barplot",
                         offLabel = "dotplot"),plotOutput("cnetMF")),
                    
                    box(width = 12, title = "Cellular Component", collapsible = TRUE, status = "primary",
                        plotOutput("barplot_CC"), plotOutput("dotplot_CC"),switchInput(
                          inputId = "switchCC",
                          value = TRUE, #true is barplot
                          onLabel = "barplot",
                          offLabel = "dotplot"),plotOutput("cnetCC"))
                  )),
                  tabPanel(value= "tab2","KEGG enrichment",
                           withSpinner(plotOutput("barplot_ekegDEGs"), type = 5,color="#0dc5c1"),
                           plotOutput("barplot_ekegMDGEs"),
                  ),        
                  
                  
                  tabPanel(value= "tab3","Msigdbr"),
                  tabPanel(value= "tab4","Reactome",
                           withSpinner(plotOutput("barplot_ekePDEGs"),type = 5,color="#0dc5c1"),
                           plotOutput("dotplot_ekePDEGs"))
                
      ))
      )
    )
  )
)
              

ui <- dashboardPage(header = header, 
                    sidebar = sidebar, 
                    body = body)

# SERVER ----------------------------------------------------------------

server <- function(input, output, session) {
    data(geneList, package="DOSE") #used here for the Universe
  ### Retrieve list of genes given by user  ----
  geneList = reactive({
    if(is.null(input$file)){
      req(input$manualEntry) #check if something is written before opening it
      df <- read.csv(text=input$manualEntry,header=FALSE)
      print(df)
      geneList = as.character(df[,1])
    } else {
      inFile <- input$file
      df <- read.csv(inFile$datapath,header=FALSE)
      print(df)
      geneList = as.character(df[,1])
      # geneLog = as.character(df[,2])
      geneList = sort(geneList, decreasing = TRUE)
      return(geneList)
    }
  })
  #  observe(print(geneList())) #print value in terminal

  #Show a message window to tell the user to go to the next step when gene list is given 
  observeEvent({ 
    input$file
    input$manualEntry
    1
  }, 
  {     if (!is.null(input$file) & !is.null(input$manualEntry))
    showNotification(paste("Now choose your enrichment in the sidebar"), duration = 5, type = "message") } )

  observe({text_reactive()}) #text_reactive is always updated
  text_reactive <- eventReactive( input$submit, {
    ### Isolate the reactive values to store them and manipulate them ----
    isolate({
      ###Format of the IDs ?
      req(geneList()) #if exists
      all_types = c(ENTREZ = "ENTREZID",
                    ENSEMBL = "ENSEMBL",
                    SYMBOL = "SYMBOL")
      
      #Starting to observe if tabs are clicked...
      observeEvent(input$tabset,{
        
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
        genesENTREZ= geneConv.df[["ENTREZID"]]

#        ...And execute special command for each type of enrichment/tab
        ### Geneontology enrichment ----
        if(input$tabset == "tab1"){
          ### GO enrichments ------------
          observe(print(class(input$ontology)))
          for(i in 1: length(input$ontology)){
            #For each radiobutton selected during the "choose your enrichment" step, plots will be created
            if(input$ontology[i] == "MF"){
              ggoDEGs_MF <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_MF <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneList), keyType = "ENSEMBL", readable = TRUE)
              
              #All the different outputs (cnetplot, barplot and dotplot)
              output$cnetMF = renderPlot({
                cnetplot(egoDEGs_MF, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs MF")
              })
              output$barplot_MF <-renderPlot({
                barplot(ggoDEGs_MF, showCategory = 30,  title = "GroupGO DEGs MF")})
              output$dotplot_MF <-renderPlot({
                dotplot(egoDEGs_MF, title = "GO enrichment DEGs MF")
              })
              observeEvent(input$switchMF,{  #Show the barplot if switch button is on Barplot mode
              if (input$switchMF == TRUE){
                show("barplot_MF")
                hide("dotplot_MF")  
              }
              else{
                show("dotplot_MF")
                hide("barplot_MF")  
                
              }})
              
            }else if (input$ontology[i] == "BP"){
              ggoDEGs_BP <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_BP <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneList), keyType = "ENSEMBL", readable = TRUE)
              output$cnetBP = renderPlot({
                cnetplot(egoDEGs_BP, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs BP")
              })
              output$barplot_BP <-renderPlot({
                barplot(ggoDEGs_BP, showCategory = 30,  title = "GroupGO DEGs BP")
              })
              output$dotplot_BP <-renderPlot({
                dotplot(egoDEGs_BP, title = "GO enrichment DEGs BP")
              })
              
              observeEvent(input$switchBP,{  #Show the barplot if switch button is on Barplot mode
                if (input$switchBP == TRUE){
                  show("barplot_BP")
                  hide("dotplot_BP")  
                }
                else{
                  show("dotplot_BP")
                  hide("barplot_BP")  
                  
                }})
            }else if (input$ontology[i] == "CC"){
              ggoDEGs_CC <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_CC <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneList), keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
              output$cnetCC = renderPlot({
                cnetplot(egoDEGs_CC, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs CC")
                
              })
              output$barplot_CC <-renderPlot({
                barplot(ggoDEGs_CC, showCategory = 30,  title = "GroupGO DEGs CC")
              })
              output$dotplot_CC <-renderPlot({
                dotplot(egoDEGs_CC, title = "GO enrichment DEGs CC")
              })
              observeEvent(input$switchCC,{  #Show the barplot if switch button is on Barplot mode
                if (input$switchCC == TRUE){
                  show("barplot_CC")
                  hide("dotplot_CC")  
                }
                else{
                  show("dotplot_CC")
                  hide("barplot_CC")  
                  
                }})
            }
            
          }
        }
        
        
       
        # Convert ENSEMBL IDs to ENTREZ <===== WHY NOT USE BITR ?
        # genesENTREZ <- as.character(mapIds(input$organism, gene_list, to ="ENTREZID", from= "ENSEMBL"))
        
        ### KEGG pathway enrichment ----
        if(input$tabset == "tab2"){
          #all organisms, here the organisms might be the one selected by the user during the "choose your enrichment" step
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
        ### Msigdbr ----
        if(input$tabset == "tab3"){
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
            }
        
        ### Reactome ----
        if(input$tabset == "tab4"){
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

# SHINY APP CALL ----
shinyApp(ui = ui, server = server)