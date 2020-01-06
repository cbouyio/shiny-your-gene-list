
# IMPORT PACKAGES ---------
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
library(msigdbr) #msigdbr analysis
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
              menuItem("Analysis", tabName = "analysis", icon = icon("clipboard-check")),
              menuItem("About", tabName = "about", icon = icon("info-circle"))
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
                             withSpinner(plotOutput("barplot_BP"), type = 5,color="#0dc5c1"), plotOutput("dotplot_BP"),switchInput(
                               inputId = "switchBP",
                               value = TRUE, #true is barplot
                               onLabel = "barplot",
                               offLabel = "dotplot"),plotOutput("cnetBP")),
                         
                         box(width = 12, title ="Molecular Function", collapsible = TRUE, status = "primary",
                             withSpinner(plotOutput("barplot_MF"), type = 5,color="#0dc5c1"), plotOutput("dotplot_MF"), switchInput(
                               inputId = "switchMF",
                               value = TRUE, #true is barplot
                               onLabel = "barplot",
                               offLabel = "dotplot"),plotOutput("cnetMF")),
                         
                         box(width = 12, title = "Cellular Component", collapsible = TRUE, status = "primary",
                             withSpinner(plotOutput("barplot_CC"), type = 5,color="#0dc5c1"), plotOutput("dotplot_CC"),switchInput(
                               inputId = "switchCC",
                               value = TRUE, #true is barplot
                               onLabel = "barplot",
                               offLabel = "dotplot"),plotOutput("cnetCC"))
                       )),
                       tabPanel(value= "tab2","KEGG enrichment",fluidRow(
                         box(width = 12, title ="KEGG enrichment", collapsible = TRUE, status = "primary",
                             withSpinner(plotOutput("barplot_ekegDEGs"), type = 5,color="#0dc5c1"),
                             plotOutput("barplot_ekegMDGEs")),
                       )),        
                       
                       
                       tabPanel(value= "tab3","Msigdbr", fluidRow(
                         box(width = 12, title ="MSIGDB enrichment", collapsible = TRUE, status = "primary",
                             sliderInput("sliderMinMax", label="Size of genes, by default min = 50 and max = 500 ( Plot will update automatically)", min = 50, max = 1000, value=c(50,500), step = NULL),
                             withSpinner(plotOutput("barplot_esigDEGs"), type = 5,color="#0dc5c1"),plotOutput("dotplot_esigDEGs"),
                             switchInput(
                               inputId = "switchMsigdbr",
                               value = TRUE, #true is barplot
                               onLabel = "barplot",
                               offLabel = "dotplot")
                             ),
                         box(width = 12, title ="Gene Set enrichment", collapsible = TRUE, status = "primary",
                             plotOutput("dotplot_esigsDEGs")
                         ))),
                       tabPanel(value= "tab4","Reactome",fluidRow(
                         box(width = 12, title ="", collapsible = TRUE, status = "primary",
                             withSpinner(plotOutput("barplot_ekePDEGs"),type = 5,color="#0dc5c1"),
                             plotOutput("dotplot_ekePDEGs"))))
                       
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
  # Retrieve list of genes given by hand by user  
  geneListUser = reactive({
    if(is.null(input$file)){
      req(input$manualEntry) #check if something is written before opening it
      df <- read.csv(text=input$manualEntry,header=FALSE,col.names = c("id","log"),colClasses = c("character","numeric"))
      geneListUser=df
    } else {
      inFile <- input$file
      df <- read.csv(inFile$datapath,header=FALSE,col.names = c("id","log"),colClasses = c("character","numeric"))
      geneListUser=df
    }
    return(geneListUser)
  })
  #  observe(print(geneListUser())) #print value in terminal
  
  #Show a message window to tell the user to go to the next step when gene list is given 
  observeEvent({ 
    input$file
    input$manualEntry
    1
  }, 
  {     if (!is.null(input$file) && !is.null(input$manualEntry))
    showNotification(paste("Now choose your enrichment in the sidebar"), duration = 5, type = "message") } )
  
  observe({text_reactive()}) #text_reactive is always updated
  text_reactive <- eventReactive( input$submit, {
    # Isolate the reactive values to store them and manipulate them 
    isolate({
      ###Format of the IDs ?
      req(geneListUser()) #if exists
      
      gene_list= geneListUser()$id
      
      #If there is a one column gene list (unordered) or two column (ordered with a logFoldChange)
      if (sum(is.na(geneListUser()$log))<1){
        log_fold = geneListUser()$log
        all_gl = geneListUser()
        geneListGSEA = all_gl[,2]
        names(geneListGSEA) = as.character(all_gl[,1])
        geneListGSEA = sort(geneListGSEA, decreasing = TRUE)
        print(geneListGSEA)
        sorted_genelist= all_gl[order(all_gl$log,decreasing = TRUE) , ]
        print(all_gl)
        print(sorted_genelist)
      }
      
      
      all_types = c(ENTREZ = "ENTREZID",
                    ENSEMBL = "ENSEMBL",
                    SYMBOL = "SYMBOL")
      
      #Starting to observe if tabs are clicked...
      observeEvent(input$tabset,{
        
        
        geneConv.df = bitr(gene_list, fromType = input$idType,
                           toType = c(all_types[!all_types %in% input$idType]),
                           OrgDb = input$organism)
        #/!\ input$idType might be different than ENSEMBL !!
        rownames(geneConv.df) <- geneConv.df[["ENSEMBL"]]
        
        #remove ensembl column
        geneConv.df$ENSEMBL <- NULL
        #add only symbol elements
        geneListUserSYMB <- geneConv.df[["SYMBOL"]]
        genesENTREZ= geneConv.df[["ENTREZID"]]
        
        
        # geneConvGSEA.df = bitr(sorted_genelist, fromType = input$idType,
        #                    toType = c(all_types[!all_types %in% input$idType]),
        #                    OrgDb = input$organism)
        # print(geneConvGSEA.df)
        # 
        # rownames(geneConvGSEA.df) <- geneConvGSEA.df[["ENSEMBL"]]
        # 
        # #remove ensembl column
        # geneConvGSEA.df$ENSEMBL <- NULL
        # #add only symbol elements
        # geneListUserGSEASYMB <- geneConv.df[["SYMBOL"]]
        
        
        #        ...And execute special command for each type of enrichment/tab
        ### Geneontology enrichment ----
        if(input$tabset == "tab1"){
          ### GO enrichments ------------
          for(i in 1: length(input$ontology)){
            #For each radiobutton selected during the "choose your enrichment" step, plots will be created
            if(input$ontology[i] == "MF"){
              ggoDEGs_MF <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_MF <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneListUser), keyType = "ENSEMBL", readable = TRUE)
              
              #All the different outputs (cnetplot, barplot and dotplot)
              output$cnetMF = renderPlot({
                cnetplot(egoDEGs_MF, foldChange = geneListUserSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs MF")
              })
              output$barplot_MF <-renderPlot({
                barplot(ggoDEGs_MF, showCategory = 30,  title = "GroupGO DEGs MF")})
              output$dotplot_MF <-renderPlot({
                dotplot(egoDEGs_MF, title = "GO enrichment DEGs MF")
              })
              observeEvent(input$switchMF,{  #Show the barplot if switch button is on Barplot mode
                toggle(id = "barplot_MF", condition = input$switchMF == TRUE)
                toggle(id = "dotplot_MF", condition = input$switchMF == FALSE)
                })
              
            }else if (input$ontology[i] == "BP"){
              ggoDEGs_BP <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_BP <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneListUser), keyType = "ENSEMBL", readable = TRUE)
              output$cnetBP = renderPlot({
                cnetplot(egoDEGs_BP, foldChange = geneListUserSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs BP")
              })
              output$barplot_BP <-renderPlot({
                barplot(ggoDEGs_BP, showCategory = 30,  title = "GroupGO DEGs BP")
              })
              output$dotplot_BP <-renderPlot({
                dotplot(egoDEGs_BP, title = "GO enrichment DEGs BP")
              })
              
              observeEvent(input$switchBP,{  #Show the barplot if switch button is on Barplot mode
                toggle(id = "barplot_BP", condition = input$switchBP == TRUE)
                toggle(id = "dotplot_BP", condition = input$switchBP == FALSE)
                })
            }else if (input$ontology[i] == "CC"){
              ggoDEGs_CC <- groupGO(gene = gene_list, OrgDb = input$organism, ont = input$ontology[i], level = 2, keyType = "ENSEMBL", readable = TRUE)
              egoDEGs_CC <- enrichGO(gene = gene_list, OrgDb = input$organism, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = names(geneListUser), keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
              output$cnetCC = renderPlot({
                cnetplot(egoDEGs_CC, foldChange = geneListUserSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs CC")
                
              })
              output$barplot_CC <-renderPlot({
                barplot(ggoDEGs_CC, showCategory = 30,  title = "GroupGO DEGs CC")
              })
              output$dotplot_CC <-renderPlot({
                dotplot(egoDEGs_CC, title = "GO enrichment DEGs CC")
              })
              observeEvent(input$switchCC,{  #Show the barplot if switch button is on Barplot mode
                toggle(id = "barplot_CC", condition = input$switchCC == TRUE)
                toggle(id = "dotplot_CC", condition = input$switchCC == FALSE)
                })
            }
            #FIX ME : Deal with the case were no ontology is selected
            # else if (input$ontology[i] == NULL) {
            #   break
            # }
            
          }
        }
        
        # Convert ENSEMBL IDs to ENTREZ <===== WHY NOT USE BITR ?
        # genesENTREZ <- as.character(mapIds(input$organism, gene_list, to ="ENTREZID", from= "ENSEMBL"))
        
        ### KEGG pathway enrichment (FIXME) ----
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
          #FIXME problem "no gene can be mapped" 
          ekegMDGEs <- enrichMKEGG(gene = as.character(genesENTREZ), organism = allOrganisms_KEGG[[input$organism]], pvalueCutoff = 0.05)
          # observe(print(class(genesENTREZ))) #class is character
          observe(print(genesENTREZ)) # genes are FIXF
          output$barplot_ekegMDGEs <-renderPlot({
            barplot(ekegMDGEs, title = "DEGs KEGG modules enrichment")  # Only one "ribosome"
          })
          
        }
        ### Msigdbr ----
        if(input$tabset == "tab3"){
          data(geneList, package="DOSE") #used here for the Universe
          if (sum(is.na(geneListUser()$log))>1){
            showModal(modalDialog(
              title = "About Gene Set Enrichment",
              "GSEA cannot work with an unordered list, please change the input genes list",
              easyClose = TRUE
            ))
          }
          #geneConv2 <- bitr(names(geneList), fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
          #observe(print(head(geneConv2)))
          #univPPglu2 <- geneConv2[["SYMBOL"]]
          
          ## MSIGDB enrichment
          allOrganisms_Msigdbr <- c("Homo sapiens" ,"Mus musculus","Danio rerio","Drosophila melanogaster","Caenorhabditis elegans","Saccharomyces cerevisiae")
          names(allOrganisms_Msigdbr)= c("org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db","org.Dm.eg.db","org.Ce.eg.db","org.Sc.sgd.db")
          
          m_df <- msigdbr(species = allOrganisms_Msigdbr[[input$organism]])
          m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!
          
          # Gene enrichment
          #By default universe for enricher will be all the genes for the selected species
          observeEvent(input$sliderMinMax,{
            esigDEGs <- enricher(geneListUserSYMB, TERM2GENE = m_df, minGSSize = input$sliderMinMax[1], maxGSSize = input$sliderMinMax[2])
            output$barplot_esigDEGs = renderPlot({
              barplot(esigDEGs, showCategory = 50, title = "BarPlot DEGs MsiGDB enrichment")
            } )
            output$dotplot_esigDEGs = renderPlot ({
              dotplot(esigDEGs, showCategory = 50, title = "DotPlot DEGs MsiGDB enrichment")
            })
          })# One can play around with the set sizes to obtain something meaningfull.
          
          observeEvent(input$switchMsigdbr,{  #Show the barplot if switch button is on Barplot mode
            toggle(id = "barplot_esigDEGs", condition = input$switchMsigdbr == TRUE)
            toggle(id = "dotplot_esigDEGs", condition = input$switchMsigdbr == FALSE)
            })
          
          # Gene set enrichment
          # var = names(geneListGSEA)
          # geneConvGSEA <- bitr(var, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
          # 
          # print(geneConvGSEA)
          #          geneList should be a decreasing sorted vector... 
          # geneGSEASYMB = as.character(geneConvGSEA[["SYMBOL"]])
          # esigsDEGs <- GSEA( geneGSEASYMB, minGSSize = 20, TERM2GENE = m_df)
          # output$dotplot_esigsDEGs = renderPlot ({
          #   dotplot(esigsDEGs, showCategory = 50, title = "DotPlot MsiGDB DEGS GSEA")
          # })
        }
        # 
        ### Reactome ----
        if(input$tabset == "tab4"){
          # Enrich REACTOME Pathways
          ###/!\ A.thaliana and E.coli are not in the list, here human and yeast were written twice but will give a result error
          allOrganisms_enrichKEGG = c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly","yeast","human")
          names(allOrganisms_enrichKEGG)= c("org.Hs.eg.db","org.Mm.eg.db","org.Xl.eg.db","org.Dr.eg.db","org.Dm.eg.db","org.Ce.eg.db","org.Sc.sgd.db", "org.EcK12.eg.db","org.At.tair.db")
          
          ekePDEGs <- enrichPathway(gene = genesENTREZ, organism = allOrganisms_enrichKEGG[[input$organism]], pvalueCutoff = 0.05)
          output$barplot_ekePDEGs <-renderPlot({
            barplot(ekePDEGs, showCategory = 30, title = "BarPlot DEGs REACTOME Pathways enrichment")
          })
          output$dotplot_ekePDEGs <-renderPlot({
            dotplot(ekePDEGs, showCategory = 20, title = "DotPlot DEGs REACTOME Pathways enrichment")
          })
        }
      })
    })
  })
}

# SHINY APP CALL ----
shinyApp(ui = ui, server = server)