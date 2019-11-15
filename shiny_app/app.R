library(shiny)

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
      textAreaInput("manualEntry", "...Or write your list of genes", "geneID1, geneID2, etc.", width = "200px"),
      verbatimTextOutput("value"),

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
                   selected = ",")
      ),
      # radioButtons("quote", "Quote",
      #              choices = c(None = "",
      #                          "Double Quote" = '"',
      #                          "Single Quote" = "'"),
      #              selected = '"'),
      
      #),
      # Horizontal line ----
      tags$hr(),
      selectInput("organism", "Select organism(s)", list(
        "Vertebrate"  = c("Human (Homo sapiens)"="human","Mouse (Mus musculus)"="mouse",
                          "Western clawed frog (Xenopus tropicalis)"="xenopus","Zebrafish (Danio rerio)"="danio"
                          ),
        "Eukaryotic" = c("Fruit fly (Drosophila melanogaster)"="droso",
                          "Nematode worm (Caenorhabditis elegans)"="celegans","Yeast (Saccharomyces cerevisiae)"="saccharomyces"),
        "Prokaryotic" = c("Bacteria (Escherichia coli)"="ecoli")
        
      ), selected = NULL, multiple = TRUE)),
    # Main panel for displaying outputs ----
      mainPanel(
      h3("This is the main title of the page"),
      )
      # Output: Data file ----
      #tableOutput("contents")
      
    
  )
)



####SERVER####################

# Define server logic to read selected file ----
server <- function(input, output) {
  output$value <- renderText({input$manualEntry })
  output$file_listener <- reactive({
    return(!is.null(input$file))
  })
  #Listener is always active 
  outputOptions(output, "file_listener", suspendWhenHidden = FALSE)
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(input$file)
    #Check if the user wants to upload a file or paste a list of genes
    if(is.null(input$file) == FALSE){
      df <- read.csv(input$file$datapath,
                     header = input$header,
                     sep = input$sep)
                     #quote = input$quote)
      
    }
    else{
      output$cleanValue = as.list(strsplit(output$value, ",")[[1]])
      #renderPrint(output$cleanValue)
    }
  })
}
# Run the app ----
shinyApp(ui, server)