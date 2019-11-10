#write shinyapp will automaticaly load
library(shiny)

ui <- fluidPage(
                sliderInput(inputId = "num", 
                            label = "Choose a number",
                            value = 25,
                            min = 1,
                            max = 100),
                plotOutput("hist"),
                plotOutput("hist2")
  
)

#Build the output 
server <- function(input, output, session) {
  #save output$X  and plotOutput("X")
  output$hist = renderPlot({
  title = "title"
  hist(rnorm(input$num), main =title)})
  output$hist2 = renderPlot({
  title = "title"
  hist(dnorm(input$num), main =title)})
}

shinyApp(ui, server)
