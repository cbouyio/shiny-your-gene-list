#server.R
library(shiny)
function(input, output, session) {
  #save output$X  and plotOutput("X")
  output$hist = renderPlot({
    title = "title"
    hist(rnorm(input$num), main =title)})
}

