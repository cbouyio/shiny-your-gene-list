library(shiny)
library(shinydashboard)

body <- dashboardBody(
  fluidRow(
    tabBox(
      title = NULL, width = 12, height = 6,
      # The id lets us use input$tabset1 on the server to find the current tab
      id = "tabset1", height = "250px",
      tabPanel("Victim", "Victim tab", tabBox(title ="toto",plotOutput("s1"))),
      tabPanel("Trafficker", "Trafficker tab")
    )
  ),
  fluidRow(infoBoxOutput("tabset1Selected"))
)

shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "Human Trafficking"),
    dashboardSidebar(disable = TRUE),
    body
  ),
  server = function(input, output) {
    # The currently selected tab from the tab box
    output$tabset1Selected <- renderInfoBox({
      infoBox("Selected Tab", input$tabset1, icon = icon("info-circle"))
    })
    output$sl <- renderUI({
      sliderInput(
        "obs",
        "Number of observations:",
        min = 0, max = 1000, value = 500
      )
    })
  }
)