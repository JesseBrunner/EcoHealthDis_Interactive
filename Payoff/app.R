#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)


Repay <- function(InitialBalance, APR, MonthlyPayment){
  MPR <- APR/12
  # Starting values
  Time <- 0
  Balance <- InitialBalance
  Interest <- 0

  while(last(Balance) > 0 & last(Balance) < 10^8){
    Time[Time+2] <- Time + 1
    Interest[last(Time)+1] <- last(Balance)*MPR
    Balance[last(Time)+1] <- last(Balance) + last(Interest) - MonthlyPayment
  }
  df <- tibble(Time, Interest, Balance, "Cumulative Interest" = cumsum(Interest))
  df$Balance[df$Balance < 0] <- 0
  return(df)

}


# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Paying off debts, one month at a time"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      numericInput("Init",
                   "Initial amount:",
                   value = 1000),
      numericInput("APR",
                   "Annual percentage rate (as a %):",
                   value = 10),

      numericInput("Payment",
                   "Monthly payment amount:",
                   value = 100),

      br(),
      actionButton("button", "Run"),

      br(),
      p("Note: this is a simple calculator where each month the payment is applied at the same time the interest (APR/12) is applied. There are no daily running balances or anything tricky some loans or credit cards involve.")
    ),

    # Show a plot of the generated distribution
    mainPanel(

      plotOutput("outPlot"),
      dataTableOutput("outTable")

    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  observeEvent(input$button, print(""))

  df <- eventReactive(input$button, {
    Repay(InitialBalance = input$Init,
              APR = input$APR/100,
              MonthlyPayment = input$Payment)
  })


    output$outPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
      df() %>%
        pivot_longer(cols = c(Balance, `Cumulative Interest`),
                     values_to = "Amount") %>%
        ggplot(., aes(x=Time, y = Amount, color = name)) +
        geom_step() +
        scale_x_continuous("Months") +
        scale_color_brewer("", palette = "Dark2") +
        theme(legend.position = "bottom")
    })

    output$outTable <- renderDataTable(round(df(),2))
}

# Run the application
shinyApp(ui = ui, server = server)
