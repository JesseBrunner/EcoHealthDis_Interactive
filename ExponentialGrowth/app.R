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
library(scales)

# Pre-Shiny stuff
makeDF <- function(r){
    tibble(Time = 0:50,
           N = 1*exp(r*Time))
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Exponential growth"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("r",
                        "value of r (intrinsic growth rate):",
                        min = 0.05,
                        max = 0.5,
                        value = 0.15),
            radioButtons("yaxis",
                         "Y axis",
                        c("Linear"="liny",
                          "Logarithmic"="logy")),
            withMathJax(p("This graph shows the dynamics of exponential population growth over time. The underlying equestion is:")),
            p("$$ N(t) = N(0)e^{rt},$$"),
            withMathJax(p("where \\(r\\) is the intrinsic growth rate and \\(N(0)\\) is the initial population size")),
            p("Try adjusting \\(r\\) so that the line hits the point")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("expPlot"),
           tableOutput("expTable")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {



    output$expPlot <- renderPlot({
        df <- makeDF(input$r)

        P <- ggplot(df, aes(Time, N)) +
            geom_line() +
            geom_point(data=tibble(Time=40, N=10^5))+
            theme_minimal() +
            coord_cartesian(ylim=c(1, max(10^5, max(df$N)) ))

        switch(input$yaxis,
               "liny" = plot(P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))),
               "logy" = plot(P+
                                 scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                               breaks = 10^c(0:15)))
               )

    })

    output$expTable <- renderTable({
        arrange( makeDF(input$r), desc(Time))
        },
        digits=1)
}

# Run the application
shinyApp(ui = ui, server = server)
