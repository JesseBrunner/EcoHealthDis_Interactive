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
makeDF <- function(r, t_intervene=50){
    df <- tibble(Time = 0:50,
                 N = 1*exp(r*Time),
                 N_intervention =
                     1*exp(r*t_intervene) +
                     ifelse(Time < t_intervene,
                            exp(r*t_intervene)*(exp(r*(Time-t_intervene))-1),
                            exp(r*t_intervene)*(exp((r/2)*(Time-t_intervene))-1)
                     )
    )

    return(df)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Exponential growth"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("r",
                        "Value of r (intrinsic growth rate):",
                        min = 0.05,
                        max = 0.5,
                        value = 0.25),
            radioButtons("yaxis",
                         "Y axis",
                         c("Linear"="liny",
                           "Logarithmic"="logy")),
            withMathJax(p("This graph shows the dynamics of exponential population growth over time. The underlying equestion is:")),
            p("$$ N(t) = N(0)e^{rt},$$"),
            withMathJax(p("where \\(r\\) is the intrinsic growth rate and \\(N(0)\\) is the initial population size"))
        ),

        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Exponential growth",
                                 p("Try adjusting \\(r\\) so that the line hits the point."),
                                 plotOutput("expPlot"),
                                 tableOutput("expTable")
                        ),
                        tabPanel("Intervening",
                                 p("Imagine that somewhere along the time series there is an intervention that reduces \\(r\\) by half. See how changing when the intervention occurs changes the"),
                                 sliderInput("t_intervene",
                                             "An intervention cutting r in half starts at:",
                                             min = 5,
                                             max = 50,
                                             value = 50),
                                 plotOutput("intPlot"),
                                 tableOutput("intTable")

                        )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {



    output$expPlot <- renderPlot({
        df <- makeDF(input$r, input$t_intervene)

        P <- ggplot(df, aes(Time, N)) +
            geom_line() +
            geom_point(data=tibble(Time=40, N=10^5))+
            theme_minimal() +
            coord_cartesian(ylim=c(1, max(10^5, max(df$N)) ))

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))},
               "logy" = {P <- P+
                                 scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                               breaks = 10^c(0:15))}
        )

        plot(P)

    })

    output$intPlot <- renderPlot({
        df <- makeDF(input$r, input$t_intervene)

        P <- ggplot(df, aes(Time, N)) +
            geom_line() +
            geom_line(aes(y=N_intervention), color = "blue") +
            geom_point(data=tibble(Time=40, N=10^5))+
            theme_minimal() +
            coord_cartesian(ylim=c(1, max(10^5, max(df$N)) ))

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))},
               "logy" = {P <- P+
                   scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                 breaks = 10^c(0:15))}
        )

        plot(P)

    })

    output$expTable <- renderTable({
        df <- makeDF(input$r, input$t_intervene)
        # format table to use commas in numbers
        df$N <- comma(df$N, accuracy = 0.1)
        arrange( df, desc(Time))
    },
    digits=1, align = "r")

    output$intTable <- renderTable({
        df <- makeDF(input$r, input$t_intervene)
        # format table to use commas in numbers
        df$N <- comma(df$N, accuracy = 0.1)
        df$N_intervention <- comma(df$N_intervention, accuracy = 0.1)
        arrange( df, desc(Time))
    },
    digits=1, align = "r")
}

# Run the application
shinyApp(ui = ui, server = server)
