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

theme_set(theme_minimal())


# Data import
Congo <- read_csv("data/Ebola_Congo_1995.csv")
Congo <- Congo %>%
    mutate(Cumulative = cumsum(Cases)) %>%
    filter(Day >=7)
Uganda <- read_csv("data/Ebola_Uganda_2000.csv")
Uganda <- Uganda %>%
    mutate(Cumulative = cumsum(Cases))

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
                           "Logarithmic (base 10)"="logy",
                           "Logarithmic (natural log)" = "lny")),
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

                        ),


                        tabPanel("Ebola outbreaks",
                                 p("These are data on the cumulative onset of Ebola disease in two outbreaks. In both cases there was an intervention at some point in time."),
                                 p("This panel fits a linear regression to the data up to the intervention start date (value in slider) and then projects that line forward."),
                                 p("See if you can guess when that intervention started. (Note that unlike our simple example in the previous tab, interventions do not happen instantaneously everywhere, so it may not look sudden.)"),
                                 radioButtons("dataset",
                                              "Data set",
                                              c("Congo 1995"="Congo",
                                                "Uganda 2000"="Uganda")),
                                 sliderInput("cutoff",
                                             "Intervention started at day:",
                                             min = 20,
                                             max = 100,
                                             value = 50),
                                 plotOutput("ebolaPlot")
                        )

            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    # Tab 1: Exponential growth
    output$expPlot <- renderPlot({
        df <- makeDF(input$r, input$t_intervene)

        P <- ggplot(df, aes(Time, N)) +
            geom_line() +
            geom_point(data=tibble(Time=40, N=10^5))+
            coord_cartesian(ylim=c(1, max(10^5, max(df$N)) ))

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))},
               "logy" = {P <- P+
                                 scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                               breaks = 10^c(0:15))},
               "lny" = {P <- P + scale_y_continuous("N(t)",
                                                    trans = log_trans(),
                                                    labels = label_comma(accuracy=1),
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

    # Tab 2: Intervention in exponential growth
    output$intPlot <- renderPlot({
        df <- makeDF(input$r, input$t_intervene)

        P <- ggplot(df, aes(Time, N)) +
            geom_line() +
            geom_line(aes(y=N_intervention), color = "blue") +
            geom_point(data=tibble(Time=40, N=10^5))+
            coord_cartesian(ylim=c(1, max(10^5, max(df$N)) ))

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))},
               "logy" = {P <- P+
                   scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                 breaks = 10^c(0:15))},
               "lny" = {P <- P + scale_y_continuous("N(t)", trans = log_trans(),
                                                    labels = label_comma(accuracy=1),
                                                    breaks = 10^c(0:15))}
        )
        plot(P)
    })

    output$intTable <- renderTable({
        df <- makeDF(input$r, input$t_intervene)
        # format table to use commas in numbers
        df$N <- comma(df$N, accuracy = 0.1)
        df$N_intervention <- comma(df$N_intervention, accuracy = 0.1)
        arrange( df, desc(Time))
    },
    digits=1, align = "r")


    # Tab 3 Ebola outbreaks
    output$ebolaPlot <- renderPlot({

        switch(input$dataset,
               "Congo" = {df <- Congo},
               "Uganda" = {df <- Uganda})

        # get right y-values
        df <- df %>%
            mutate(y = Cumulative)
        if(input$yaxis == "logy") {
            df <- df %>%
            mutate(y = log(Cumulative))
        }

        cutoff <- input$cutoff # choose cutoff

        # fit regression & get slope and standard error
        regr <- summary(
            lm(y ~ Day, data = filter(df, Day <= cutoff))
            )$coefficients["Day", 1:2]
        regr <- round(regr, 3)



        P <- ggplot(df, aes(x=Day, y=Cumulative)) +
            geom_point() +
            geom_smooth(method="lm",
                        data=filter(df, Day <= cutoff),
                        fullrange=TRUE, se=F, linetype=2) +
            geom_smooth(method="lm",
                        data=filter(df, Day <= cutoff)) +
            scale_x_continuous("Days since first case", breaks = 10*0:16) +
            annotate("text", x=10, y=max(df$Cumulative),
                     label = bquote(slope == .(regr[1]) %+-% .(regr[2])), hjust=0 )

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P +
                   scale_y_continuous("Cumulative number of cases",
                                      labels = label_comma(accuracy=1),
                                      limits=c(0, max(df$Cumulative)*1.2)
                   )},
               "logy" = {P <- P +
                   scale_y_log10("Cumulative number of cases",
                                 labels = label_comma(accuracy=1),
                                 breaks = 10^c(0:15),
                                 limits=c(1, max(df$Cumulative)*3)
                   )},
               "lny" = {P <- P + scale_y_continuous("Cumulative number of cases",
                                                    trans = log_trans(),
                                                    labels = label_comma(accuracy=1),
                                                    breaks = 10^c(0:15),
                                                    limits=c(1, max(df$Cumulative)*3))}
        )

        plot(P)

    })


}

# Run the application
shinyApp(ui = ui, server = server)
