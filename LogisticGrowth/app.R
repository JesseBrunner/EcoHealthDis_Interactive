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
library(grid)
library(gridExtra)

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
makeDF <- function(r, k, timeSpan = 50, N0=1){
    df <- expand.grid(Time = seq(0,timeSpan, length.out=201), N0=N0)
    df <- df %>% mutate(
        Exponential = N0*exp(r*Time),
        Logistic = k/(1+(k/N0 -1)*exp(-r*Time))
    )

    df <- df %>%
        gather(key="Model", value="N", Exponential, Logistic)

    return(df)
}

# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Two models of population growth"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("r",
                        "Value of r (intrinsic growth rate):",
                        min = 0.05,
                        max = 0.5,
                        value = 0.25,
                        step = 0.005),
            sliderInput("k",
                        "Value of k (\"carrying capacity\"):",
                        min = 100,
                        max = 5*10^3,
                        value = 10^3,
                        step = 10),
            checkboxGroupInput("N0s",
                               "Initial population sizes (N0)",
                                c("1"="1", "100" = "100", "500"="500",
                                  "1,000"="1000", "5,000"="5000"),
                               selected = "1"),
            sliderInput("timeSpan",
                        "Time over which to run model:",
                        min = 20,
                        max = 200,
                        value = 50),
            radioButtons("yaxis",
                         "Y axis",
                         c("Linear"="liny",
                           "Logarithmic"="logy")),
            withMathJax(p("This graph shows the dynamics of exponential and logistic population growth over time.")),
            p("The exponential model is: $$\\frac{dN}{dt} = rN,$$ with a solution of:
              $$N(t) = N(0)e^{rt},$$"),
            p("The logistic model is: $$\\frac{dN}{dt} = rN \\left( 1-\\frac{N}{k} \\right) ,$$ with a solution of:$$ N(t) = \\frac{k}{1+\\left( \\frac{k}{N_0} - 1 \\right) e^{-rt}},$$"),
            withMathJax(p("where \\(r\\) is the intrinsic growth rate, \\(k\\) is the \"carrying capacity\", and \\(N(0)\\) is the initial population size."))
        ),

        mainPanel(
            tabsetPanel(type = "tabs",


                        tabPanel("Logistic vs. exponential",
                                 p("There are three plots. The first shows the population size through time while the second and third show population and per capita growth rates as a function of population size. Be sure to read the axes!"),
                                 p("See how \\(r\\) and \\(k\\) affect each of these graphs. Do you see why?"),
                                 plotOutput("modPlots", height="800px",width="500px")
                        ),


                        tabPanel("Ebola outbreaks",
                                 p("These are data on the cumulative onset of Ebola disease in two outbreaks. In both cases there was an intervention at some point in time."),
                                 p("See if you can get the logistic model to fit this time series, visually. It may not be perfect, but can you capture the general dynamics?"),
                                 radioButtons("dataset",
                                              "Data set",
                                              c("Congo 1995"="Congo",
                                                "Uganda 2000"="Uganda")),
                                 plotOutput("ebolaPlot")
                        )

            )
        )
    )
)

# Define server logic
server <- function(input, output) {

    # Tab 1: Exponential & logistic growth
    output$modPlots <- renderPlot({
        df <- makeDF(r=input$r, k=input$k, timeSpan=input$timeSpan, N0=as.numeric(input$N0s))

        # Calcuate growth rates
        df <- df %>%
            group_by(Model) %>%
            mutate(Growth = (lead(N)-N)/(lead(Time)-Time),
                   PerCapita = Growth/((N+lead(N))/2)
            )

        # Population size by time
        P <- ggplot(df, aes(Time, N, color = Model, group=interaction(Model, N0))) +
            geom_line() +
            coord_cartesian(ylim=c(1, 7500) )

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous("N(t)", labels = label_comma(accuracy=1))},
               "logy" = {P <- P+
                   scale_y_log10("N(t)", labels = label_comma(accuracy=1),
                                 breaks = c(1,5,10,50,100,500,1000,5000))}
        )
        # plot(P)

        # Population growth by N

        # find points to plot arrows along curves... 21 points through time
        timePoints <-  df$Time[10*0:20 + 1] # depends on length.out = 201 in data frame
        df_points <- filter(df, Time %in% timePoints)



        df2 <- tibble(N = seq(0, 7500, length.out = 201),
                      Exponential = input$r*N,
                      Logistic = input$r*N*(1-N/input$k)) %>%
            gather(key="Model", value="Growth", Exponential, Logistic)

        Q <- ggplot(df2, aes(N, Growth, color = Model)) +
            geom_hline(yintercept=0, color="darkgrey") +
            geom_line() +
            geom_segment(data = df_points, aes(x=N, xend=lead(N),
                                               y=Growth, yend=lead(Growth)),
                         arrow=arrow(length = unit(2, "mm")),
                         size=1/2, alpha = 1/2) +
            # geom_text(data = filter(df, Time %in% timePoints),
            #           aes(label=round(Time)),
            #           nudge_x = 1, nudge_y = 1) +
            scale_x_continuous("N(t)", labels = comma,
                               lim=c(0,5000)) +
            scale_y_continuous("Population growth rate (dN/dt)",
                               labels = comma,
                               lim=c(-25,750)) +
            labs(caption = "Arrows show trajectory over time in 20 steps. (Not all may show up)")

        # plot(Q)
        df3 <- tibble(N = seq(0, 7500, length.out = 201),
                      Exponential = input$r,
                      Logistic = input$r*(1-N/input$k)) %>%
            gather(key="Model", value="PerCapita", Exponential, Logistic)

        R <- ggplot(df3, aes(N, PerCapita, color = Model)) +
            geom_hline(yintercept=0, color="darkgrey") +
            geom_line() +
            geom_segment(data = df_points, aes(x=N, xend=lead(N),
                                               y=PerCapita, yend=lead(PerCapita)),
                         arrow=arrow(length = unit(2, "mm")),
                         size=1/2, alpha = 1/2) +
            scale_x_continuous("N(t)", labels = comma,
                               lim=c(0,5000)) +
            scale_y_continuous("Per capita growth rate ([dN/dt]/N)",
                               labels = comma,
                               lim=c(-0.1, 0.5)) +
            labs(caption = "Arrows show trajectory over time in 20 steps. (Not all may show up)")

        grid.arrange(P, Q, R, ncol=1)
    })


    # Tab 2 Ebola outbreaks
    output$ebolaPlot <- renderPlot({

        switch(input$dataset,
               "Congo" = {df_ebola <- Congo},
               "Uganda" = {df_ebola <- Uganda})

        # get the model results
        df_mod <- makeDF(r=input$r, k=input$k, timeSpan=input$timeSpan, N0=df_ebola$Cumulative[1])

        P <- ggplot(df_ebola, aes(x=Day, y=Cumulative)) +
            geom_point() +
            geom_line(data=filter(df_mod, Model=="Logistic"),
                      aes(x=Time, y=N),
                      color="blue") +
            scale_x_continuous("Days since first case",
                               lim=c(0, max(df_ebola$Day)),
                               breaks = 10*0:16)

        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P +
                   scale_y_continuous("Cumulative number of cases",
                                      labels = label_comma(accuracy=1),
                                      limits=c(0, max(df_ebola$Cumulative)*1.2)
                   )},
               "logy" = {P <- P +
                   scale_y_log10("Cumulative number of cases",
                                 labels = label_comma(accuracy=1),
                                 breaks = 10^c(0:15),
                                 limits=c(1, max(df_ebola$Cumulative)*3)
                   )}
        )

        plot(P)


    })


}

# Run the application
shinyApp(ui = ui, server = server)
