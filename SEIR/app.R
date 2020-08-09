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
library(simecol)


theme_set(theme_minimal())

# SEIR model

# The model
SEIR <- odeModel(
    main = function (time, init, parms, ...) {
        # (1) unpack variables
        S <- init["S"]
        E <- init["E"]
        I <- init["I"]
        R <- init["R"]
        C <- init["C"]

        N <- S+E+I+R

        with(as.list(parms),  {
            # Model from Chowell et al. 2004 The basic reproductive number of Ebola and the effects of public health measures: the cases of Congo and Uganda. Journal of Theoretical Biology 229:119-126.

            dS <- -beta*S*I/N
            dE <-  beta*S*I/N  -k*E
            dI <-               k*E -gamma*I
            dR <-                    gamma*I
            dC <-               k*E
            list(c(dS, dE, dI, dR, dC))
        })
    },
    solver = "rk4" # the function that does the solving
)




# Define UI f
ui <- fluidPage(

    # Application title
    titlePanel("Compartment epidemic models"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("beta",
                        "Value of transmission parameter, \\(\\beta\\):",
                        min = 0.1,
                        max = 0.5,
                        value = 0.25),
            sliderInput("incubation",
                        "Incubation time (1/k):",
                        min = 0,
                        max = 10,
                        value = 5),

            sliderInput("infectPeriod",
                        "Infectious period (1/gamma):",
                        min = 0,
                        max = 10,
                        value = 5),
            sliderInput("N0",
                        "Population size:",
                        min = 100,
                        max = 100000,
                        value = 1000),
            sliderInput("I0",
                        "Initial number infected:",
                        min = 1,
                        max = 10,
                        value = 1),
            sliderInput("timespan",
                        "Time over which to simulate:",
                        min = 10,
                        max = 100,
                        value = 50),


        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("seirPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$seirPlot <- renderPlot({
        # choose inits, pars, etc.
        inits <- c(S=input$N0-input$I0, E=0, I=input$I0, R=0, C=1)
        timespan <- c(from=0, to=input$timespan, by=1)
        pars <- c(beta = input$beta, k=1/input$incubation, gamma = 1/input$infectPeriod)


        # set parameters
        parms(SEIR) <- pars
        init(SEIR) <- inits
        times(SEIR) <- timespan

        SEIR <- sim(SEIR)

        out(SEIR) %>% # get the output from the model
            filter(time >=0) %>%
            select(-S) %>%
            gather(key=Box, value=Number, -time) %>%
            ggplot(., aes(x=time, y=Number, color = Box)) + # construct the plot
            geom_line()

    })
}

# Run the application
shinyApp(ui = ui, server = server)
