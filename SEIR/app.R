#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



#
# THINGS TO HAVE THEM DO WITH SHINY MODEL
#
# * Change parameters and initial variable to see if they can make epidemic not happen (Fizzle) or be more severe (fewer R at end) by playing with contact rate, recovery rate, and virulence; S0
# * At the end ask why, if epidemics eventually die-out, do we see persistent, endemic disease.


library(shiny)
library(tidyverse)
library(scales)
library(ggrepel)
library(grid)
library(gridExtra)
library(simecol)
library(diagram)

theme_set(theme_minimal())

# SIR model
SIR <- odeModel(
    main = function (time, init, parms, ...) {
        # (1) unpack variables
        S <- init["S"]
        I <- init["I"]
        R <- init["R"]
        C <- init["C"]

        N <- S+I+R

        with(as.list(parms),  {
            # Model from Chowell et al. 2004 The basic reproductive number of Ebola and the effects of public health measures: the cases of Congo and Uganda. Journal of Theoretical Biology 229:119-126.

            dS <- -beta*S*I/N
            dI <-  beta*S*I/N              -gamma*I
            dR <-                    (1-CM)*gamma*I
            dD <-                       CM *gamma*I
            dC <-               beta*S*I/N
            list(c(dS, dI, dR, dD, dC))
        })
    },
    solver = "rk4" # the function that does the solving
)


# The SEIR model
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
            dI <-               k*E        -gamma*I
            dR <-                    (1-CM)*gamma*I
            dD <-                       CM *gamma*I
            dC <-               k*E
            list(c(dS, dE, dI, dR, dD, dC))
        })
    },
    solver = "rk4" # the function that does the solving
)


# function to plot epidemics
plotEpi <- function(df){
    max_time <- max(df$time)

    ggplot(df, aes(x=time, y=Number, color = Box)) + # construct the plot
        geom_line() +
        scale_y_continuous("Number in each class",
                           labels = label_comma(accuracy=1)) +
        geom_text_repel(data=filter(df, time == max_time), aes(label = paste(Box, "=", round(Number))),
                        hjust = 0, nudge_x = 0.5,
                        direction = "y") +
        theme(legend.position = "null") +
        coord_cartesian(xlim = c(0, max_time*1.1))
}


# Define UI f
ui <- fluidPage(
    withMathJax(),
    # Application title
    titlePanel("Compartment epidemic models"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(

        sidebarPanel(

            sliderInput("beta",
                        "Transmission parameter \\(\\beta\\)",
                        min = 0.1,
                        max = 0.5,
                        value = 0.25),
            sliderInput("incubation",
                        "Incubation period in days (1/\\(k\\), in SEIR):",
                        min = 0.5,
                        max = 10,
                        value = 0.5),
            sliderInput("infectPeriod",
                        "Infectious period in days (1/\\(\\gamma\\)):",
                        min = 0.5,
                        max = 10,
                        value = 5),
            sliderInput("caseMort",
                        "Case fatality rate (proportion):",
                        min=0, max=1,
                        value = 0.1),

            shinyWidgets::sliderTextInput("N0",
                                          "Population size (\\(N_0\\)):",
                                          choices=c(10^2, 10^3, 10^4, 10^5, 10^6),
                                          selected=1000, grid = T),
            sliderInput("I0",
                        "Initial number infected (\\(I_0\\)):",
                        min = 1,
                        max = 10,
                        value = 1),
            sliderInput("fracV",
                        "Fraction vaccinated (moved from \\(S\\)->\\(R\\)):",
                        min = 0,
                        max = 1,
                        value = 0),
            sliderInput("timespan",
                        "Time over which to simulate:",
                        min = 100,
                        max = 500,
                        value = 20),
            checkboxInput("includeS",
                          "Plot the susceptible class?",
                          value = FALSE),

            withMathJax( p("\\(R_0\\) = ")),
            textOutput("R0")


        ),

        # Show plots
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("SIR model",
                                 p("The classic \"compartment\" model with frequency-dependent transmission."),
                                 helpText("Try various \"interventions\" to tame the epidemic or prevent it (i.e., keep the rate of transmission from increasing)? What fraction must be vaccinated to prevent the epidemic from taking off?"),

                                 plotOutput("sirPlots",
                                            height = "600px"),
                                 plotOutput("sirDiagram", height = "250px")
                        ),

                        tabPanel("SEIR model",
                                 p("This model includes an Exposed but not yet infectious class."),
                                 helpText("Try various \"interventions\" to tame the epidemic or prevent it (i.e., keep the rate of transmission from increasing)? What fraction must be vaccinated to prevent the epidemic from taking off?"),

                                 plotOutput("seirPlots",
                                            height = "600px"),
                                 plotOutput("seirDiagram", height = "250px")
                        )


            )
        )
    )

)


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$R0 <- renderText({
        input$beta*input$infectPeriod
    })

    output$sirDiagram <- renderPlot({
        openplotmat(asp=1/2)
        # pos <- coordinates(c(7))
        straightarrow(c(0.17, 0.65), c(0.4, 0.65), arr.pos = 1)
        # curvedarrow(c(0.45, 0.5), c(0.16, 0.5), arr.pos = 0.9, curve=0.8, endhead = TRUE)
        straightarrow(c(0.6, 0.65), c(0.75, 0.65), arr.pos = 1)
        straightarrow(c(0.5, 0.5), c(0.5, 0.2), arr.pos = 1)


        textrect(c(0.07, 0.65), lab ="S", radx = 0.05, rady = 0.1, cex = 2)
        textrect(c(0.5, 0.65), lab ="I", radx = 0.05, rady = 0.1, cex = 2)
        textrect(mid=c(0.91, 0.65), lab ="R",  radx = 0.05, rady = 0.1, cex = 2)
        textplain(mid=c(0.5, 0.06), lab ="Dead",  cex = 1.2)

        textplain(mid=c(0.3, 0.85),
                  lab = expression(beta*frac(I,N)*S), cex=1)
        textplain(mid=c(0.7, 0.85), lab = expression((1-phi)*gamma*I), cex=1)
        textplain(mid=c(0.6, 0.3), lab = expression(phi*gamma*I), cex=1)
    })

    output$seirDiagram <- renderPlot({
        openplotmat(asp=1/2)
        # pos <- coordinates(c(7))
        straightarrow(c(0.1, 0.65), c(0.28, 0.65), arr.pos = 1)
        straightarrow(c(0.4, 0.65), c(0.58, 0.65), arr.pos = 1)
        # curvedarrow(c(0.45, 0.5), c(0.16, 0.5), arr.pos = 0.9, curve=0.8, endhead = TRUE)
        straightarrow(c(0.75, 0.65), c(0.87, 0.65), arr.pos = 1)
        straightarrow(c(0.66, 0.5), c(0.66, 0.2), arr.pos = 1)


        textrect(c(0.05, 0.65), lab ="S", radx = 0.05, rady = 0.1, cex = 2)
        textrect(c(0.36, 0.65), lab ="E", radx = 0.05, rady = 0.1, cex = 2)
        textrect(c(0.66, 0.65), lab ="I", radx = 0.05, rady = 0.1, cex = 2)
        textrect(mid=c(0.95, 0.65), lab ="R",  radx = 0.05, rady = 0.1, cex = 2)
        textplain(mid=c(0.66, 0.06), lab ="Dead",  cex = 1.2)

        textplain(mid=c(0.2, 0.85),
                  lab = expression(beta*frac(I,N)*S), cex=1)
        textplain(mid=c(0.5, 0.85),
                  lab = expression(k*E), cex=1)
        textplain(mid=c(0.8, 0.85),
                  lab = expression((1-phi)*gamma*I), cex=1)
        textplain(mid=c(0.7, 0.35),
                  lab = expression(phi*gamma*I), cex=1)
    })


    output$sirPlots <- renderPlot({

        # get/set inits, pars, etc.
        init(SIR)  <- c(S=(input$N0-input$I0)*(1-input$fracV),
                        I=input$I0,
                        R=(input$N0-input$I0)*input$fracV,
                        D=0,
                        C=input$I0)
        times(SIR) <- c(from=0, to=input$timespan, by=1/4)
        parms(SIR) <- c(beta = input$beta,
                         gamma = 1/input$infectPeriod,
                         CM = input$caseMort)

        df <- out(sim(SIR)) %>% # get the output from the model
            filter(time >=0) # in case start before zero


        # Graph of epidemic, boxes by time
        df_epi <- df %>%
            gather(key=Box, value=Number, -time) %>%
            mutate(Box = factor(Box,
                                labels = c("Cases", "Dead", "Infected", "Recovered", "Susceptible")))

        if(!input$includeS) {df_epi <- filter(df_epi, Box != "Susceptible")}

        # get plot of epidemic
        P_epi <- plotEpi(df = df_epi)

        # graph of rates during epidemic through time
        df_rates <- df %>%
            mutate(Transmission = input$beta*S*I/(S+I+R),
                   Recovery = (1-input$caseMort)*I/input$infectPeriod,
                   Death = input$caseMort*I/input$infectPeriod
            ) %>%
            select(time, Transmission, Recovery, Death) %>%
            gather(key = "Process", value = "Rate", -time)

        df_labr <- filter(df_rates, time == max(time)) # for rate labels

        P_rate <- ggplot(df_rates, aes(x=time, y=Rate, color = Process)) +
            geom_line() +
            geom_text_repel(data=df_labr, aes(label = Process),
                            hjust = 0, nudge_x = 0.5,
                            direction = "y") +
            theme(legend.position = "null") +
            coord_cartesian(xlim = c(0, max(df_rates$time)*1.1))

        grid.arrange(P_epi, P_rate, ncol=1,
                     layout_matrix = rbind(c(1,1), c(1,1), c(1,1), c(2,2)))

    })

    output$seirPlots <- renderPlot({

        # get/set inits, pars, etc.
        init(SEIR)  <- c(S=(input$N0-input$I0)*(1-input$fracV),
                         E=0,
                         I=input$I0,
                         R=(input$N0-input$I0)*input$fracV,
                         D=0,
                         C=input$I0)
        times(SEIR) <- c(from=0, to=input$timespan, by=1/4)
        parms(SEIR) <- c(beta = input$beta,
                         k=1/input$incubation,
                         gamma = 1/input$infectPeriod,
                         CM = input$caseMort)

        df <- out(sim(SEIR)) %>% # get the output from the model
            filter(time >=0) # in case start before zero


        # Graph of epidemic, boxes by time
        df_epi <- df %>%
            gather(key=Box, value=Number, -time) %>%
            mutate(Box = factor(Box,
                                labels = c("Cases", "Dead", "Exposed", "Infected", "Recovered", "Susceptible")))

        if(!input$includeS) {df_epi <- filter(df_epi, Box != "Susceptible")}

        # get plot of epidemic
        P_epi <- plotEpi(df = df_epi)

        # graph of rates during epidemic through time
        df_rates <- df %>%
            mutate(Transmission = input$beta*S*I/(S+E+I+R),
                   Recovery = (1-input$caseMort)*I/input$infectPeriod,
                   Death = input$caseMort*I/input$infectPeriod
            ) %>%
            select(time, Transmission, Recovery, Death) %>%
            gather(key = "Process", value = "Rate", -time)

        df_labr <- filter(df_rates, time == max(time)) # for rate labels

        P_rate <- ggplot(df_rates, aes(x=time, y=Rate, color = Process)) +
            geom_line() +
            geom_text_repel(data=df_labr, aes(label = Process),
                            hjust = 0, nudge_x = 0.5,
                            direction = "y") +
            theme(legend.position = "null") +
            coord_cartesian(xlim = c(0, max(df_rates$time)*1.1))

        grid.arrange(P_epi, P_rate, ncol=1,
                     layout_matrix = rbind(c(1,1), c(1,1), c(1,1), c(2,2)))

    })
}

# Run the application
shinyApp(ui = ui, server = server)
