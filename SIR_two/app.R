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
library(deSolve)

theme_set(theme_minimal())


# SIR model
SIR <- odeModel(
  main = function (time, init, parms, ...) {
    # (1) unpack variables
    Sa <- init["Sa"]
    Ia <- init["Ia"]
    Sb <- init["Sb"]
    Ib <- init["Ib"]

    Ra <- init["Ra"]
    Da <- init["Da"]

    Rb <- init["Rb"]
    Db <- init["Db"]

    Ca <- init["Ca"]
    Cb <- init["Cb"]

    N <- Sa+Sb+Ia+Ib+Ra+Rb

    with(as.list(parms),  {


      dSa <- -(beta_aa*Ia + beta_ba*Ib)*Sa/N
      dIa <-  (beta_aa*Ia + beta_ba*Ib)*Sa/N - gamma_a*Ia
      dSb <- -(beta_ab*Ia + beta_bb*Ib)*Sb/N
      dIb <-  (beta_ab*Ia + beta_bb*Ib)*Sb/N - gamma_b*Ib
      dRa <-                    (1-CM_a)*gamma_a*Ia
      dDa <-                       CM_a *gamma_a*Ia
      dRb <-                    (1-CM_b)*gamma_b*Ib
      dDb <-                       CM_b *gamma_b*Ib
      dCa <-   (beta_aa*Ia + beta_ba*Ib)*Sa/N
      dCb <-   (beta_ab*Ia/N + beta_bb*Ib/N)*Sb
      list(c(dSa, dIa, dSb, dIb, dRa, dDa, dRb, dDb, dCa, dCb))
    })
  },
  solver = "rk4" # the function that does the solving
)



# function to plot epidemics
plotEpi <- function(df){
  max_time <- max(df$time)

  ggplot(df, aes(x=time, y=Number, color = Box, linetype=Population)) + # construct the plot
    geom_line() +
    scale_y_continuous("Number in each class",
                       labels = label_comma(accuracy=1)) +
    geom_text_repel(data=filter(df, time == max_time), aes(label = paste(Box, "=", round(Number))),
                    hjust = 0, nudge_x = 0.5,
                    direction = "y") +
    scale_linetype_manual(values = c(1,2)) +
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

      # mixing
      p("Mixing parameters"),

      sliderInput("beta_ab",
                  "Transmission from A to B \\(\\beta_{A\\rightarrow B}\\)",
                  min = 0,
                  max = 1,
                  value = 0),

      sliderInput("beta_ba",
                  "Transmission from B to A \\(\\beta_{B\\rightarrow A}\\)",
                  min = 0,
                  max = 1,
                  value = 0),

      # Population A
      tags$hr(style="border-color: black;"),
      p("Population A: solid lines"),
      sliderInput("beta_aa",
                  "Transmission within population A \\(\\beta_{A\\rightarrow A}\\)",
                  min = 0.1,
                  max = 1,
                  value = 0.5),
      sliderInput("infectPeriod_a",
                  "Infectious period in days (1/\\(\\gamma_A\\)):",
                  min = 0.5,
                  max = 10,
                  value = 5),
      sliderInput("caseMort_a",
                  "Case fatality rate (proportion of A that die):",
                  min=0, max=1,
                  value = 0.1),

      tags$hr(style="border-color: black;"),
      p("Population B: dashed lines"),
      sliderInput("beta_bb",
                  "Transmission within population B \\(\\beta_{B\\rightarrow B}\\)",
                  min = 0.1,
                  max = 1,
                  value = 0.5),
      sliderInput("infectPeriod_b",
                  "Infectious period in days (1/\\(\\gamma_B\\)):",
                  min = 0.5,
                  max = 10,
                  value = 5),
      sliderInput("caseMort_b",
                  "Case fatality rate (proportion of B that die):",
                  min=0, max=1,
                  value = 0.1),

      tags$hr(style="border-color: black;"),
      p("Initial conditions"),
      sliderInput("I0_a",
                  "Initial number of A infected (\\(I_{A_0}\\)):",
                  min = 1,
                  max = 10,
                  value = 1),

      shinyWidgets::sliderTextInput("N0_a",
                                    "Size of population A (\\(N_{A_0}\\)):",
                                    choices=c(10^2, 500, 10^3, 5000, 10^4),
                                    selected=1000, grid = T),

      shinyWidgets::sliderTextInput("N0_b",
                                    "Size of population B (\\(N_{B_0}\\)):",
                                    choices=c(10^2, 500, 10^3, 5000, 10^4),
                                    selected=1000, grid = T),

      sliderInput("timespan",
                  "Time over which to simulate:",
                  min = 100,
                  max = 500,
                  value = 200)

    ),

    # Show plots
    mainPanel(

      p("The classic \"compartment\" model with frequency-dependent transmission within and (potentially) between two populations."),
      helpText("The infection starts in population A. How much transmission from A to B (or vice versa) does it take before population B is strongly affected?"),
      helpText("Does it matter of the populations are really different in their transmission within the population, infectious period, or case fatality rate?"),
      helpText("What if A is really infectious to B, or vice versa?"),
      plotOutput("sirPlots",
                 height = "600px")
      ,
      plotOutput("sirDiagram", height = "400px")
    )
  )

)


# Define server logic required to draw a histogram
server <- function(input, output) {


  output$sirDiagram <- renderPlot({
    openplotmat(asp=1/2)
    # pos <- coordinates(c(7))
    straightarrow(c(0.1, 0.8), c(0.45, 0.8), arr.pos = 1)
    # curvedarrow(c(0.45, 0.5), c(0.16, 0.5), arr.pos = 0.9, curve=0.8, endhead = TRUE)
    straightarrow(c(0.6, 0.8), c(0.8, 0.8), arr.pos = 1)
    straightarrow(c(0.55, 0.7), c(0.7, 0.55), arr.pos = 1)


    textrect(c(0.05, 0.8), lab =expression(S[A]), radx = 0.05, rady = 0.1, cex = 2)
    textrect(c(0.55, 0.8), lab =expression(I[A]), radx = 0.05, rady = 0.1, cex = 2)
    textrect(mid=c(0.91, 0.8), lab =expression(R[A]),  radx = 0.05, rady = 0.1, cex = 2)
    textplain(mid=c(0.77, 0.5), lab ="Dead",  cex = 1.2)

    textplain(mid=c(0.3, 0.9),
              lab = expression(bgroup("(", beta[AA]*frac(I[A],N) + beta[BA]*frac(I[B],N), ")")*S[A]), cex=1)
    textplain(mid=c(0.72, 0.9), lab = expression((1-phi[A])*gamma[A]*I), cex=1)
    textplain(mid=c(0.7, 0.65), lab = expression(phi[A]*gamma[A]*I[A]), cex=1)

    ## Group B

    straightarrow(c(0.1, 0.2), c(0.45, 0.2), arr.pos = 1)
    # curvedarrow(c(0.45, 0.5), c(0.16, 0.5), arr.pos = 0.9, curve=0.8, endhead = TRUE)
    straightarrow(c(0.6, 0.2), c(0.8, 0.2), arr.pos = 1)
    straightarrow(c(0.55, 0.2), c(0.7, 0.44), arr.pos = 1)


    textrect(c(0.05, 0.2), lab =expression(S[B]), radx = 0.05, rady = 0.1, cex = 2)
    textrect(c(0.55, 0.2), lab =expression(I[B]), radx = 0.05, rady = 0.1, cex = 2)
    textrect(mid=c(0.91, 0.2), lab =expression(R[B]),  radx = 0.05, rady = 0.1, cex = 2)

    textplain(mid=c(0.3, 0.1),
              lab = expression(bgroup("(", beta[BB]*frac(I[B],N) + beta[AB]*frac(I[A],N), ")")*S[B]), cex=1)
    textplain(mid=c(0.72, 0.1), lab = expression((1-phi[B])*gamma[B]*I), cex=1)
    textplain(mid=c(0.7, 0.3), lab = expression(phi[B]*gamma[B]*I[B]), cex=1)

  })


  output$sirPlots <- renderPlot({

    # get/set inits, pars, etc.
    init(SIR)  <- c(Sa=(input$N0_a-input$I0_a),
                    Ia=input$I0_a,
                    Sb=input$N0_b,
                    Ib=0,

                    Ra=0,
                    Da=0,

                    Rb=0,
                    Db=0,
                    Ca=input$I0_a,
                    Cb=0
    )

    times(SIR) <- c(from=0, to=input$timespan, by=1/5)
    parms(SIR) <- c(beta_aa = input$beta_aa,
                    beta_bb = input$beta_bb,
                    beta_ab = input$beta_ab,
                    beta_ba = input$beta_ba,
                    gamma_a = 1/input$infectPeriod_a,
                    gamma_b = 1/input$infectPeriod_b,
                    CM_a = input$caseMort_a,
                    CM_b = input$caseMort_b)

    df <- out(sim(SIR)) %>% # get the output from the model
      filter(time >=0) # in case start before zero


    # Graph of epidemic, boxes by time
    df_epi <- df %>%
      gather(key=Box, value=Number, -time) %>%
      mutate(Population = ifelse(str_ends(Box, "a"), "A", "B"),
             Box = factor(str_sub(Box, 1, 1), labels = c("Cases", "Deaths", "Infecteds", "Recovereds", "Susceptibles"))
      )


    # get plot of epidemic
    P_epi <- plotEpi(df = df_epi)

    # graph of rates during epidemic through time
    df_rates <- df %>%
      mutate(Transmission = ( (input$beta_aa*Ia + input$beta_ba*Ib)*Sa +
                                (input$beta_ab*Ia + input$beta_bb*Ib)*Sb )/(Sa+Sb+Ia+Ib+Ra+Rb) ,
             Recovery = (1-input$caseMort_a)*Ia/input$infectPeriod_a +
               (1-input$caseMort_b)*Ib/input$infectPeriod_b,
             Death = input$caseMort_a*Ia/input$infectPeriod_a +
               input$caseMort_b*Ib/input$infectPeriod_b
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
