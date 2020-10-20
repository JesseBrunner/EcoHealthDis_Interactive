#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.


library(shiny)
library(googlesheets4)
library(tidyverse)
library(scales) # for prettier axis labels
library(ggrepel) # for non-overlapping labels

theme_set(theme_minimal())

# prev <- read_csv("LifetimePreventable.csv")
gs4_deauth()
prev <- read_sheet("1xyH3-JmurBdUAqTeion9B2O7DXpVoKB4jKqHZig52Is")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Odds of death by..."),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            radioButtons("yaxis",
                         "Y axis",
                         c("Linear"="liny",
                           "Logarithmic"="logy")),
            uiOutput("choose_vars")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot", height = "600px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$choose_vars <- renderUI({


        Vars <- unique(prev$Event)
        Vars <- Vars[order(Vars)]

        # Create the radio buttons and select the first by default
        checkboxGroupInput("vars", "Choose the events to display",
                     choices  = Vars,
                     selected = c("Heart disease", "Lightning"))

    })

    output$distPlot <- renderPlot({

        P <- prev %>%
            filter(Event %in% input$vars) %>%
            ggplot(., aes(reorder(Event, -Odds), Odds, color=Per)) +
            geom_point() +
            geom_text_repel(aes(label = Event),
                            hjust = 0, nudge_x = 0.5,
                            direction = "y"
            ) +
            scale_x_discrete("") +
            scale_color_brewer(palette = "Dark2") +
            # facet_grid(.~Per, labeller = label_both) +
            theme_minimal() +
            theme(axis.text.x = element_blank())


        # get the right y-axis
        switch(input$yaxis,
               "liny" = {P <- P+scale_y_continuous(
                   breaks = 1/c(2,6,10,20,50,100, 5000),
                   labels=label_math(1:.x,
                                     format = function(x){round(1/x)}
                   )
               )
               },
               "logy" = {P <- P+scale_y_log10(
                   labels=label_math(1:.x,
                                     format = function(x){round(1/x)})
               )
               }
        )
        plot(P)
    })
}

# Run the application
shinyApp(ui = ui, server = server)
