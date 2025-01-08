#
library(shiny)
library(leaflet)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata)

### pre stuff
obs_data<- read_csv("worms_ll.csv")
obs_data$Selected <- "No"
bioclim_data_sm <- rast("bioclim_data_sm.tif")

# Determine geographic extent of our data
max_lat <- ceiling(max(obs_data$latitude))
min_lat <- floor(min(obs_data$latitude))
max_lon <- ceiling(max(obs_data$longitude))
min_lon <- floor(min(obs_data$longitude))


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Sketch out the realized niche of hookworms"),
  h4("Locations with ≥ 5% hookworm prevalence"),
  p("Data from",  a("Fleitas et al. (2022) Parasites & Vectors 15: 197", 
                    href="https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-022-05284-w")),
  leafletOutput("map"),
  
  p("Click on locations in the map to add their conditions to the niche space graph"),
  h4("Niche space"),
  
  plotOutput("plot")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  selected <- 0
  
  sel <- reactive({
    selected <<- c(selected, input$map_shape_click[["id"]])
    # print(selected)
    return(selected)
  })
  
  pal <- colorNumeric("RdBu", reverse = TRUE, domain = c(-20,31), na.color = NA) 
  
  output$map <- renderLeaflet({
    leaflet(obs_data %>% mutate(Selected = "No")) %>% addTiles() %>%
      setMaxBounds(lng1=-180, lat1=-60, lng2=180, lat2=60) %>% 
      fitBounds(lng1=min_lon, lat1=min_lat, lng2=max_lon, lat2=max_lat) %>% 
      addRasterImage(bioclim_data_sm[["wc2.1_2.5m_bio_1"]], opacity = 0.8, colors = pal) %>%
      addLegend(title = "Ave.\ntemp.\n(°C)", pal = pal, values = 5*(0:6)) %>% 
      addCircles(layerId = obs_data$ID,
                 weight = 10,
                 popup = ~htmltools::htmlEscape(paste0(
                   #round(latitude, 2), "°, ", round(longitude, 2), "° \n", 
                   round(Temperature,1), "°C; ", 
                   round(Precipitation), "mm")
                 )
      ) 
  })
  
  observe({
    selected <- sel()
    # print(selected)
    df <- obs_data %>%
      mutate(Selected = ifelse(ID %in% selected, "Yes", "No"),
             Selected = ifelse(ID == last(selected), "Last", Selected)
      )
    
    leafletProxy("map", data = df) %>%
      clearShapes() %>%
      addCircles(layerId = obs_data$ID,
                 color = ifelse(df$Selected == "Last",
                                "red",
                                ifelse(df$Selected == "Yes",
                                       "black",
                                       "blue")),
                 fill = FALSE,
                 fillColor = "white",
                 radius = 5,
                 weight = ifelse(df$Selected == "Yes", 5, 10),
                 popup = ~htmltools::htmlEscape(paste0(round(Temperature,1), "°C; ",
                                                       round(Precipitation), "mm")
                 ),
                 # popupOptions = c(closeOnClick = TRUE) # <- causes problems with Input to asJSON(keep_vec_names=TRUE) is a named vector. In a future version of jsonlite, this option will not be supported, and named vectors will be translated into arrays instead of objects. If you want JSON object output, please use a named list instead. See ?toJSON.
      )
  })
  
  output$plot <- renderPlot({
    selected <- sel()
    df <- obs_data %>% 
      filter(ID %in% selected) %>%
      mutate(sel = if_else(ID == last(selected), "No", "Yes")
      )
    ggplot(df, aes(Temperature, Precipitation)) + 
      geom_polygon(data=df[chull(df$Temperature, df$Precipitation),], 
                   color = "darkgray", fill = NA,
      ) + 
      geom_point(aes(size = sel, color = sel), alpha = 0.5) +
      scale_x_continuous("Average annual temperature (°C)", 
                         limits = c(min(obs_data$Temperature)-5, 
                                    max(obs_data$Temperature)+5)) + 
      scale_y_continuous("Average annual precipitation (mm)", 
                         limits = c(min(obs_data$Precipitation)-500, 
                                    max(obs_data$Precipitation)+500)) + 
      scale_color_manual(values = c("red", "black"), guide = "none") + 
      scale_size_manual(values = c(5,3), guide = FALSE) + 
      theme_bw() +
      # labs(subtitle = "Niche space") + 
      theme(text = element_text(size = 20))
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
