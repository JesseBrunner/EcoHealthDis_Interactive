library(terra)
library(geodata)

# Data on hookworm distribution from 
# Mapping the global distribution of Strongyloides stercoralis and hookworms by ecological niche modeling
# https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-022-05284-w
# Additional file 2: Table S2: Surveys used for the ecological niche model of hookworms.

obs_data <- read_csv("worms.csv")

obs_data <- obs_data %>% filter(Prevalence > 5) %>% 
  select(longitude, latitude) %>% 
  distinct()

# summary(obs_data)

# following this SDM guide: https://jcoliver.github.io/learn-r/011-species-distribution-models.html

bioclim_data <- worldclim_global(var = "bio",
                                 res = 2.5,
                                 path = "data/")
# see https://www.worldclim.org/data/bioclim.html for codes
# 1 is annual mean temperature
# 12 is annual precipitation
# and https://www.worldclim.org/data/worldclim21.html for all historical data


bioclim_extract <- extract(x = bioclim_data,
                           y = obs_data[, c("longitude", "latitude")],
                           ID = TRUE) # No need for an ID column

df <- bind_cols(
  obs_data, 
  bioclim_extract %>% 
    select(ID, wc2.1_2.5m_bio_1, wc2.1_2.5m_bio_12) %>% 
    rename(Temperature = wc2.1_2.5m_bio_1, Precipitation = wc2.1_2.5m_bio_12)
)
write_csv(df, "worms_ll.csv")


s <- rast(nrows = 1000, ncols = 2500, xmin=-180, xmax=180, ymin=-70, ymax=70)

bioclim_data_sm <- terra::resample(bioclim_data[[c(1,12)]], s)
writeRaster(bioclim_data_sm, "bioclim_data_sm.tif", overwrite = TRUE)



             