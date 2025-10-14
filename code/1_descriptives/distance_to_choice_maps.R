rm(list=ls())
library(geosphere)
library(tidycensus)
library(dplyr)
library(sf)
library(ggplot2)
library(haven)
library(magick)
library(pdftools)
# ---- SETUP ------------------------------------------------------------------

# 1. Get all census blocks in Los Angeles County with geometry
#    (If you haven't done so, set your API key with census_api_key)

blocks20_raw <- get_decennial(
  geography = "block",
  variables = "P1_001N",  # or any other decennial variable
  year = 2020,
  state = "CA",
  county = "Los Angeles",
  geometry = TRUE
)

# 2. Read in your distance data
geodir <- "/Volumes/lausd/magnet/data/raw/geodata/"


# 3. Read in the attendance zone shapefiles and transform them
attendance_zone_ms_key_shp <- paste0(geodir, "LAUSD_Attendance_Boundary_(Middle_Schools).shp")
attendance_zones_ms <- st_read(attendance_zone_ms_key_shp) %>% st_transform(crs = 4326)


setwd("/Volumes/lausd/decentralized_choice/output/figures/")
# This color scale can be reused in each loop:
my_scale <- scale_fill_viridis_c(
  option="rocket",
  na.value = "gray50",
  limits = c(0, 5),            # fix the color scale from 0 to 15
  breaks = c(1, 2, 5),        # breaks at 5, 10, 15
  oob = scales::squish         # values over 10 get "squished" to 10
)

my_scale <- scale_fill_viridis_c(
  option = "rocket",
  begin  = 0.2,     # skip the darkest portion (near black)
  end    = 0.9,     # capture more of the brighter reds/oranges
  direction = 1,    # or -1, if you want to invert the gradient
  na.value = "gray50",
  limits   = c(0, 5),
  breaks   = c(1, 2, 5),
  oob      = scales::squish
)

my_scale <- scale_fill_gradientn(
  colours = c("black", "#800000", "#B22222", "#CD5C5C", "#C2B280", "#D2B48C", "#F5DEB3", "#DEB887" ),
  na.value = "gray50",
  limits = c(0, 5),
  breaks = c(1, 2, 5),
  oob = scales::squish
)

my_scale <- scale_fill_gradientn(
  colours = c("#000000", "#800000", "#CD5C5C", "#F5DEB3"),
  na.value = "gray50",
  limits = c(0, 5),
  breaks = c(1, 2, 3, 5),
  oob = scales::squish
)

my_scale2 <- scale_fill_gradientn(
  colours = c("#000000", "#800000", "#CD5C5C", "#F5DEB3"),
  na.value = "gray50",
  limits = c(0, 3),
  breaks = c(0.5,1.5, 2.5),
  oob = scales::squish
)


# ---- LOOP THROUGH YEARS: MIDDLE SCHOOL --------------------------------------
years <- c(2004,2023)
distance_dir <- "/Volumes/lausd/build/rawdata/aux_data/"
distance_dat <- read_dta(paste0(distance_dir, "relative_choice_distances_grade5.dta"))

for (yr in years) {
  distance_temp <- distance_dat %>%
    filter(endyear == yr, schooltype == "ms") %>%
    select(censusblockid, min_choice_dist, schooltype,
           rel_choice_distance)
  
  block_temp <- blocks20_raw %>%
    left_join(distance_temp, by = c("GEOID" = "censusblockid")) %>%
    mutate(
      centroid = st_centroid(geometry),
      lon = st_coordinates(centroid)[, 1],
      lat = st_coordinates(centroid)[, 2]
    ) %>%
    filter(
      lat > 33.5,  lat < 34.4,
      lon > -118.68, lon < -118.1
    )
  
  p <- ggplot(block_temp) +
    geom_sf(aes(fill = min_choice_dist), color = NA) +
    geom_sf(data = attendance_zones_ms, color = "white", fill = NA) +
    my_scale +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text = element_blank(),   # remove axis labels
      axis.ticks = element_blank(),   # remove axis ticks
      axis.title = element_blank()    # remove axis titles
    ) +
    labs(fill = "Miles")
  
  
  out_file <- paste0("choice_distance_map_ms_min_choice_dist_", yr, ".pdf")
  ggsave(filename = out_file, plot = p, width = 6, height = 8)
 
}

library(dplyr)
library(tidyr)
# Create a map showing changes 
distance_combined <- distance_dat %>%
  # keep only middle‐school obs in 2004 or 2023
  filter(endyear %in% c(2004, 2023), schooltype == "ms") %>%
  # keep the key vars
  select(censusblockid, year = endyear, min_choice_dist, rel_choice_distance) %>%
  # spread 2004 vs 2023 into separate columns
  pivot_wider(
    names_from   = year,
    values_from  = c(min_choice_dist, rel_choice_distance),
    names_glue   = "{.value}_{year}"
  )
distance_combined$diff <- distance_combined$rel_choice_distance_2023 - distance_combined$rel_choice_distance_2004
block_combined <- blocks20_raw %>%
  left_join(distance_combined, by = c("GEOID" = "censusblockid")) %>%
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  ) %>%
  filter(
    lat > 33.5,  lat < 34.4,
    lon > -118.68, lon < -118.1
  )


my_scale_diff <- scale_fill_gradientn(
  colours = c("#000000", "#800000", "#CD5C5C", "#F5DEB3"),
  na.value = "gray50",
  limits = c(-2.5, 2.5),
  breaks = c(-1.5, -0.5, 0, 0.5, 1.5),
  oob = scales::squish
)

ggplot(block_combined) +
  geom_sf(aes(fill = diff), color = NA) +
  geom_sf(data = attendance_zones_ms, color = "white", fill = NA) +
  my_scale_diff + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text = element_blank(),   # remove axis labels
    axis.ticks = element_blank(),   # remove axis ticks
    axis.title = element_blank()    # remove axis titles
  ) +
  labs(fill = "Miles")
out_file_diff <- "choice_distance_map_ms_min_choice_dist_diff.pdf"
ggsave(filename = out_file_diff, width = 6, height = 8)
