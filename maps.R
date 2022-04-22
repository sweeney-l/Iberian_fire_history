#Sweeney, Harrison and Vander Linden 2022. 
#Assessing anthropogenic influence on fire history during the Holocene in the Iberian Peninsula
# ---------------------------------------------------------

#There are seven seperate scripts associated with this research:
# 1. data_import.R 
# 2. charcoal_analysis.R 
# 3. radiocarbon_analysis.R
# 4. correlation_analysis.R
# 5. SEA_analysis.R
# 6. neolithic_analysis.R (this script)
# 7. maps.R


# ---------------------------------------------------------


# Maps
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script uses the data generated within the
# charcoal_analyis.R and radiocarbon_analysis.R scrips to 
# produce maps of the locations of used radiocarbon and
# charcoal sites. 



# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Load in data and combine
# 3. Map plots 



# ---------------------------------------------------------


# 1. Packages
# ---------------------------------------------------------
# ---------------------------------------------------------

# Install
# ---------------------------------------------------------
# install.packages("c14bazAAR", repos = c(ropensci = "https://ropensci.r-universe.dev"))
# install.packages("changepoint")
# install.packages("fable")
# install.packages("fabletools")
# install.packages("feasts")
# install.packages("ggpubr")
# install.packages("gridExtra")
# install.packages("gstat")
# install.packages("locfit")
# install.packages("maps")
# install.packages("naniar")
# install.packages("rnaturalearth")
# install.packages("raster")
# install.packages("rcarbon")
# install.packages("rgdal")
# install.packages("rgeos")
# install.packages("rio")
# install.packages("sf")
# install.packages("spatstat")
# install.packages("tidyverse")
# install.packages("tseries")
# install.packages("zoo")


# Library
# ---------------------------------------------------------
library(c14bazAAr)
library(changepoint)
library(fable)
library(fabletools)
library(feasts)
library(ggpubr)
library(gridExtra)
library(gstat)
library(locfit)
library(maps)
library(naniar)
library(rnaturalearth)
library(raster)
library(rcarbon)
library(rgdal)
library(rgeos)
library(rio)
library(sf)
library(spatstat)
library(tidyverse)
library(tseries)
library(zoo)


# Paths
# ---------------------------------------------------------
# The script currently runs with the following sub-folders:
# /data: any data that is used for the analysis
# /figs:  figures outputted from the analysis
#   /rpd/cor            : correlation plots
#   /rpd/maps           : maps
#   /pop/SPD/           : SPD plots
#   /pop/SPD/Neo/Kriging: Neolithic surface plots
#   /pop/SPD/Neo/Vario  : Neolithic variogram plots
#   /rpd/rpd            : charcoal composite plots
#   /rpd/neo_rpd        : Neolithic analysis plots
#   /rpd/sea_rpd        : SEA analysis plots
# /other_output         : other types of output
#    /rpd               : intermediate output
#       /rpd_debug      : location for debug files (charcoal)
#       /rpd_influx     : generated influx csvs (charcoal)
#       /z_trans_rpd    : transformed csvs (charcoal)
#       /stats_rpd      : statistics relating to records (charcoal)
#    /pop               : intermediate outpu


# If the folders have not yet been created, the following code
# will generate these from the working directory, with the exception of
# the code file, where this and the other 6 scripts should be saved:

# dir.create("data")
# dir.create("figs")
# dir.create("figs/cor")
# dir.create("figs/maps")
# dir.create("figs/pop")
# dir.create("figs/pop/SPD")
# dir.create("figs/pop/Neo")
# dir.create("figs/pop/Neo/Vario")
# dir.create("figs/pop/Neo/Kriging")
# dir.create("figs/rpd")
# dir.create("figs/sea_rpd")
# dir.create("figs/neo_rpd")
# dir.create("figs/sea_rpd")
# dir.create("other_output")
# dir.create("other_output/pop")
# dir.create("other_output/rpd")
# dir.create("other_output/rpd/rpd_debug")
# dir.create("other_output/rpd/rpd_influx")
# dir.create("other_output/rpd/stats_rpd")
# dir.create("other_output/rpd/z_trans_rpd")

# ---------------------------------------------------------




# RPD data to import
# ---------------------------------------------------------
rpd_map_data <- rio::import("other_output/rpd/rpd_data_input.csv") %>%
  dplyr::group_by(entity_name) %>% 
  dplyr::summarise(longitude = mean(longitude), latitude = mean(latitude))



# Pop data to import (based on output from radiocarbon analysis)
# ---------------------------------------------------------
pop_map_data <- rio::import("./other_output/pop/radiocarbon_data_cal.csv")  %>%
  dplyr::mutate(char_pop = "Archaeological site") %>%
  dplyr::filter(dplyr::between(median_caldates, 3500, 10000))  %>% 
  dplyr::select(longitude, latitude, char_pop) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>% 
  sf::st_transform(crs = 25830) %>% 
  sf::st_transform(crs = 4326)



# ---------------------------------------------------------

# 2. Load in data and combine
# ---------------------------------------------------------
# ---------------------------------------------------------

# Map data
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")      #World map with countries

# Combine rpd and pop data
combined_map_data <- rpd_map_data %>%
  dplyr::mutate(char_pop = "Charcoal site") %>% 
  dplyr::select(longitude, latitude, char_pop) %>%  
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant") %>%   
  rbind(pop_map_data)


# ---------------------------------------------------------

# 3. Map plots 
# ---------------------------------------------------------
# ---------------------------------------------------------
map_char_pop <- ggplot2::ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(text = element_text(size = 10), plot.margin = unit(c(4,1,1,1), "mm")) +
  geom_sf(data = combined_map_data, aes(fill = char_pop, shape = char_pop), size = 2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(21, 23)) +
  coord_sf(xlim = c(-9.7, 3.55), ylim = c(35.95, 43.8), expand = FALSE)+
  labs(fill = "Site type", shape = "Site type")+
  theme(legend.position=c(0.865, 0.07))+
  theme(text = element_text(family = "Times"))+
  theme(legend.title = element_text(size=10), legend.text = element_text(size=8))

ggplot2::ggsave("figs/maps/Iberia_all_char_and_pop.pdf", plot = map_char_pop,  height = 8, width = 11.28, units = "cm")

map_char <- ggplot2::ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(text = element_text(size = 10, family = "Times"), plot.margin = unit(c(4,1,1,1), "mm")) +
  geom_sf(data = dplyr::filter(combined_map_data, char_pop == "Charcoal site"), aes(fill = char_pop, shape = char_pop), size = 2) +
  scale_fill_manual(values = c("blue")) +
  scale_shape_manual(values = c(23)) +
  coord_sf(xlim = c(-9.7, 3.55), ylim = c(35.95, 43.8), expand = FALSE)+
  labs(x = "Longitude", y = "Latitude") + 
  annotate(geom = "text", label = "A", x = -9.3, y = 43.5, family = "Times", fontface = "bold")+
  guides(fill = F, shape = F)

ggplot2::ggsave("figs/maps/Iberia_all_char.pdf", plot = map_char,  height = 8, width = 11.28, units = "cm")

map_pop <-  ggplot2::ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  theme(text = element_text(size = 10, family = "Times"), plot.margin = unit(c(4,1,1,1), "mm")) +
  geom_sf(data = dplyr::filter(combined_map_data, char_pop == "Archaeological site"), aes(fill = char_pop, shape = char_pop), size = 2) +
  scale_fill_manual(values = c("red")) +
  scale_shape_manual(values = c(23)) +
  coord_sf(xlim = c(-9.7, 3.55), ylim = c(35.95, 43.8), expand = FALSE)+
  labs(x = "Longitude", y = "Latitude") + 
  annotate(geom = "text", label = "B", x = -9.3, y = 43.5, family = "Times", fontface = "bold")+
  guides(fill = F, shape = F)

ggplot2::ggsave("figs/maps/Iberia_all_pop.pdf", plot = map_pop, height = 8, width = 11.28, units = "cm")


maps_char_pop <- gridExtra::grid.arrange(map_char, map_pop, nrow = 2)

ggplot2::ggsave("figs/maps/Iberia_all_pop_char_dual.pdf", plot = maps_char_pop, height = 16, width = 11.28, units = "cm")


