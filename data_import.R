#Sweeney, Harrison and Vander Linden 2021. 
#Assessing anthropogenic influence on fire history during the Holocene in the Iberian Peninsula
# ---------------------------------------------------------

#There are seven seperate scripts associated with this research:
# 1. data_import.R (this script)
# 2. charcoal_analysis.R 
# 3. radiocarbon_analysis.R
# 4. correlation_analysis.R
# 5. SEA_analysis.R
# 6. neolithic_analysis.R
# 7. maps.R


# ---------------------------------------------------------



# RPD charcoal data and radiocarbon data import
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script details how to import and treat the data
# used in subsequent analysis




# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Import and treat charcoal data
# 3. Import and treat radiocarbon data 

# ---------------------------------------------------------


# 1. Packages
# ---------------------------------------------------------
# ---------------------------------------------------------

# Install
# ---------------------------------------------------------
# install.packages("c14bazAAr")
# install.packages("changepoint")
# install.packages("fable")
# install.packages("fabletools")
# install.packages("gridExtra")
# install.packages("gstat")
# install.packages("locfit")
# install.packages("maps")
# install.packages("naniar")
# install.packages("rnaturalearth")
# install.packages("raster")
# install.packages("rcarbon")
# install.packages("rgeos")
# install.packages("rio")
# install.packages("sf")
# install.packages("spatstat")
# install.packages("tidyverse")
# install.packages("zoo")


# Library
# ---------------------------------------------------------
library(changepoint)
library(fable)
library(fabletools)
library(gridExtra)
library(gstat)
library(locfit)
library(maps)
library(naniar)
library(rnaturalearth)
library(raster)
library(rcarbon)
library(rgeos)
library(rio)
library(sf)
library(spatstat)
library(tidyverse)
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
#    /pop               : intermediate output


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




# 2. Import and treat charcoal data
# ---------------------------------------------------------
# ---------------------------------------------------------

# Import data from Iberia database https://researchdata.reading.ac.uk/294/ and save csv in datafile to file data/
# iberia_charcoal_records <- rio::import("data/iberia_charcoal_records.csv")
iberia_charcoal_records <- rio::import("data/iberia_charcoal_records.xlsx") %>%   #temp using data from Roberto
  dplyr::filter(!is.na(INTCAL2020_median))

# Filter and treat data 
# ---------------------------------------------------------
minmax_entity <- iberia_charcoal_records %>% #Establish first and last date per record
  dplyr::group_by(entity_name) %>% 
  dplyr::summarise(max_age = max(INTCAL2020_median), min_age = min(INTCAL2020_median))

char_data <- iberia_charcoal_records %>% 
  dplyr::left_join(minmax_entity, by = "entity_name") %>% 
  dplyr::filter(max_age >= 3500 & min_age < 10000) %>% #Filter to records that provide data within age range
  dplyr::filter(!entity_name %in% c("PRD1", "PRD2", "PRD3", "PRD5")) %>%    #Remove multiple cores from same site
  dplyr::filter(!entity_name %in% c("Pena Negra core", 
                                    "Vilamora P01-5_100minus", 
                                    "Vilamora P01-5_100plus", 
                                    "Besos core Riera-Mora",
                                    "Las Vinuelas core_micro")) %>%    #Remove records with limited data (less than 5 samples)
  dplyr::filter(!entity_name %in% c("El Brezosa core_micro",
                                    "Hoya del Castillo N-CAS macro",
                                    "N-GUA_micro",
                                    "PORTALET_micro")) %>%  #Remove macro_micro records with less sample data than macro or micro equivalent
  dplyr::mutate(macro_micro = dplyr::if_else(entity_name %in% c("Abi 05_07_100minus",
                                                                "ADP 01_06_100minus",
                                                                "Baza section",
                                                                "BSM08",
                                                                "Canada de la Cruz core",
                                                                "Candieira (Charco da Candieira) core",
                                                                "Cha das Lameiras soil profile_micro",
                                                                "El Carrizal core",
                                                                "Espinosa de Cerrato core",
                                                                "Gador core",
                                                                "Hinojos Marsh_core S1_micro",
                                                                "Hoya del Castillo N-CAS",
                                                                "Marbore composite",
                                                                "NAVARRE3",
                                                                "Ojos del Tremendal core",
                                                                "Pena da Cadela core",
                                                                "PRD4",
                                                                "Siles Lake core",
                                                                "Tubilla del Lago core",
                                                                "VdL PB2_100minus",
                                                                "Villaverde core"),
                                             "micro", "macro")) %>%  #Specify whether records are macro or micro charcoal
  dplyr::mutate(multiple_entity_site = dplyr::if_else(site_name %in% c("Alvor Estuary Ribeira do Farelo Ribeira da Torre",
                                                                       "Armacao de Pera Ribeira de Alcantarilha",
                                                                       "Cha das Lameiras",
                                                                       "Hinojos Marsh", 
                                                                       "Valle do Lobo Ribeira de Carcavai"),
                                                      "y", "n")) %>%  #For remaining sites, identify records from same site but micro or macro
  dplyr::group_by(site_name) %>% #Add unique site ID
  dplyr::mutate(ID_SITE = dplyr::cur_group_id()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(entity_name) %>% #Add unique entity ID
  dplyr::mutate(ID_ENTITY = dplyr::cur_group_id()) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(sample_depth = avg_depth..cm.) %>% #Change variable names
  dplyr::rename(quantity = charcoal.quantity) %>% #Change variable names
  dplyr::mutate(ID_SAMPLE = dplyr::row_number()) #Add sample ID

# Exclude certain records (nb to generate different composite curves, change these)
rpd_data_input <- char_data %>% 
  dplyr::filter(elevation < 2000) %>% #Filter sites above a certain elevation
  # dplyr::filter(elevation > 2000) %>% #Select sites above a certain elevation
  dplyr::filter(multiple_entity_site == "n" | (multiple_entity_site == "y" & macro_micro == "macro")) #Filter micro sites if multi record site
  # dplyr::filter(multiple_entity_site == "n" | (multiple_entity_site == "y" & macro_micro == "micro")) #Filter macro sites if multi record site

length(unique(rpd_data_input$entity_name)) #The number of individual records
rio::export(rpd_data_input, "other_output/rpd/rpd_data_input.csv") #Export data for subsequent analysis

# ---------------------------------------------------------




# 3. Import and treat radiocarbon data 
# ---------------------------------------------------------
# ---------------------------------------------------------

# The radiocarbon dataset used for this analysis has been 
# uploaded to http://dx.doi.org/10.17864/1947.000340. To replicate the analysis, import data
# from Iberia database and save csv in datafile to file
# data/
radiocarbon_data <- rio::import("other_output/pop/radiocarbon_data.csv")






# The following code was used to generate this dataset.

# Radiocarbon data from the following datasets:
#   Balsera et al. (2015) - https://doi.org/10.1016/j.quaint.2015.06.022
#   Drake et al. (2017) - https://doi.org/10.1007/s10816-016-9286-y
#   McLaughlin et al. (2021) - https://doi.org/10.1098/rstb.2019.0724
#   Pardo-Gordó et al. (2019) - https://doi.org/10.5334/joad.49
# and through R package c14bazAAr (Schmid et al., 2019)
#   Capuzzo et al. (2014) - https://doi.org/10.2458/56.17453
#   d’Errico et al. (2011) - https://doi.org/doi:10.4207/PA.2011.ART40
#   Hinz et al. (2012) - https://doi.org/10.12766/jna.2012.65
#   Kniesel et al. (2014)  - no DOI. see http://radon-b.ufg.uni-kiel.de.
#   Manning et al. (2016) - https://doi.org/10.1016/j.quascirev.2014.07.003
#   Vermeersch (2020) - https://doi.org/10.1016/j.dib.2020.105793

# GTOPO30 Elevation maps from https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-30-arc-second-elevation-gtopo30?qt-science_center_objects=0#qt-science_center_objects
# Save maps gt30w020n40 and gt30w020n90 in datapath

#First download and save the data from each source in csv to data path

# Import the raw data
balsera_data_raw <- rio::import("data/Balsera_dataset_v4.csv") #NB this data is available as pdf, so needs to be converted first to csv. Note, be careful with the conversion of long lats - ensure these are numeric in csv
drake_data_raw <- rio::import("data/Drake_dataset.csv")
pardo_data_raw <- rio::import("data/Pardo_dataset.csv")
# c14baz_data_raw <- c14bazAAR::get_c14data(c("eubar", "euroevol", "radon", "radonb", "pacea", "14cpalaeolithic")) #import data and then save
# rio::export(c14baz_data_raw, "data/c14baz_dataset.csv") # Saved file so that don't have to keep downloading
c14baz_data_raw <- rio::import("data/c14baz_dataset.csv")
mclaughlin_data_raw<- rio::import("data/McLaughlin_dataset.xlsx")

# Import elevation information
box <- raster::extent(-10, 5, 35, 45)
elevation_map1 <- raster::raster("data/gt30w020n40.tif")
elevation_map2 <- raster::raster("data/gt30w020n90.tif")
elevation_map <- raster::merge(elevation_map1, elevation_map2) %>% #Merge maps and reduce to Iberia
  raster::crop(box)
elevation_point <- raster::rasterToPoints(elevation_map)  # Get point elevation data for later filtering
elevations <- data.frame(elevation_point)
colnames(elevations) = c("lon", "lat", "asl")  # Rename output

# Generate owin window
iberia_map <- maps::map("world", c("Spain(?!:)", "Portugal", "Andorra"), fill = T, resolution = 0) #Download map of Iberia
iberia_map <- sf::st_as_sf(iberia_map) #Convert to sf
iberia_map_25830 <- sf::st_transform(iberia_map, crs = 25830) #Change the crs to flat projection
iberia_owin_25830 <- spatstat.geom::as.owin(iberia_map_25830) #Generate owin
iberia_owin_25830_exp <- spatstat.core::expand.owin(iberia_owin_25830, distance = 2000) #Expand owin to try and ensure that samples are not excluded due to coord inaccuracies 
iberia_map_exp_25830 <-  sf::st_as_sf(iberia_owin_25830_exp) #Convert back to sf
sf::st_crs(iberia_map_exp_25830) <- 25830 #Set crs
iberia_map_exp <- sf::st_transform(iberia_map_exp_25830, crs = 4258) #Change crs back to ellipsoid for raster



# Clean data/ standardise
# ---------------------------------------------------------
balsera_data <- balsera_data_raw %>% 
  dplyr::arrange(Date, SD, Site, `Lab #`) %>% 
  dplyr::rename(c14_age = Date, c14_std = SD, site = Site, material = Material, longitude = Lon_ETRS89, latitude = Lat_ETRS89) %>%  
  dplyr::mutate(lab = stringr::word(`Lab #`, 1, sep = "[:digit:]")) %>% #Set up lab numbers
  dplyr::mutate(lab = stringr::word(lab, 1, sep = "\\/")) %>% #Set up lab numbers
  dplyr::mutate(lab1 = stringr::str_extract(`Lab #`, "\\d+\\.*\\d*")) %>% #Set up lab numbers
  dplyr::mutate(lab_nr = dplyr::if_else(is.na(lab1), "Unreported", paste0(lab, "-", lab1))) %>%  #Set up lab numbers
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "\\*", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = trimws(lab_nr)) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "--", "-")) %>%
  dplyr::mutate(db_code = paste0("b", row_number())) %>%  #Add unique item code
  dplyr::mutate(period = "") %>% #Add cultural period variable
  dplyr::mutate(source_db = "balsera") %>% 
  dplyr::mutate(longitude = as.numeric(longitude)) %>% 
  dplyr::mutate(latitude = as.numeric(latitude)) %>% 
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Circe", "DSA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "BIRM", "Birm")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "DHS", "DSA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gak", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GaK", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gx", "GX")) %>%#Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "OxAV", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Rome", "R")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UCIAMS", "UCI")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGA", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAMS", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGaMS", "UGa")) %>%   #Clean lab codes
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "CSIC-626" & c14_age == 4630, "CSIC-626A", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "CSIC-626" & c14_age == 4630, "CSIC-626B", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "CSIC-637" & c14_age == 7200, "CSIC-637R", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "I-16079" & c14_age == 4660, "I-16079C", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "MC-1112" & c14_age == 4700, "MC-1112B", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "OxA-2360" & c14_age == 3946, "OxA-2360-15", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "OxA-2360" & c14_age == 4062, "OxA-2360-23", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "OxA-2360" & c14_age == 6399, "OxA-2360-25", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "UBAR-464" & c14_age == 3450, "UBAR-464B", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "UBAR-464" & c14_age == 3570, "UBAR-464A", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "UGRA-76" & c14_age == 3890, "UGRA-76B", lab_nr)) %>% #Different lab nr
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "Gif-8079", "Gif-8079A", lab_nr)) %>% #Same lab nr as another, both look correct, so adjust
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "ICEN-76" & c14_age == 7810, "ICEN-76A", lab_nr)) %>% #Same lab nr as another, both look correct, so adjust
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "GifA-99112" & c14_age == 5480, "GifA-99113", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "GifA-99112" & c14_age == 5580, "GifA-99114", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "GrN-7009" & c14_age == 3980, "GrN-7008", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "ICEN-1159" & c14_age == 4460, "ICEN-1149", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "IRPA-1062" & c14_age == 3390, "IRPA-1063", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "KIK-3487" & c14_age == 4715, "KIA-37691", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "Ua-35665" & c14_age == 3830, "Ua-35655", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "Ua-37894" & c14_age == 3480, "Ua-37895", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "Ua-40763" & c14_age == 3681, "Ua-40761", lab_nr)) %>% #Error in database following check
  dplyr::mutate(lab_nr = dplyr::if_else(lab_nr == "Ua-4820", "Ua-4821", lab_nr)) %>% #Error in database following check
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "GAK-8959" & c14_age == 6480, 180L, c14_std)) %>% #Error in database following check
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "Ua-35665" & c14_age == 4370, 35L, c14_std)) %>% #Error in database following check
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "UBAR-274", 50L, c14_std)) %>% #Error in database following check
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "Wk-25162", 30L, c14_std)) %>% #Error in database following check
  naniar::replace_with_na(replace = list(lab_nr = c("KNI-200", "KNJ-115", "KNJ-117"))) %>% #Clean lab codes
  dplyr::mutate(longitude = if_else(longitude == 1.838611111, -1.838611111, longitude)) %>% #Error in database following check
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

c14baz_data <- c14baz_data_raw %>% 
  dplyr::mutate(db_code = paste0("c", row_number())) %>%   #Add unique item code
  dplyr::filter(!(is.na(lat) | is.na(lon))) %>% #Remove missing coordinate entries
  dplyr::rename(c14_age = c14age, c14_std = c14std, source_db = sourcedb, lab_nr = labnr) %>% 
  dplyr::mutate(longitude = lon) %>% 
  dplyr::mutate(latitude = lat) %>%
  dplyr::mutate(country = maps::map.where(x =  longitude, y = latitude)) %>%   
  dplyr::filter(country %in% c("Spain","Portugal", "Andorra") | is.na(country)) %>% #Filter to area, or where uncertain
  dplyr::filter(between(latitude, 36, 45)) %>% 
  dplyr::filter(between(longitude, -10, 5)) %>%
  dplyr::filter(c14_age <=20000) %>% #Limit c14 date range
  dplyr::mutate(longitude = dplyr::if_else(lab_nr == "Beta-297104", -4.687220, longitude)) %>% #Wrong longitude, checked
  dplyr::mutate(longitude = dplyr::if_else(site == "Casinha Derribada", -7.863333, longitude)) %>% #Wrong longitude, checked
  dplyr::mutate(longitude = dplyr::if_else(site == "Casinha Derribada 3", -7.863333, longitude)) %>% #Wrong longitude, checked
  dplyr::mutate(latitude = dplyr::if_else(site == "Casinha Derribada", 40.71944, latitude)) %>% #Wrong latitude, checked
  dplyr::mutate(latitude = dplyr::if_else(site == "Casinha Derribada 3", 40.71944, latitude)) %>% #Wrong latitude, checked
  dplyr::filter(lab_nr != "Ly-706") %>% #In France
  dplyr::filter(site != "Montpellier Richemont") %>% #In France
  dplyr::filter(site != "La Fangade") %>% #In France
  dplyr::mutate(material = dplyr::if_else(lab_nr == "AA-57439", "bone / bos", material)) %>% #Duplicate code wwith Balsera, but different ste. This is correct, use matreial details from Balsera
  dplyr::mutate(material = dplyr::if_else(lab_nr == "Beta-166228", "carbon", material)) %>% #Duplicate code with Pardo, but different ste. This is correct, use matreial details from Pardo
  dplyr::mutate(c14_age = dplyr::if_else(lab_nr == "OxA-5508", 4050L, c14_age)) %>% #Error in database
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "OxA-5508", 60L, c14_std)) %>% #Error in database
  dplyr::mutate(c14_std = dplyr::if_else(lab_nr == "UGRA-185", 120L, c14_std)) %>% #Error in database
  dplyr::mutate(c14_age = dplyr::if_else(lab_nr == "UGRA-185", 3930L, c14_age)) %>% #Error in database
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "\\*", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, " ", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "AA-8647.T.461", "AA-8647")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Bera", "Beta")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Bea", "Beta")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "BIRM", "Birm")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "BM-1453/I-10736", "BM-1453")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "CAMP", "CAMS")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "CAMS-9918/Beta-67949", "CAMS-9918")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Col", "COL")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "EHT", "ETH")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gak", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GaK", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "COL2", "COL-2")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "COL4", "COL-4")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GAK-6448/6460", "GAK-6448")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gra", "GrA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GdrA", "GrA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GRN", "GrN")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gx", "GX")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "HD", "Hd")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Kia", "KIA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Kn", "KN")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Lod", "LOD")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "LYON", "Ly")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Lyon", "Ly")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Ly-49(OxA)", "Ly-49")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "NzA", "NZA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "OA", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Oxa", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "OXA", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Pioz", "Poz")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Pz", "Poz")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "PxA", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Rome", "R")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Sanu", "SANU")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "SU", "Su")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "SuA", "SUA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGA", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGaMS", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAMS", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAM", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "WK", "Wk")) %>% #Clean lab codes
  naniar::replace_with_na(replace = list(lab_nr = c("?", "No det.", "I-?", #Clean lab codes
                                                    "8649X-404", "1Tucson", "Ariz.Univ",
                                                    "Beta-", "Gif-?", "Hv-?",
                                                    "Gif-A-II.4", "Gif-", "GrN-",
                                                    "LUNK-0725", "LUNK-0762", "LUNK-0859",
                                                    "LUNK-0860", "LUNK-0867", "LUNK-0867",
                                                    "n/a-n/a", "Nodet.", "OxA-",
                                                    "SMU-?", "Ua-", "Ua-CabezoJuré"))) %>%
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code) %>%
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db) 


drake_data <- drake_data_raw %>% 
  dplyr::mutate(db_code = paste0("d", row_number())) %>% #Add unique item code
  dplyr::mutate(source_db = "drake") %>% 
  dplyr::rename(c14_age = "14C Year", c14_std = Sigma, lab_nr = Sample, site = Site, material = Material) %>% 
  dplyr::mutate(longitude = Longitude) %>% 
  dplyr::mutate(latitude = Latitude) %>% 
  dplyr::mutate(period = "") %>%
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "\\*", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "=", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = trimws(lab_nr)) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, " ", "-")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "CSIC-150B", "CSIC-150")) %>% #Based on duplicate Balsera
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "DSH", "DSA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGA", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Mams", "MAMS"))  %>% #Clean lab codes
  naniar::replace_with_na(replace = list(lab_nr = c("0", "", "?-17697", "UNK-1"))) %>% #Clean lab codes
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db) 

pardo_data <- pardo_data_raw %>%
  dplyr::mutate(db_code = paste0("p", row_number())) %>% #Add unique item code
  dplyr::mutate(source_db = "pardo") %>% 
  dplyr::rename(c14_age = "FechaBP", c14_std = Desviación, lab_nr = IDmuestra, site = Yacimiento, material = Material, period = "Adscripcion cronológica") %>% 
  dplyr::mutate(longitude = Long) %>% 
  dplyr::mutate(latitude = Lat) %>% 
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "\\*", "")) %>% #Tidy up lab number code
  dplyr::mutate(lab_nr = trimws(lab_nr)) %>% #Tidy up lab number code
  dplyr::mutate(c14_age = dplyr::if_else(db_code == "p598", 5170L, c14_age)) %>% #Wrong date
  dplyr::mutate(c14_age = dplyr::if_else(db_code == "p319", 6220L, c14_age)) %>% #Wrong date
  dplyr::mutate(material = dplyr::if_else(db_code == "p371", "Human bone", material)) %>% #Wrong material
  dplyr::mutate(lab_nr = dplyr::if_else(db_code == "p538", "OxA-2360-25", lab_nr)) %>% #Prob should be this rather than OxA-236025 in line with balsera
  dplyr::mutate(lab_nr = dplyr::if_else(db_code == "p1539", "Wk-27462", lab_nr)) %>% #Checked, wrong lab-nr for this item
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "AA-15499", "AA-16499")) %>% #Checked and should be same as Balsera 16499 (15499 not found)
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Beta-59996", "Beta-59998")) %>% #Checked and should be same as Drake 59998 (59996 not found)
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Beta-166232", "Beta-166231")) %>% #Checked and should be same as Drake 166231 (166232 not found)
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Beta-89225", "Beta-89285")) %>% #Checked and should be same as Balsera 89285 (89225 not found)
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "BLN", "Bln")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Cams", "CAMS")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "CCDS-1", "GAK-8950")) %>% #Based on manual checck of data with c14baz
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Col", "COL")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "DSH", "DSA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gak", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GaK", "GAK")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GiF/LSM-11037", "Gif/LSM-11037")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Gr-15369", "GrA-15369")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GRN", "GrN")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "GrN-19282", "GrN-18282")) %>% #Checked and should be same as Balsera 18282 (92 relates to record in Ireland)
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Mams", "MAMS")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Oxa", "OxA")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGA", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAMS", "UGa")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAMS", "UGaMS")) %>% #Clean lab codes
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UtC-9329", "UtC-9320")) %>% #Checked and should be same as Balsera 9320 (9329 not found)
  dplyr::mutate(c14_std = dplyr::if_else(db_code == "p299", 55L, c14_std)) %>% #Error in database
  dplyr::mutate(c14_std = dplyr::if_else(db_code == "p239", 75L, c14_std)) %>% #Error in database
  naniar::replace_with_na(replace = list(lab_nr = c("rev-"))) %>% #Clean lab codes
  dplyr::mutate(longitude = dplyr::if_else(longitude<(-1000), longitude/1000000, longitude)) %>% #Correct errors in dataset
  dplyr::mutate(latitude = dplyr::if_else(latitude >1000, latitude/1000000, latitude)) %>% #Correct errors in dataset
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period,  db_code) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

mclaughlin_data_4258 <- mclaughlin_data_raw %>% 
  dplyr::mutate(db_code = paste0("m", row_number())) %>% #Add unique item code
  dplyr::rename(period = Cultural.attribution) %>% 
  dplyr::rename(site = Site) %>% 
  dplyr::rename(lab_nr = LabCode) %>% 
  dplyr::rename(c14_age = BP) %>% 
  dplyr::rename(c14_std = Std) %>% 
  dplyr::rename(material = Sample) %>% 
  dplyr::mutate(source_db = "mclaughlin") %>%  
  dplyr::mutate(longitude = UTM_E) %>% 
  dplyr::mutate(latitude = UTM_N) %>% 
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "Q-(AM85B2b)", "QAM-85B2b")) %>% #To standardise
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "ICEN-889", "ICEN-899")) %>% #Checked and should be same as Pardo 158 
  dplyr::mutate(lab_nr = stringr::str_replace_all(lab_nr, "UGAMS", "UGa")) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 25830, agr = "constant") %>%
  sf::st_transform(., crs = 4258) 

mclaughlin_data_lon_lat <- sf::st_coordinates(mclaughlin_data_4258) %>% 
  dplyr::as_tibble()

mclaughlin_data <- mclaughlin_data_4258 %>% #Convert to lon lat projection
  sf::st_set_geometry(NULL) %>% 
  cbind(mclaughlin_data_lon_lat) %>% 
  dplyr::rename(longitude = X, latitude = Y) %>% 
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

data_order <- c("balsera", 
                "pardo",
                "mclaughlin",
                "drake", 
                "bird", 
                "radon", 
                "radonb", 
                "eubar", 
                "euroevol", 
                "14cpalaeolithic",
                "pacea")

# Combine datasets
# ---------------------------------------------------------

# Combine data, check coordinates, search for missing labs and replace with "unknown_#"
c14_data <- dplyr::bind_rows(balsera_data, mclaughlin_data, drake_data, pardo_data, c14baz_data, .id = NULL) %>% 
  dplyr::mutate(lab_check = stringr::str_count(lab_nr, "\\-")) %>%  #ID lab nrs without lab and nr format
  dplyr::mutate(lab_nr = ifelse(lab_check < 1, NA, lab_nr)) %>% #Make incomplete lab nrs. NA
  dplyr::mutate(source_db = factor(source_db, levels = data_order)) %>% 
  dplyr::filter(c14_std <= 200) %>% #Filter out SD greater than 200
  dplyr::filter(!(longitude > 2 & latitude < 41)) %>% #Remove data from Balearics
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

c14_data_lon_lat <- c14_data %>% #Lon lats of data
  dplyr::select(longitude, latitude) %>% 
  as.matrix()

c14_elevation <- raster::extract(elevation_map, c14_data_lon_lat)  # Using topo raster, extract elevations for pop data

combined_data <- c14_data %>% 
  dplyr::mutate(elevation = c14_elevation) %>% 
  dplyr::mutate(elevation = dplyr::if_else(is.na(elevation), 0, elevation)) %>% 
  dplyr::filter(elevation <= 2000) #Limit to less than 2000m

na_lab_count <- sum(is.na(combined_data$lab_nr)) #How many samples do not have complete lab codes
unknown_labs <- ifelse(is.na(combined_data$lab_nr),  "unknown-", combined_data$lab_nr) %>% #Generate tibble with NA converted to numbered unknown no.
  dplyr::as_tibble() %>%
  dplyr::rename(lab_nr = value) %>% 
  dplyr::mutate(unknown_nr = row_number()) %>% 
  dplyr::mutate(lab_nr = ifelse(lab_nr == "unknown-", paste0("unknown-", unknown_nr), lab_nr)) 

population_data <- combined_data %>% 
  dplyr::mutate(lab_nr = coalesce(combined_data$lab_nr, unknown_labs$lab_nr)) %>%  #Replace NA values with unknown lab values
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db) %>% 
  dplyr::mutate(id = row_number())

# Make sure that period tags are included where duplicate lab nrs are in the data
period_population_data <- population_data %>%  #Extract period data from other potentially excluded
  dplyr::filter(period != "")

dist_period_population_data <- dplyr::distinct(period_population_data, lab_nr, .keep_all = T)

no_period_population_data <- population_data %>%  
  dplyr::filter(period == "")  %>%  
  dplyr::left_join(dplyr::select(dist_period_population_data, period, lab_nr), by = "lab_nr") %>%  
  dplyr::mutate(period = period.y) %>% 
  dplyr::select(-period.x, -period.y) %>% 
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code, lab_check, elevation, id)

pop_data_tag <- period_population_data %>% 
  rbind(no_period_population_data) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

lab_codes <- unique(sub("-.*", "", pop_data_tag$lab_nr)) %>% 
  as_tibble()

# Ensure material information isn't lost where dublicates are in the data
material_population_data <- pop_data_tag %>%  #Extract material data from other potentially excluded
  naniar::replace_with_na(replace = list(material = c("Unreported", "Unreported / \"\"bulk sediment\"\"",
                                                      "n.d", "Otros","", "miscellaneous", "No det.", "from level 21 ivd",
                                                      "?", "from layer 7 rosa", "layer ?", "from -2.19m", "possibly intrusive (?)from B2-B3",
                                                      "from layer IB", "from layer IIB", "from surface near art 188, Puenta",
                                                      "global sample from IV", "line on panl nø 14", "from level 1c1",
                                                      "from leel 1be", "line below a deer", "from level 1c2", "PAI B12.285, layer IV",
                                                      ".", "from Er-1790", "from level 1c3", "e Pai F10.269, layer 4", 
                                                      "from Er-30-75", "triangular sign nr 57A", "PAI F10.269, layer 4", "from layer IV",
                                                      "from PAI F14.296, layer 4", "from the base", "Hogar from level 8/s", "possible from level III"))) %>%
  dplyr::filter(!is.na(material)) 

dist_material_population_data <- dplyr::distinct(material_population_data, lab_nr, .keep_all = T) #Needed otherwise duplicate lab_nrs can end up adding rows

no_material_population_data <- pop_data_tag %>% 
  naniar::replace_with_na(replace = list(material = c("Unreported", "Unreported / \"\"bulk sediment\"\"",
                                                      "n.d", "Otros","", "miscellaneous", "No det.", "from level 21 ivd",
                                                      "?", "from layer 7 rosa", "layer ?", "from -2.19m", "possibly intrusive (?)from B2-B3",
                                                      "from layer IB", "from layer IIB", "from surface near art 188, Puenta",
                                                      "global sample from IV", "line on panl nø 14", "from level 1c1",
                                                      "from leel 1be", "line below a deer", "from level 1c2", "PAI B12.285, layer IV",
                                                      ".", "from Er-1790", "from level 1c3", "e Pai F10.269, layer 4", 
                                                      "from Er-30-75", "triangular sign nr 57A", "PAI F10.269, layer 4", "from layer IV",
                                                      "from PAI F14.296, layer 4", "from the base", "Hogar from level 8/s", "possible from level III"))) %>%
  dplyr::filter(is.na(material)) %>%    
  dplyr::left_join(dplyr::select(dist_material_population_data, material, lab_nr), by = "lab_nr") %>% 
  dplyr::mutate(material = material.y) %>% 
  dplyr::select(-material.x, -material.y) %>% 
  dplyr::select(c14_age, c14_std, source_db, longitude, latitude, lab_nr, site, material, period, db_code, lab_check, elevation, id)

pop_data <- material_population_data %>% 
  rbind(no_material_population_data) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)


# Test for duplicates
# ---------------------------------------------------------
# First pass: lab nr, age and error
duplicate_testing_lab <- vector("integer", nrow(pop_data))

d1_pop_data <- pop_data %>% #Ensure data in the right order for testing for duplicates
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

for (i in 2:nrow(d1_pop_data)){#Generate vector to see whether duplicate rows based on lab_nr
  duplicate_testing_lab [[i]] <- if_else((d1_pop_data$lab_nr[[i]] == d1_pop_data$lab_nr[[(i-1)]]
                                          & d1_pop_data$c14_age[[i]] == d1_pop_data$c14_age[[(i-1)]] 
                                          & d1_pop_data$c14_std[[i]] == d1_pop_data$c14_std[[(i-1)]])
                                         , 1, 0 )
}

pop_data_check_dr1 <- d1_pop_data %>% #insert duplicates, visually check this!
  dplyr::mutate(duplicate = duplicate_testing_lab)

pop_data_dr1 <- pop_data_check_dr1 %>%  #Remove duplicate labs
  dplyr::filter(duplicate == 0) %>% 
  dplyr::select(-duplicate) %>% 
  dplyr::filter(c14_std > 0) %>% 
  dplyr::mutate(latitude2 = round(latitude, 2)) %>% 
  dplyr::mutate(longitude2 = round(longitude, 2))

# 2nd pass: coordinates, age and error
duplicate_testing_coords <- vector("integer", nrow(pop_data_dr1))

d2_pop_data <- pop_data_dr1 %>% #Ensure data in the right order for testing for duplicates
  dplyr::arrange(latitude2, longitude2, c14_age, c14_std, lab_nr, source_db)

for (i in 2:nrow(d2_pop_data)){#Generate vector to see whether duplicate rows based on coords
  duplicate_testing_coords [[i]] <- if_else((d2_pop_data$longitude2[[i]] == d2_pop_data$longitude2[[(i-1)]]
                                             & d2_pop_data$latitude2[[i]] == d2_pop_data$latitude2[[(i-1)]]   
                                             & d2_pop_data$c14_age[[i]] == d2_pop_data$c14_age[[(i-1)]] 
                                             & d2_pop_data$c14_std[[i]] == d2_pop_data$c14_std[[(i-1)]]
                                             & (is.na(d2_pop_data$lab_check[[i]]) | is.na(d2_pop_data$lab_check[[i-1]])))
                                            , 1, 0 )
}

pop_data_check_dr2 <- d2_pop_data %>% #insert duplicates, visually check this!
  dplyr::mutate(duplicate = duplicate_testing_coords)

pop_data_dr2 <- pop_data_check_dr2 %>%  #Remove duplicate coords, etc.
  dplyr::filter(duplicate == 0) %>% 
  dplyr::select(-duplicate) 

# 3rd pass: error, visual check
duplicate_testing_error <- vector("integer", nrow(pop_data_dr2))

d3_pop_data <- pop_data_dr2 %>% #Ensure data in the right order for testing for duplicates
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)

for (i in 2:nrow(d3_pop_data)){#Generate vector to see whether duplicate rows based on std and age only
  duplicate_testing_error [[i]] <- if_else((d3_pop_data$c14_std[[i]] == d3_pop_data$c14_std[[(i-1)]]
                                            & d3_pop_data$c14_age[[i]] == d3_pop_data$c14_age[[(i-1)]])
                                           , 1, 0 )
}

pop_data_check_dr3 <- d3_pop_data %>% #insert duplicates, visually check this!
  dplyr::mutate(duplicate = duplicate_testing_error)

pop_data_dr3 <- pop_data_check_dr3 %>%  #Remove duplicate coords, etc.
  dplyr::filter(!(lab_nr == "CSIC-201B" & c14_age == 2570 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "KN-I.-201" & c14_age == 2770 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "IAB-110" & c14_age == 2980 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "Ua-2310" & c14_age == 3000 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "ICEN-84" & c14_age == 3000 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-534" & c14_age == 3080 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "unknown-647" & c14_age == 3170 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "unknown-671" & c14_age == 3180 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "P-4069A" & c14_age == 3200 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "I-7706" & c14_age == 3230 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "unknown-768" & c14_age == 3245 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UtC-1436" & c14_age == 3280 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Ua-394892" & c14_age == 3287 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-1032" & c14_age == 3330 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-1150" & c14_age == 3365 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "UtC-1354" & c14_age == 3370 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UBAR-673" & c14_age == 3370 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "ICEN-11069" & c14_age == 3470 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-1572" & c14_age == 3475 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "UtC-1433" & c14_age == 3480 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-1625" & c14_age == 3490 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-11070" & c14_age == 3520 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UtC-2289" & c14_age == 3530 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UtC-1437" & c14_age == 3530 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-91583" & c14_age == 3570 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "KIA-22256" & c14_age == 3580 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UtC-1439" & c14_age == 3580 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UBAR-661" & c14_age == 3630 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "KIA-22257" & c14_age == 3630 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "ICEN-16352" & c14_age == 3640 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "H-2048/1458" & c14_age == 3650 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "C-341" & c14_age == 3680 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "B-6331" & c14_age == 3680 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "ICEN-16063" & c14_age == 3680 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UtC-2284" & c14_age == 3700 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Ua-26013" & c14_age == 3705 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UtC-2292" & c14_age == 3720 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "I-19306" & c14_age == 3830 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "I-15319" & c14_age == 3870 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UBAR-297" & c14_age == 3890 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "CSIC-30647" & c14_age == 3900 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "OxA-236015" & c14_age == 3946 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-3226" & c14_age == 3950 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-7007C" & c14_age == 3950 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GrN-15511" & c14_age == 3990 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "unknown-3373" & c14_age == 3990 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-3405" & c14_age == 3995 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-125862" & c14_age == 4000 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-601" & c14_age == 4010 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "UtC-2630" & c14_age == 4040 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "OxA-236023" & c14_age == 4062 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-193744" & c14_age == 4130 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "AA-4238" & c14_age == 4220 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GAK-10943" & c14_age == 4220 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "H-2049/148" & c14_age == 4260 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "H-204-247" & c14_age == 4295 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "unknown-4382" & c14_age == 4295 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-4610" & c14_age == 4420 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "R-1768" & c14_age == 4515 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-172" & c14_age == 4540 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-5506" & c14_age == 4600 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "MC-1112B" & c14_age == 4600 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "unknown-5169" & c14_age == 4790 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Gif-4956" & c14_age == 4800 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-242781" & c14_age == 4890 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Bln-5540" & c14_age == 4892 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "I-15349" & c14_age == 4920 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-5383" & c14_age == 4940 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "DESCONOCIDO-17694" & c14_age == 4950 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "DESCONOCIDO-17693" & c14_age == 4960 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-5430" & c14_age == 4965 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-5833/40816" & c14_age == 5000 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-5532" & c14_age == 5010 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-5533" & c14_age == 5010 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrN-26226" & c14_age == 5045 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-5850" & c14_age == 5135 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-5851" & c14_age == 5135 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-5932" & c14_age == 5175 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-6009" & c14_age == 5210 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-156686" & c14_age == 5220 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-184182" & c14_age == 5230 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-4381/32340" & c14_age == 5245 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GifA-99113" & c14_age == 5330 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "unknown-6365" & c14_age == 5380 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UtC-7318" & c14_age == 5404 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Gif-11037" & c14_age == 5460 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrA-13624" & c14_age == 5480 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-166179" & c14_age == 5630 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-5832/40815" & c14_age == 5635 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Ua-16205" & c14_age == 5640 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-5860/41134" & c14_age == 5645 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-5834/40817" & c14_age == 5685 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-180981" & c14_age == 5690 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "I-17789" & c14_age == 5700 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-6828" & c14_age == 5710 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GAK-15223" & c14_age == 5710 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "KIK/KIA-5785/40878" & c14_age == 5715 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Wk-16148" & c14_age == 5831 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Ua-4821" & c14_age == 5960 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Q-(AM85B2b)" & c14_age == 5990 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "UGa-7565" & c14_age == 6120 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-74212" & c14_age == 6130 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "KIA-6790" & c14_age == 6144 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-7740" & c14_age == 6184 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-7741" & c14_age == 6184 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "AA-78527" & c14_age == 6203 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-166182" & c14_age == 6240 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "H-1754/1208" & c14_age == 6265 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-22340" & c14_age == 6270 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-7987" & c14_age == 6270 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "CSIC-114A" & c14_age == 6320 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "OxA-V-2392-26" & c14_age == 6341 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "OxA-2360" & c14_age == 6389 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-V-2360-22" & c14_age == 6389 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "OxA-V-2360-25" & c14_age == 6399 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-189076" & c14_age == 6400 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-V-26075" & c14_age == 6430 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-180721" & c14_age == 6440 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-8345" & c14_age == 6440 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-13157" & c14_age == 6590 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-233092" & c14_age == 6660 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "A-1" & c14_age == 6710 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GrN-26400" & c14_age == 6710 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-74220" & c14_age == 6730 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-8712" & c14_age == 6760 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-8714" & c14_age == 6760 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Gif-2368" & c14_age == 6780 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "unknown-8764" & c14_age == 6800 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "KIA-34796" & c14_age == 6810 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-141050" & c14_age == 6910 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "unknown-8936" & c14_age == 6970 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "TO-354" & c14_age == 6970 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "UGaMS-7196" & c14_age == 6990 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "H-2119" & c14_age == 7080 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "H-2119/1546" & c14_age == 7080 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "CSIC-637R" & c14_age == 7200 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-171910" & c14_age == 7280 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "TO-11819" & c14_age == 7300 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "TO-11819-R" & c14_age == 7300 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-7093" & c14_age == 7360 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Gif-3741" & c14_age == 7620 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-7495" & c14_age == 7710 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrA-22806" & c14_age == 8250 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GrN-18115" & c14_age == 9260 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-261" & c14_age == 9430 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "ICEN-92" & c14_age == 9530 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "GrN-17785" & c14_age == 9715 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "unknown-10061" & c14_age == 10190 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "Gif-95617" & c14_age == 10260 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "UGRA-147?" & c14_age == 12060 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "Gif-95630" & c14_age == 12240 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "unknown-10478" & c14_age == 12282 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "unknown-10601" & c14_age == 12896 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "GifA-96096" & c14_age == 13210 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "unknown-10662" & c14_age == 13300 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GifA-96139" & c14_age == 13320 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GrN-17255" & c14_age == 14020 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "unknown-10788" & c14_age == 14040 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "GifA-91130" & c14_age == 14250 & source_db == "pacea")) %>%
  dplyr::select(-duplicate)


# 4th pass: same lab no but different data or se, visual check
duplicate_testing_error <- vector("integer", nrow(pop_data_dr3))

d4_pop_data <- pop_data_dr3 %>% #Ensure data in the right order for testing for duplicates 
  dplyr::arrange(lab_nr, c14_age, c14_std, latitude, longitude, source_db)

for (i in 2:nrow(d4_pop_data)){#Generate vector to see whether duplicate rows based on coords
  duplicate_testing_error [[i]] <- if_else((d4_pop_data$lab_nr[[i]] == d4_pop_data$lab_nr[[(i-1)]])
                                           , 1, 0 )
}

pop_data_check_dr4 <- d4_pop_data %>% #insert duplicates, visually check this!
  dplyr::mutate(duplicate = duplicate_testing_error)

pop_data_dr4 <- pop_data_check_dr4 %>% #Manually checked all of these duplicates
  dplyr::filter(!(lab_nr == "AA-29648" & c14_age == 3565 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "AA-57439" & c14_age == 4604 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "AA-59519" & c14_age == 7256 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "AA-8647" & c14_age == 9988 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Beta-123555" & c14_age == 3650 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-124523" & c14_age == 4460 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-124525" & c14_age == 4040 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-124540" & c14_age == 3500 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-125110" & c14_age == 7320 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-126686" & c14_age == 5460 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-127449" & c14_age == 7120 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-136676" & c14_age == 6900 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-136677" & c14_age == 7000 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-142289" & c14_age == 6510 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-162879" & c14_age == 17850 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Beta-164901" & c14_age == 4540 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-165793" & c14_age == 6350 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-166181" & c14_age == 5810 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-166228" & c14_age == 5020 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-166229" & c14_age == 4250 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-166231" & c14_age == 6010 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-169264" & c14_age == 4450 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-169546" & c14_age == 4430 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-176897" & c14_age == 4280 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-179900" & c14_age == 5980 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-180980" & c14_age == 3869 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-186855" & c14_age == 3850 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-188925" & c14_age == 3310 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-188926" & c14_age == 3360 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-189075" & c14_age == 5170 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-191083" & c14_age == 6850 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-193745" & c14_age == 4110 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-193760" & c14_age == 7000 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-197384" & c14_age == 6100 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-197385" & c14_age == 6380 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-202343" & c14_age == 5100 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-208132" & c14_age == 6120 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-208133" & c14_age == 6150 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-208134" & c14_age == 6320 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-220914" & c14_age == 6110 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-222342" & c14_age == 6590 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-222444" & c14_age == 4000 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-223091" & c14_age == 5880 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-225217" & c14_age == 4710 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-225218" & c14_age == 5080 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-225224" & c14_age == 5010 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-229791" & c14_age == 3920 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-231876" & c14_age == 5890 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-231880" & c14_age == 6600 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-232340" & c14_age == 6020 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-232341" & c14_age == 6800 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-232342" & c14_age == 6780 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-232484" & c14_age == 5880 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-235487" & c14_age == 3450 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-235488" & c14_age == 3960 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-236821" & c14_age == 3310 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-240897" & c14_age == 5010 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-247406" & c14_age == 5340 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-248523" & c14_age == 6020 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-250404" & c14_age == 5000 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-252263" & c14_age == 5120 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-256325" & c14_age == 3910 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Beta-258466" & c14_age == 3340 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-278255" & c14_age == 5262 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-278256" & c14_age == 5129 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-288933" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288934" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288935" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288936" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288937" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288938" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-288939" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-298124" & c14_age == 6290 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-298125" & c14_age == 6270 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-298126" & c14_age == 6180 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299101" & c14_age == 3430 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "Beta-299302" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299303" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299307" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299308" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299309" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299311" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-299312" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301219" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301220" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301221" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301222" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301223" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301224" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301225" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-301226" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-307795" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-307796" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-307797" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-315717" & c14_age == 3980 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-321414" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-321415" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-321416" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-321418" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-321419" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-321420" & c14_age == 3550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-330866" & c14_age == 10520 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Beta-58933" & c14_age == 8790 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-58934" & c14_age == 6189 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-61490" & c14_age == 5580 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Beta-61490" & c14_age == 5880 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-64940" & c14_age == 4100 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-72552" & c14_age == 5000 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-72553" & c14_age == 5100 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-74311" & c14_age == 6180 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-74313" & c14_age == 6130 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-75211" & c14_age == 3710 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Beta-75214" & c14_age == 5790 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-75214" & c14_age == 5970 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Beta-79492" & c14_age == 4790 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-80602" & c14_age == 5320 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-83126" & c14_age == 4029 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-89286" & c14_age == 6060 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Beta-90884" & c14_age == 5920 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Bln-5050" & c14_age == 312 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "BM-1450" & c14_age == 12282 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "BM-2275R" & c14_age == 6440 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "BM-2276r" & c14_age == 7840 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "BM-2347" & c14_age == 4027 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "BM-2363" & c14_age == 7890 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CAMS-77427" & c14_age == 4720 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CAMS-83631" & c14_age == 5440 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CAMS-88195" & c14_age == 6050 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CNA-38" & c14_age == 3427 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CNA-554" & c14_age == 6212 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-1264" & c14_age == 5412 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-140" & c14_age == 4050 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-153" & c14_age == 4330 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-153B" & c14_age == 7950 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-157" & c14_age == 3390 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "CSIC-165" & c14_age == 3020 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-1664" & c14_age == 4939 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-173" & c14_age == 9800 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "CSIC-1822" & c14_age == 3788 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-31" & c14_age == 4070 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-494" & c14_age == 3420 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-57" & c14_age == 5400 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-629" & c14_age == 3970 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-637" & c14_age == 7180 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "CSIC-775" & c14_age == 5680 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "CSIC-823" & c14_age == 4970 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "CSIC-901" & c14_age == 3140 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "CSIC-927" & c14_age == 3921 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GAK-11395" & c14_age == 5890 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GAK-15222" & c14_age == 7030 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GAK-8959" & c14_age == 3120 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GAK-8968" & c14_age == 5440 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "GAK-8975" & c14_age == 7130 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Gif-10028" & c14_age == 9830 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Gif-9708" & c14_age == 8620 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "GifA-95146" & c14_age == 11270 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GifA-95360" & c14_age == 11630 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GifA-95362" & c14_age == 13700 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GifA-95370" & c14_age == 14260 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GifA-98156" & c14_age == 14750 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "GifA-99112" & c14_age == 5480 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GifA-99112" & c14_age == 5580 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-18296" & c14_age == 6730 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrA-21550" & c14_age == 7350 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-22825" & c14_age == 5220 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-24791" & c14_age == 5880 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-27278" & c14_age == 7955 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrA-45763" & c14_age == 7035 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrA-9312" & c14_age == 4250 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrA-9789" & c14_age == 6220 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-10744" & c14_age == 4000 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GrN-10971" & c14_age == 4230 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GrN-10973" & c14_age == 3960 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GrN-14599" & c14_age == 5380 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GrN-15760" & c14_age == 3340 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-16073" & c14_age == 5270 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrN-17698" & c14_age == 5230 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-18667" & c14_age == 5170 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GrN-19048" & c14_age == 4820 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-19569" & c14_age == 6755 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrN-19569" & c14_age == 6775 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "GrN-21353" & c14_age == 4990 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "GrN-22442" & c14_age == 7920 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-26225" & c14_age == 4970 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "GrN-26398" & c14_age == 7340 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "GrN-5598" & c14_age == 3833 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "GrN-5992" & c14_age == 5285 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GrN-5993" & c14_age == 17970 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "GrN-6634" & c14_age == 3525 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GrN-7004" & c14_age == 3955 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "GX-25854" & c14_age == 11950 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "HAR-146" & c14_age == 4090 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "I-11364" & c14_age == 3220 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "I-16442" & c14_age == 4690 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "I-17763" & c14_age == 5660 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "I-9239" & c14_age == 9560 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "I-9867" & c14_age == 5715 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-103" & c14_age == 6930 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-103" & c14_age == 6930 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "ICEN-1049" & c14_age == 4310 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-107" & c14_age == 7750 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-114" & c14_age == 4340 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-1149" & c14_age == 4460 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-1227" & c14_age == 6970 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "ICEN-1242" & c14_age == 4100 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "ICEN-146" & c14_age == 6970 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "ICEN-146" & c14_age == 6970 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-150" & c14_age == 7010 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-151" & c14_age == 7350 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-153" & c14_age == 7960 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-211" & c14_age == 7590 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-211" & c14_age == 7970 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-213" & c14_age == 7520 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-215" & c14_age == 7500 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-266" & c14_age == 10020 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "ICEN-267" & c14_age == 10090 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "ICEN-273" & c14_age == 6730 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "ICEN-277" & c14_age == 6760 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-278" & c14_age == 6720 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-303" & c14_age == 5210 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "ICEN-405" & c14_age == 2920 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "ICEN-416" & c14_age == 6940 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-417" & c14_age == 6980 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-448" & c14_age == 4140 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-60" & c14_age == 4200 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-64" & c14_age == 4490 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-645" & c14_age == 6420 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-677" & c14_age == 5960 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "ICEN-718" & c14_age == 7210 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-720" & c14_age == 7530 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-729" & c14_age == 7140 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-740" & c14_age == 4460 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-748" & c14_age == 10400 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "ICEN-76" & c14_age == 6060 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "ICEN-77" & c14_age == 7850 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-80" & c14_age == 9590 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-81" & c14_age == 9410 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-826" & c14_age == 3870 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "ICEN-85" & c14_age == 3330 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "ICEN-873" & c14_age == 6540 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-899" & c14_age == 7110 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-935" & c14_age == 5860 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "ICEN-939" & c14_age == 4630 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "ICEN-979" & c14_age == 3720 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "KIA-16581" & c14_age == 7075 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "KIA-30181" & c14_age == 6610 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "KIA-36268" & c14_age == 2795 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "KIA-7258" & c14_age == 3891 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "KIA-7260" & c14_age == 4234 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "KIK-56" & c14_age == 3790 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Ly-1198" & c14_age == 7750 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ly-2850" & c14_age == 6550 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ly-4420" & c14_age == 5929 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "MC-1112" & c14_age == 4700 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "MC-1112/MC1112b" & c14_age == 4650 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "MC-1470" & c14_age == 5100 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "MC-1470" & c14_age == 5180 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "MC-1474" & c14_age == 5200 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "MC-1477" & c14_age == 5470 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "MC-2140" & c14_age == 11050 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "MC-2297" & c14_age == 3990 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "OxA-10050" & c14_age == 10088 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-10191" & c14_age == 6275 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-10994" & c14_age == 3895 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "OxA-1131" & c14_age == 7010 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-17731" & c14_age == 4500 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "OxA-2360-25" & c14_age == 6399 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-26069" & c14_age == 6458 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "OxA-26075" & c14_age == 6430 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-4963" & c14_age == 3775 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-5507" & c14_age == 4410 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-5508" & c14_age == 4410 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-5514" & c14_age == 4370 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-5525" & c14_age == 7720 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-5534" & c14_age == 4010 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-5535" & c14_age == 4605 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-6580" & c14_age == 5840 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-7150" & c14_age == 6870 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-7300" & c14_age == 7516 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-7433" & c14_age == 4590 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "OxA-8571" & c14_age == 8686 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "OxA-8659" & c14_age == 10880 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "OxA-8882" & c14_age == 6140 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "OxA-9017" & c14_age == 15440 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "OxA-972" & c14_age == 12930 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Poz-19929" & c14_age == 3320 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Pta-9163" & c14_age == 6260 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Q-2493" & c14_age == 6660 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Q-2495" & c14_age == 6470 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Q-2496" & c14_age == 6050 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Q-2497" & c14_age == 6350 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Q-2499" & c14_age == 6730 & source_db == "mclaughlin")) %>%
  dplyr::filter(!(lab_nr == "RCD-2110" & c14_age == 3520 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Sac-1458" & c14_age == 3020 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Sac-1459" & c14_age == 6560 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Sac-1479" & c14_age == 5420 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Sac-1558" & c14_age == 6740 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Sac-1558" & c14_age == 6760 & source_db == "mclaughlin")) %>%
  dplyr::filter(!(lab_nr == "Sac-1560" & c14_age == 7200 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Sac-1587" & c14_age == 8620 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Sac-1594" & c14_age == 6070 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Sac-1608" & c14_age == 6200 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Sac-1794" & c14_age == 5690 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Sac-1893" & c14_age == 5420 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Sac-2069" & c14_age == 3930 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Sac-2102" & c14_age == 6250 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Sac-2168" & c14_age == 3640 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "TO-11219" & c14_age == 5980 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "TO-11863" & c14_age == 6650 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "TO-11864" & c14_age == 6890 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "TO-354" & c14_age == 6960 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "TO-358" & c14_age == 4170 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "TO-359a" & c14_age == 6960 & source_db == "mclaughlin")) %>%
  dplyr::filter(!(lab_nr == "TO-359a" & c14_age == 6970 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Ua-10272" & c14_age == 3425 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-10375" & c14_age == 5235 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Ua-12467" & c14_age == 5785 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-12467" & c14_age == 5885 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-19444" & c14_age == 6190 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-19499" & c14_age == 3065 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-20005" & c14_age == 5055 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Ua-20011" & c14_age == 5190 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-24423" & c14_age == 5495 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-24426" & c14_age == 5230 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-24651" & c14_age == 10360 & source_db == "14cpalaeolithic")) %>%
  dplyr::filter(!(lab_nr == "Ua-26019" & c14_age == 4070 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "Ua-3126" & c14_age == 4856 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-3128" & c14_age == 4955 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-34552" & c14_age == 5135 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-34708" & c14_age == 4080 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-35662" & c14_age == 4125 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-35665" & c14_age == 3830 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Ua-36023" & c14_age == 3660 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Ua-36206" & c14_age == 5625 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-36209" & c14_age == 6085 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-37834" & c14_age == 6090 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Ua-3903" & c14_age == 3530 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-39481" & c14_age == 3511 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-39867" & c14_age == 4415 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Ua-4301" & c14_age == 15615 & source_db == "pacea")) %>%
  dplyr::filter(!(lab_nr == "Ua-4821" & c14_age == 6010 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UBAR-100" & c14_age == 5100 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UBAR-209" & c14_age == 4860 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "UBAR-258" & c14_age == 2380 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UBAR-274" & c14_age == 5280 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "UBAR-671" & c14_age == 3305 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UBAR-673" & c14_age == 3360 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "UBAR-698" & c14_age == 3590 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UBAR-760" & c14_age == 6405 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "UBAR-760" & c14_age == 6450 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "UBAR-865" & c14_age == 3040 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-103" & c14_age == 3470 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-107" & c14_age == 3190 & source_db == "eubar")) %>%
  dplyr::filter(!(lab_nr == "UGRA-109" & c14_age == 3440 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UGRA-155" & c14_age == 3360 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UGRA-155" & c14_age == 3450 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-164" & c14_age == 3921 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UGRA-185" & c14_age == 3930 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-19" & c14_age == 3260 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-254" & c14_age == 6160 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "UGRA-260" & c14_age == 3530 & source_db == "radonb")) %>%
  dplyr::filter(!(lab_nr == "UGRA-306" & c14_age == 3480 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-309" & c14_age == 3990 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-326" & c14_age == 3050 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "UGRA-327" & c14_age == 5160 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Wk-13682" & c14_age == 6185 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Wk-13692" & c14_age == 6712 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-14793" & c14_age == 6737 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-14794" & c14_age == 6821 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-14797" & c14_age == 6860 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Wk-14797" & c14_age == 6862 & source_db == "mclaughlin")) %>%
  dplyr::filter(!(lab_nr == "Wk-27995" & c14_age == 4739 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Wk-28006" & c14_age == 4775 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-28008" & c14_age == 4932 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-28010" & c14_age == 4765 & source_db == "pardo")) %>%
  dplyr::filter(!(lab_nr == "Wk-28635" & c14_age == 5441 & source_db == "balsera")) %>%
  dplyr::filter(!(lab_nr == "Wk-8939" & c14_age == 8580 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Wk-8950" & c14_age == 8640 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Wk-8951" & c14_age == 8400 & source_db == "radon")) %>%
  dplyr::filter(!(lab_nr == "Wk-9213" & c14_age == 7370 & source_db == "drake")) %>%
  dplyr::filter(!(lab_nr == "Wk-9744" & c14_age == 5753 & source_db == "drake")) %>%
  dplyr::select(-duplicate) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr, latitude, longitude, source_db)


gen_site_id <-  dplyr::group_indices(pop_data_dr4, latitude2, longitude2) #necessary for later binning of sites in SPD, used rounded coords

pop_data_dr4<- pop_data_dr4%>% 
  dplyr::mutate(site_id = gen_site_id) %>% 
  dplyr::filter(c14_std>0)

cleaned_data_all <- pop_data_dr4 %>% #Filter to approx dates required
  dplyr::filter(between(c14_age, 1000, 15000)) %>% 
  dplyr::mutate(calibration_curve = dplyr::if_else(material %in% c("Littorina l.",
                                                                   "Malacofauna",
                                                                   "mollusc",
                                                                   "P. vulgata shell B1 (Middle Mahgdalenian)",
                                                                   "shel from level 1Cm",
                                                                   "shell",
                                                                   "Shell",
                                                                   "Shell . Obusata",
                                                                   "shell from level 1 B M",
                                                                   "shell from level 1B P",
                                                                   "Shell, Littorina littorea",
                                                                   "shell, Patella vulgata"),
                                                   "marine20", "intcal20")) %>% 
  dplyr::mutate(calibration_curve = dplyr::if_else(is.na(material), "intcal20", calibration_curve))

shell_data <- cleaned_data_all %>% #Lots of manual searches here for offset values
  dplyr::filter(material %in% c("Littorina l.",
                                "Malacofauna",
                                "mollusc",
                                "P. vulgata shell B1 (Middle Mahgdalenian)",
                                "shel from level 1Cm",
                                "shell",
                                "Shell",
                                "Shell . Obusata",
                                "shell from level 1 B M",
                                "shell from level 1B P",
                                "Shell, Littorina littorea",
                                "shell, Patella vulgata")) %>%  
  dplyr::left_join(dplyr::select(mclaughlin_data_4258, lab_nr, Reservior, ReserviorError), by = "lab_nr") %>% 
  dplyr::rename(reservoir_offset = Reservior, reservoir_std = ReserviorError) %>% 
  dplyr::select(!geometry)  %>% 
  dplyr::mutate(offset_source = dplyr::if_else(!is.na(reservoir_offset), 
                                               "McLaughlin, T.R., Gómez-Puche, M., Cascalheira, J., Bicho, N., Fernández-López de Pablo, J., 2021. Late Glacial and Early Holocene human demographic responses to climatic and environmental change in Atlantic Iberia. Philos. Trans. R. Soc. B Biol. Sci. 376, 20190724. https://doi.org/10.1098/rstb.2019.0724",
                                               "")) %>%  
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Santa Sofia", 95, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Santa Sofia", 15, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "Santa Sofia", "Soares, A. M. M.; Arruda, A. M., 2017.Cronologia de radiocarbono para a Idade do Ferro Orientalizante no território português. Uma leitura crítica dos dados arqueométricos e arqueológicos. In Barcleó, J. A.; Bogdanovic, I.; Morell, B., eds. - IberCrono 2016 Cronometrías Para la Historia de la Península Ibérica (Chronometry for the History of the Iberian Peninsula). Actas del Congrso de Cronometrias para la Historia de la Península Ibérica, Barcelona, 17-19 de octubre 2016. Barcelona: CEUR. pp. 235-259. http://hdl.handle.net/10451/30122",
                                               offset_source)) %>%   #Mannual search
  dplyr::mutate(reservoir_offset = dplyr::if_else(site %in% c("Alcalar (settlement)", "Alcalar 7", "Alcalar Monument 7") , 69, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site %in% c("Alcalar (settlement)", "Alcalar 7", "Alcalar Monument 7"), 17, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site %in% c("Alcalar (settlement)", "Alcalar 7", "Alcalar Monument 7"), "Morán H., María E. 2015. El asentamiento prehistórico de Alcalar (Portimão, Portugal): la organización del territorio y el proceso de formación de un estado prístino en el tercer milenio A.N.E. (Tesis Doctoral Inédita). Universidad de Sevilla. https://hdl.handle.net/11441/73261",
                                               offset_source)) %>% #Mannual search
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Castro de Chibanes" , 95, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Castro de Chibanes", 15, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "Castro de Chibanes", "Silva, C.T. and Soares, J., 2014. O Castro de Chibanes (Palmela) e o tempo social do III milénio BC na Estremadura. Setúbal Arqueológica-II Encontro de Arqueologia da Arrábida. Homenagem a AI Marques da Costa, 15, pp.105-172.http://hdl.handle.net/10451/10914",
                                               offset_source)) %>% #Mannual search
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Magoito" & is.na(reservoir_offset) , 95, reservoir_offset)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Magoito" & is.na(reservoir_std), 15, reservoir_std)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(offset_source = dplyr::if_else(site == "Magoito" & offset_source == "", "same site as within McLaughlin et al. 2021",
                                               offset_source)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "La Loma" , 22, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "La Loma", 35, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "La Loma", "Aranda Jiménez, G., Cámalich Massieu, M.D., Martín Socas, D., Morgado, A., Martínez-Sevilla, F., Lozano, J., Rodríguez Rodríguez, A., Mancilla Cabello, M.I. and Román Punzón, J., 2012. La Loma (Íllora, Granada). Un yacimiento de fosas del VI-IV milenios cal. BC. Consejería de Cultura Junta de Andalucía, Sevilla. ISBN: 978-84-9959-105-6",
                                               offset_source)) %>% #Mannual search  
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Pico Ramos" , -1, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Pico Ramos", 42, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "Pico Ramos", "Zapata, L., 2017. Level 4 of the cave of Pico Ramos (Muskiz, Bizkaia): excavation, stratigraphy, chronology and materials. In: Zapata, L., (ed.) The shell midden of Pico Ramos (Muskiz, Bizkaia): Humans on the Basque coast during the 6th and 5th millennium B.C. Trrres , Bilbao , pp. 1-514. ISBN: 978-84-617-8618-3",
                                               offset_source)) %>% #Mannual search 
  dplyr::mutate(reservoir_offset = dplyr::if_else(lab_nr %in% c("QAM-85B2b", "Wk-28050"),  140, reservoir_offset)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_std = dplyr::if_else(lab_nr %in% c("QAM-85B2b", "Wk-28050"), 40, reservoir_std)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(offset_source = dplyr::if_else(lab_nr %in% c("QAM-85B2b", "Wk-28050"), "same site as within McLaughlin et al. 2021",
                                               offset_source)) %>% #Same as McLaughlin et al. 2021 
  dplyr::mutate(reservoir_offset = dplyr::if_else(lab_nr %in% c("Wk-17029", "Wk-6075"),  -116, reservoir_offset)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_std = dplyr::if_else(lab_nr %in% c("Wk-17029", "Wk-6075"), 40, reservoir_std)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(offset_source = dplyr::if_else(lab_nr %in% c("Wk-17029", "Wk-6075"), "same site as within McLaughlin et al. 2021",
                                               offset_source)) %>% #Same as McLaughlin et al. 2021 
  dplyr::mutate(reservoir_offset = dplyr::if_else(lab_nr == "ICEN-873" , 69, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(lab_nr == "ICEN-873" , 17, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(lab_nr == "ICEN-873" , "Morán H., María E. 2015. El asentamiento prehistórico de Alcalar (Portimão, Portugal): la organización del territorio y el proceso de formación de un estado prístino en el tercer milenio A.N.E. (Tesis Doctoral Inédita). Universidad de Sevilla. https://hdl.handle.net/11441/73261",
                                               offset_source)) %>% #Mannual search 
  dplyr::mutate(reservoir_offset = dplyr::if_else(lab_nr == "ICEN-645" , 69, reservoir_offset)) %>% #same site as Morán H. and María E. 2015 
  dplyr::mutate(reservoir_std = dplyr::if_else(lab_nr == "ICEN-645" , 17, reservoir_std)) %>% #same site as Morán H. and María E. 2015 
  dplyr::mutate(offset_source = dplyr::if_else(lab_nr == "ICEN-645" , "#same site as within Morán H. and María E. 2015 ",
                                               offset_source)) %>% #same site as Morán H. and María E. 2015 
  dplyr::mutate(reservoir_offset = dplyr::if_else(lab_nr == "Beta-168461",  -110, reservoir_offset)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_std = dplyr::if_else(lab_nr == "Beta-168461", 40, reservoir_std)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(offset_source = dplyr::if_else(lab_nr == "Beta-168461", "same site as within McLaughlin et al. 2021",
                                               offset_source)) %>% #Same as McLaughlin et al. 2021
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Cierro" , -117, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Cierro", 70, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "Cierro", "Álvarez-Fernández, E., Bécares, J., Pardo, J.F.J., Agirre-Uribesalgo, A., Álvarez-Alonso, D., Aparicio, M.T., Barrera-Mellado, I., Carral, P., Carriol, R.P., Cubas, M. and Cueto, M., 2020. Palaeoenvironmental and chronological context of human occupations at El Cierro cave (Northern Spain) during the transition from the late Upper Pleistocene to the early Holocene. Journal of Archaeological Science: Reports, 29, p.102138. https://doi.org/10.1016/j.jasrep.2019.102138",
                                               offset_source)) %>% #Mannual search 
  dplyr::mutate(reservoir_offset = dplyr::if_else(site == "Cova Rosa" , -117, reservoir_offset)) %>% #Mannual search
  dplyr::mutate(reservoir_std = dplyr::if_else(site == "Cova Rosa", 70, reservoir_std)) %>% #Mannual search
  dplyr::mutate(offset_source = dplyr::if_else(site == "Cova Rosa", "Álvarez-Fernández, E., Jordá Pardo, J. F., Arias, P., Bécares, J., Martín-Jarque, S., Portero, R., Teira, L. and Douka, K. (2021) “RADIOCARBON DATES FOR THE LATE PLEISTOCENE AND EARLY HOLOCENE OCCUPATIONS OF COVA ROSA (RIBADESELLA, ASTURIAS, SPAIN),” Radiocarbon. Cambridge University Press, 63(3), pp. 1053–1072. doi: 10.1017/RDC.2021.18. https://doi.org/10.1017/RDC.2021.18",
                                               offset_source)) %>%  #Mannual search 
  dplyr::mutate(reservoir_offset = dplyr::if_else(is.na(reservoir_offset), -9999, reservoir_offset))
  
  
  
rio::export(shell_data, "other_output/marine.xlsx") #Reference information for marine offset information

shell_data_sf <- shell_data %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>% 
  dplyr::filter(reservoir_offset > -9999) 

radio_citation <- dplyr::tribble(
  ~dataset, ~dataset_citation,
  "14cpalaeolithic", "Vermeersch, P.M., 2020. Radiocarbon Palaeolithic Europe database: A regularly updated dataset of the radiometric data regarding the Palaeolithic of Europe, Siberia included. Data Br. 31, 105793. https://doi.org/10.1016/j.dib.2020.105793",
  "balsera", "Balsera, V., Díaz-del-Río, P., Gilman, A., Uriarte, A., Vicent, J.M., 2015. Approaching the demography of late prehistoric Iberia through summed calibrated date probability distributions (7000–2000 cal BC). Quat. Int. 386, 208–211. https://doi.org/10.1016/j.quaint.2015.06.022",
  "drake", "Drake, B.L., Blanco-González, A., Lillios, K.T., 2017. Regional Demographic Dynamics in the Neolithic Transition in Iberia: Results from Summed Calibrated Date Analysis. J. Archaeol. Method Theory 24, 796–812. https://doi.org/10.1007/s10816-016-9286-y",
  "eubar", "Capuzzo, G., Boaretto, E., Barceló, J.A., 2014. EUBAR: A Database of 14 C Measurements for the European Bronze Age. A Bayesian Analysis of 14 C-Dated Archaeological Contexts from Northern Italy and Southern France. Radiocarbon 56, 851–869. https://doi.org/10.2458/56.17453",
  "euroevol", "Manning, K., Colledge, S., Crema, E.R., Shennan, S., Timpson, A., 2016. The Cultural Evolution of Neolithic Europe. EUROEVOL Dataset 1: Sites, Phases and Radiocarbon Data. J. Open Archaeol. Data 5. https://doi.org/10.5334/joad.40",
  "mclaughlin", "McLaughlin, T.R., Gómez-Puche, M., Cascalheira, J., Bicho, N., Fernández-López de Pablo, J., 2021. Late Glacial and Early Holocene human demographic responses to climatic and environmental change in Atlantic Iberia. Philos. Trans. R. Soc. B Biol. Sci. 376, 20190724. https://doi.org/10.1098/rstb.2019.0724",
  "pacea", "d’Errico, F., Banks, W.E., Vanhaeren, M., Laroulandie, V., Langlais, M., 2011. PACEA geo-referenced radiocarbon database. PaleoAnthropology 2011, 1–12. https://doi.org/doi:10.4207/PA.2011.ART40",
  "pardo", "Pardo-Gordó, S., García Puchol, O., Bernabeu Aubán, J., Diez Castillo, A., 2019. Timing the Mesolithic-Neolithic Transition in the Iberian Peninsula: The Radiocarbon Dataset. J. Open Archaeol. Data 7. https://doi.org/10.5334/joad.49",
  "radon", "Hinz, M., Furholt, M., Müller, J., Rinne, C., Raetzel-Fabian, D., Sjögren, K.-G., Wotzka, H.-P., 2012. RADON - Radiocarbon dates online 2012. Central European database of 14C dates for the Neolithic and the Early Bronze Age. J. Neolit. Archaeol. 0. https://doi.org/10.12766/jna.2012.65",
  "radonb", "Kniesel, J., Hinz, M., Rinne, C., 2014. Radon-B. In: http://radon-b.ufg.uni-kiel.de."
)

cleaned_data <- cleaned_data_all %>% 
  dplyr::left_join(dplyr::select(shell_data, db_code, reservoir_offset, reservoir_std), by = "db_code") %>% 
  dplyr::mutate(reservoir_offset = ifelse(is.na(reservoir_offset), 0, reservoir_offset)) %>% 
  dplyr::mutate(reservoir_std = ifelse(is.na(reservoir_std), 0, reservoir_std)) %>%  
  dplyr::filter(reservoir_offset > -9999) %>% 
  dplyr::mutate(reservoir_offset = as.numeric(reservoir_offset)) %>% 
  dplyr::mutate(reservoir_std = as.numeric(reservoir_std)) %>% 
  dplyr::mutate(period = dplyr::if_else(period %in% c("Neolítico", #This is a bit back to front, but for the purposes of generating a file in pub, best way to go
                                                      "Ferro I A",
                                                      "Ferro I", 
                                                      "Ferro I B",
                                                      "Bronze final A",
                                                      "Bronce",
                                                      "Antic",          
                                                      "Adlerberg Gruppe",
                                                      "Calcolítico",
                                                      "Spätchalkolithikum",
                                                      "Bošáca Group",
                                                      "Chalkolithikum spät",
                                                      "Frühchalkolithikum",
                                                      "Neolithic",        
                                                      "UN"), "Neolithic", NA_character_)) %>% 
  dplyr::arrange(c14_age, c14_std, lab_nr) %>% 
  dplyr::rename(neolithic_tag = period) %>% 
  dplyr::rename(dataset = source_db) %>%  
  dplyr::left_join(radio_citation, by = "dataset")



radiocarbon_data <- cleaned_data %>% 
  dplyr::select(c14_age,
                c14_std,
                longitude,
                latitude,
                lab_nr,
                site,
                site_id,
                material,
                neolithic_tag,
                elevation,
                calibration_curve,
                reservoir_offset,
                reservoir_std,
                dataset,
                dataset_citation)

cleaned_data_sites <- length(unique(cleaned_data$site_id))

rio::export(radiocarbon_data, "other_output/pop/radiocarbon_data.csv") #This is the clean data used for analysis
  
  
  



