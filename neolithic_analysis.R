#Sweeney, Harrison and Vander Linden 2021. 
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


# Neolithic analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script uses the data generated within the
# charcoal_analyis.R and radiocarbon_analysis.R scrips to 
# assess the potential influence of the spread of farming 
# on fire history during the 10-3.5k period within mainland
# Iberia. 

# As with charcoal_analysis.R, the construction of charcoal
# composite curves are based on the code  published by Patrick 
# J. Bartlein. see "Analysis of the Global Charcoal Database
# – GCDv3". https://pjbartlein.github.io/GCDv3Analysis/
# The code  here allows for reproduction of the results
# generated within the corresponding paper



# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Extract the Neolithic date for each charcoal entity
# 3. Pre-bin Neolithic amended data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 4. Neolithic amended composite curve (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 5. Permutation work based on variable Neolithic dates
# 6. ARIMA intervention analysis




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




# 2. Extract the Neolithic date for each charcoal entity
# ---------------------------------------------------------
# ---------------------------------------------------------

# Get coordinate info for charcoal record locations
entity_list_neo <- rio::import("other_output/rpd/rpd_data_input.csv") %>% 
  dplyr::group_by(entity_name) %>% 
  dplyr::summarise(longitude = mean(longitude), latitude = mean(latitude), ID_ENTITY = mean(ID_ENTITY))
entity_lon_lat_neo <-   sf::st_as_sf(entity_list_neo, coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>%
  sf::st_transform(crs = 25830) %>%
  sf::st_coordinates()

#Extract record Neolithic date based on kriging interpolated surface
krige_neo_raster <- raster::raster("other_output/pop/kriging_neoloithic.grd") #Import surface

entity_neo_check <- raster::extract(krige_neo_raster, entity_lon_lat_neo) %>% #Extract value from surface by record location
  dplyr::as_tibble() %>% 
  dplyr::mutate(entity_name = entity_list_neo$entity_name)

if (sum(is.na(entity_neo_check$value)) > 0  ){
  NA_entity_neo <- entity_neo_check %>% #Identify NAs where marriage of coordinates not perfect with surface
    dplyr::filter(is.na(value)) %>% 
    dplyr::left_join(dplyr::select(entity_list_neo, entity_name, latitude, longitude), by = "entity_name")
  
  entity_neo_lon_lat_NA <- NA_entity_neo%>% #Get lon lat of NAs
    sf::st_as_sf(., coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>% 
    sf::st_transform(crs = 25830) %>% 
    sf::st_coordinates()
  
  entity_neo_NA <- raster::extract(krige_neo_raster, entity_neo_lon_lat_NA, buffer = 100000, fun = mean, df = T) %>% #Use buffer to get approximate values relating to NA   
    dplyr::as_tibble() %>% 
    dplyr::rename(value = 2) %>% 
    dplyr::mutate(entity_name = NA_entity_neo$entity_name) %>% 
    dplyr::select(value, entity_name) 
  
  entity_neo <- entity_neo_check %>% #Merge na and non-na values
    dplyr::filter(!(is.na(value))) %>% 
    rbind(entity_neo_NA) %>% 
    dplyr::rename(first_neo = value) %>% 
    dplyr::left_join(dplyr::select(entity_list_neo, entity_name, ID_ENTITY), by = "entity_name")
  
} else {
  entity_neo <- entity_neo_check %>% 
    dplyr::rename(first_neo = value) %>% 
    dplyr::left_join(dplyr::select(entity_list_neo, entity_name, ID_ENTITY), by = "entity_name")
}


# ---------------------------------------------------------





# 3. Pre-bin Neolithic amended data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------


datapath <- "./other_output/rpd/"
neo_sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv") #Load-in influx data
neo_window <- 800 # Set range of data to be considered

# Set baseperiod - here you will use this to select sites based on the degree of overlap with the baseperiod
basebeg <- 200
baseend <- 8200
basename <- paste0("zt", basebeg, "_", baseend, "k")

neo_dat <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>%
  dplyr::left_join(entity_neo, by = "ID_ENTITY") %>% #Insert extracted dates of first agriculture
  dplyr::mutate(EST_AGE = (EST_AGE - first_neo)) %>% #Re-centre dates based on first agriculture per site location
  dplyr::select(-c(entity_name, first_neo))

neo_min_est_age <- -neo_window # To set ages later
neo_max_est_age <- neo_window # To set ages later

length(unique(neo_dat$ID_ENTITY))  ##all entities that were transformed  in long format
neo_ent_vec2 <- unique(neo_dat$ID_ENTITY)

neo_lsover <- neo_ent_vec2 #This line of code enables subsetting of dataset if required (not utilised here)


# Binning Z scores
neo_Zbinning <- function(dt, bins, binhw){
  # Reshape data
  dat_res <- dt[,c("ID_ENTITY", "zt")] %>% 
    dplyr::group_by(ID_ENTITY) %>% 
    dplyr::mutate(number = row_number()) %>% 
    tidyr::spread(ID_ENTITY, zt) 
  dat_res <- dat_res[,-1] 
  dat_res <- as.matrix(dat_res)
  
  Age_res <- dt[,c("ID_ENTITY", "EST_AGE")] %>% 
    dplyr::group_by(ID_ENTITY) %>% 
    dplyr::mutate(number = row_number()) %>% 
    tidyr::spread(ID_ENTITY, EST_AGE)
  Age_res <- Age_res[,-1]
  Age_res <- as.matrix(Age_res)
  
  # Binning
  result <- matrix(ncol = length(Age_res[1,]), nrow = length(bins))
  colnames(result) <- colnames(dat_res) 
  for (k in 1:length(dat_res[1,])){
    if(length(dat_res[is.na(dat_res[,k]) == F, k]) != 0){
      for (i in 1:length(bins)){
        t <- na.omit(cbind(as.numeric(Age_res[,k]), as.numeric(dat_res[,k])))
        result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] < bins[i] + binhw, 2])
      }
    }
  }
  neo_BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
  return(neo_BinnedData)
}

neo_bw <- 100 #Standard for composite

neo_Zbin_dat <- neo_Zbinning(dt = neo_dat, bins = seq(neo_min_est_age, neo_max_est_age, neo_bw), binhw = neo_bw/2) 


# Give ages of each bin as a column
neo_bin_age <- seq(neo_min_est_age, neo_max_est_age,neo_bw)
neo_Zbin_dat <- cbind(neo_bin_age, neo_Zbin_dat)

# Limit this to the sites in lsover
class(neo_lsover)
neo_lsoverch <- as.character(neo_lsover)
class(neo_lsoverch)
class(colnames(neo_Zbin_dat))
neo_Zbin_dat <- neo_Zbin_dat[,c('neo_bin_age',neo_lsoverch)]
check <- colnames(neo_Zbin_dat)
check[!(check %in% neo_lsoverch)]
neo_lsoverch

# Dataframe for curve generation
neo_Zbin_dat <- as.data.frame(neo_Zbin_dat) 

# ---------------------------------------------------------




# 4. Neolithic amended composite curve (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

## Set up

# Set names
neo_basename <- paste0("zt", basebeg, "_", baseend, "k")
neo_binname <- paste0("bw", neo_bw)

# Locfit (half) window-width parameter
neo_hw <- 500 # bandwidth (smoothing parameter)

# Number of bootstrap samples/replications
nreps <- 1000

# Target ages for fitted values
neo_targbeg <- neo_min_est_age
neo_targend <- neo_max_est_age 
neo_targstep <- neo_bw 

# Array sizes
maxrecs <- 2000
maxreps <- 1000

# Plot output 
plotout <- "pdf" 



## Neo curve plotting

# Curve output path and file
neo_curvecsvpath <- "./figs/neo_rpd/"
neo_curvefile <- paste(neo_curvecsvpath,"neo_rpd",".csv", sep="") #adjust to particular data you use now. .
print(neo_curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { 
  neo_pdffile <- paste(neo_curvecsvpath,"neo_rpd", ".pdf", sep="")
  print(neo_pdffile)
}

# Read the list of sites to be processed
neo_ns <- as.numeric(colnames(neo_Zbin_dat))
neo_ns <- na.omit(neo_ns) 

# Arrays for data and fitted values
neo_age <- matrix(NA, ncol=length(neo_ns), nrow=maxrecs)
neo_influx <- matrix(NA, ncol=length(neo_ns), nrow=maxrecs)
neo_nsamples <- rep(0, maxrecs)
neo_targage <- seq(neo_targbeg,neo_targend,neo_targstep) 
neo_targage.df <- data.frame(neo_x=neo_targage)
neo_lowage <- neo_targage - neo_hw; neo_highage <- neo_targage + neo_hw
neo_ntarg <- length(neo_targage)
neo_yfit <- matrix(NA, nrow=length(neo_targage.df$neo_x), ncol=maxreps)

# Arrays for sample number and effective window span tracking
neo_ndec <- matrix(0, ncol=neo_ntarg, nrow=length(neo_ns))
neo_ndec_tot <- rep(0, neo_ntarg)
neo_xspan <- rep(0, neo_ntarg)
neo_ninwin <- matrix(0, ncol=neo_ntarg, nrow=length(neo_ns))
neo_ninwin_tot <- rep(0, neo_ntarg)


## Read data
# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0

for(i in 2:length(neo_Zbin_dat)){ 
  neo_indata <- neo_Zbin_dat[,c(i,1)] 
  neo_indata <- na.omit(neo_indata)  
  neo_nsamp <-  length(neo_indata$neo_bin_age) 
  if (neo_nsamp > 0) {
    ii <- ii+1
    neo_age[1:neo_nsamp,ii] <- neo_indata$neo_bin_age 
    neo_influx[1:neo_nsamp,ii] <- neo_indata[,1]  
    neo_nsamples[ii] <- neo_nsamp
  }
  print(i)
}

neo_nsites <- ii; neo_nsites #Number of sites with data

# Trim samples to age range
neo_influx[neo_age >= neo_targend+neo_hw] <- NA
neo_age[neo_age >= neo_targend+neo_hw] <- NA 

# Censor abs(influx) values > 10
neo_influx[abs(neo_influx) >= 10] <- NA
neo_age[abs(neo_influx) >= 10] <- NA 


##Find number of sites and samples contributing to fitted values

# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:neo_ntarg) {
  neo_agemax <- -1e32; neo_agemin <- 1e32
  for (j in 1:neo_nsites) {
    for (k in 1:neo_nsamples[j]) {
      if (!is.na(neo_age[k,j])) {
        ii <- (neo_age[k,j]-neo_targage[1])/neo_targstep + 1
        if (ii > 0 && ii <= neo_ntarg) {neo_ndec[j,ii] = 1}
        if (neo_age[k,j] >= neo_targage[i]-neo_hw && neo_age[k,j] <= neo_targage[i]+neo_hw) {
          neo_ninwin[j,i] = 1
          if (neo_agemax < neo_age[k,j]) {neo_agemax <- neo_age[k,j]}
          if (neo_agemin > neo_age[k,j]) {neo_agemin <- neo_age[k,j]}
        }
      }
    }
  }
  neo_ndec_tot[i] <- sum(neo_ndec[,i])
  neo_ninwin_tot[i] <- sum(neo_ninwin[,i])
  neo_xspan[i] <- neo_agemax - neo_agemin
}
proc.time() - ptm


## Curve-fitting and bootstrapping

# Composite curve
ptm <- proc.time()

# Reshape matrices into vectors 
neo_x <- as.vector(neo_age)
neo_y <- as.vector(neo_influx)
neo_lfdata <- data.frame(neo_x,neo_y)
neo_lfdata <- na.omit(neo_lfdata)
neo_x <- neo_lfdata$neo_x; neo_y <- neo_lfdata$neo_y

# Locfit
neo_loc01 <- locfit::locfit(neo_y ~ locfit::lp(neo_x, deg=1, h=neo_hw), maxk=800, family="qrgauss")
summary(neo_loc01)

# Get fitted values
neo_pred01 <- predict(neo_loc01, newdata=neo_targage.df, se.fit=TRUE)
neo_loc01_fit <- data.frame(neo_targage.df$neo_x, neo_pred01$fit)
neo_fitname <- paste("neo_locfit_",as.character(neo_hw), sep="")
colnames(neo_loc01_fit) <- c("neo_age", neo_fitname)
head(neo_loc01_fit)

#Detrend fitted values
fit_neo_loc01_fit <- lm(neo_locfit_500 ~ neo_age, data = neo_loc01_fit)
detrend_neo_loc01_fit <- dplyr::tibble(neo_age = neo_loc01_fit$neo_age, detrend_neo_locfit_500 = resid(fit_neo_loc01_fit)) %>% 
  as.data.frame()


proc.time() - ptm


#Bootstrap-by-site confidence intervals
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications
neo_curvecsvpath
neo_pdffile

if (plotout == "pdf") {pdf(file=neo_pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(neo_x, neo_y, xlab="Age before and after the start of the Neolithic", ylab="Detrended charcoal z-score composite", xlim=c(750,-750), ylim=c(-1,1), type="n", mgp = c(2.4, 1, 0), cex.lab = 1.2, xaxt = "n") ##make sure xlim is right.
axis(1, at = c(750, 500, 250, 0, -250, -500, -750), labels = c(-750, -500, -250, 0, 250, 500, 750))
grid()

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) { 
  print(i)
  neo_randsitenum <- sample(seq(1:neo_nsites), neo_nsites, replace=TRUE)
  neo_x <- as.vector(neo_age[,neo_randsitenum])
  neo_y <- as.vector(neo_influx[,neo_randsitenum])
  neo_lfdata <- data.frame(neo_x,neo_y)
  neo_lfdata <- na.omit(neo_lfdata)
  neo_x <- neo_lfdata$neo_x; neo_y <- neo_lfdata$neo_y
  neo_locboot <- locfit::locfit(neo_y ~ locfit::lp(neo_x, deg=1, h=neo_hw), maxk=800, maxit=20, family="qrgauss")
  neo_predboot <- predict(neo_locboot, newdata=neo_targage.df, se.fit=TRUE)
  neo_boot <- data.frame(neo_predboot$fit, neo_targage.df$neo_x) #Detrend
  fit_neo_boot <- lm(neo_predboot$fit ~ neo_targage.df.neo_x, data = neo_boot) #Detrend
  neo_yfit[,i] <- resid(fit_neo_boot)
  lines(neo_targage.df$neo_x, neo_yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}
warnings()

# Step 3 -- Plot the unresampled (initial) fit
neo_fitname <- paste("neo_locfit_",as.character(neo_hw), sep="")
colnames(neo_loc01_fit) <- c("neo_age", neo_fitname)
lines(detrend_neo_loc01_fit[,1], detrend_neo_loc01_fit[,2], lwd=2, col="red") 

# Step 4 -- Find and add bootstrap CIs
neo_yfit95 <- apply(neo_yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
neo_yfit05 <- apply(neo_yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(neo_targage.df$neo_x, neo_yfit95, lwd=1, col="red")
lines(neo_targage.df$neo_x, neo_yfit05, lwd=1, col="red")

abline(v = 0, col = "blue", lwd = 1, lty = 2)
text.default(x = -50, y = 0.7, "Neolithic start date", col = "blue", pos = 4, cex = 1)
text.default(x = 1500, y = -1, "Before", col = "blue", pos = 4, cex = 1)
text.default(x = -1000, y = -1, "After", col = "blue", pos = 4, cex = 1)

if (plotout == "pdf") {dev.off()}


# ---------------------------------------------------------








# 5. Permutation work based on variable Neolithic dates
# ---------------------------------------------------------
# ---------------------------------------------------------
# This part of the script regenerates the Neolithic surface
# x times (depending on nreps) based on amendments to the
# calibrated dates of Neolithic tagged data. As such it is
# slow.

# Set up
nreps_neo <- 1000 #how many times a date is selected
p_vario_error <- c(0) #To help catch when an error occurs

p_all_neo_curveout <- data.frame(matrix(, nrow=17, ncol=0)) #Dataframe to store output
p_all_neo_curveout_loc <- data.frame(matrix(, nrow=17, ncol=0)) #Dataframe to store output

radiocarbon_data_cal <- rio::import("./other_output/pop/radiocarbon_data_cal.csv") %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>% #Convert to sf
  sf::st_transform(crs = 25830) #Flat projection
radiocarbon_data_xy <- sf::st_coordinates(radiocarbon_data_cal) %>% #Extract coordinate information
  dplyr::as_tibble()

# Generate map window
iberia_map <- maps::map("world", c("Spain(?!:)", "Portugal", "Andorra"), fill = T, resolution = 0) #Download map of Iberia
iberia_map <- sf::st_as_sf(iberia_map) #Convert to sf
iberia_map_25830 <- sf::st_transform(iberia_map, crs = 25830) #Change the crs to flat projection
iberia_owin_25830 <- spatstat.geom::as.owin(iberia_map_25830) #Generate owin window
iberia_owin_25830_exp <- spatstat.core::expand.owin(iberia_owin_25830, distance = 2000) #Expand owin to try and ensure that samples are not excluded due to coord inaccuracies 
iberia_map_exp_25830 <-  sf::st_as_sf(iberia_owin_25830_exp) #Convert back to sf
sf::st_crs(iberia_map_exp_25830) <- 25830 #Set crs
iberia_map_exp <- sf::st_transform(iberia_map_exp_25830, crs = 4258) #Change crs back to ellipsoid for raster




#First run permutation with calibrated dates of Neolithic tagged data
for (q in 1:nreps_neo){
  
  print(q)  
  p_neo_adjust <- round(runif(length(radiocarbon_data_xy$X), min = -200, max = 200))
  
  p_radiocarbon_data_cal <- radiocarbon_data_cal %>% 
    dplyr::mutate(median_caldates = median_caldates - p_neo_adjust)
  
  
  p_neo_data <- p_radiocarbon_data_cal %>% #Remove samples with non-neolithic tags
    dplyr::mutate(x = radiocarbon_data_xy$X, y = radiocarbon_data_xy$Y) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(date_id = radiocarbon_caldates$date_id) %>% 
    dplyr::select(x, y, median_caldates, c14_age: date_id) %>% 
    dplyr::filter(neolithic_tag ==  "Neolithic") %>% #Only select dates that have Neolithic tag
    dplyr::arrange(median_caldates, x) %>% 
    dplyr::filter(between(median_caldates, 0, 7600)) %>%  #Allowable dates for first date of Neolithic
    dplyr::mutate(duplicated = dplyr::if_else(median_caldates == dplyr::lag(median_caldates) & x == dplyr::lag(x) & y == dplyr::lag(y), 1, 0)) %>% #Step need to prevent same date and location and singular vario model
    dplyr::filter(duplicated != 1)
  
  buff <- 50000 #Set buffer size
  
  p_neo_data_xy <- p_neo_data %>% #Generate x y of neo dates
    dplyr::select(x, y) 
  
  p_iberia_map_exp_25830_rast_small <- raster::raster(iberia_map_exp_25830, res = c(1000, 1000)) #Change interpolation resolution
  p_neo_data_rast <- raster::rasterize(p_neo_data_xy, p_iberia_map_exp_25830_rast_small, fun = max, field = p_neo_data$median_caldates)  #Place xy neo into raster with oldest neo date, 1km2
  p_extract_buffer <- raster::extract(p_neo_data_rast, p_neo_data_xy, buffer = buff, fun = max, df = T) #Establish buffer around each point to show max in buffer
  
  p_neo_data_xyz <- p_neo_data %>% 
    dplyr::mutate(buffer_extract = p_extract_buffer$layer) %>%   #For each sample, establish oldest date in its own buffer
    dplyr::mutate(test = dplyr::if_else(buffer_extract == median_caldates, "max", "not_max")) %>%     
    dplyr::filter(test == "max") %>%   #Only select samples that are the oldest in its own buffer
    dplyr::select(x, y, median_caldates) %>% 
    dplyr::rename(layer = median_caldates) %>% 
    sf::st_as_sf(coords = c("x", "y"), crs = 25830, agr = "constant") %>%  
    sf::as_Spatial()
  
  # Variogram
  p_vario <- gstat::variogram(layer~1, p_neo_data_xyz) #Generate variogram
  
  # Variogram model
  tryCatch( #To prevent variogram model failures stopping the loop
    p_vario_fit <- gstat::fit.variogram(p_vario, gstat::vgm(c("Exp")), fit.kappa = T), #Fit variogram
    error = function(e){
      p_vario_error <<- p_vario_error + 1
      p_vario_fit <<- p_vario_fit # Use same as previous
      p_neo_data_xyz <<- p_neo_data_xyz_previous # Use same as previous
    })
  
  p_neo_data_xyz_previous <- p_neo_data_xyz #required if need to use previous as otherwise the kriging will not produce  values
  
  # Kriging
  p_krige_neo <- gstat::krige(layer~1, p_neo_data_xyz, as(p_iberia_map_exp_25830_rast_small, "SpatialPixels"), model=p_vario_fit) #Perform kriging based on model and interpolation raster
  p_krige_neo_raster <- raster::raster(p_krige_neo) %>%  
    raster::mask(iberia_map_exp_25830)
  
  #Now extract Neolithic date for charcoal record locations
  p_entity_neo_check <- raster::extract(p_krige_neo_raster, entity_lon_lat_neo) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(entity_name = entity_list_neo$entity_name)
  
  if (sum(is.na(p_entity_neo_check$value)) > 0  ){
    p_NA_entity_neo <- p_entity_neo_check %>% #Identify NAs
      dplyr::filter(is.na(value)) %>% 
      dplyr::left_join(dplyr::select(entity_list_neo, entity_name, latitude, longitude), by = "entity_name")
    
    p_entity_neo_lon_lat_NA <- p_NA_entity_neo%>% #Get lon lat of NAs
      sf::st_as_sf(., coords = c("longitude", "latitude"), crs = 4258, agr = "constant") %>% 
      sf::st_transform(crs = 25830) %>% 
      sf::st_coordinates()
    
    p_entity_neo_NA <- raster::extract(p_krige_neo_raster, p_entity_neo_lon_lat_NA, buffer = 100000, fun = mean, df = T) %>% #Use buffer to get approximate values relating to NA   
      dplyr::as_tibble() %>% 
      dplyr::rename(value = 2) %>% 
      dplyr::mutate(entity_name = p_NA_entity_neo$entity_name) %>% 
      dplyr::select(value, entity_name) 
    
    p_entity_neo <- p_entity_neo_check %>% #Merge na and non-na values
      dplyr::filter(!(is.na(value))) %>% 
      rbind(p_entity_neo_NA) %>% 
      dplyr::rename(first_neo = value) %>% 
      dplyr::left_join(dplyr::select(entity_list_neo, entity_name, ID_ENTITY), by = "entity_name")
    
  } else {
    p_entity_neo <- p_entity_neo_check %>% #Merge na and non-na values
      dplyr::rename(first_neo = value) %>% 
      dplyr::left_join(dplyr::select(entity_list_neo, entity_name, ID_ENTITY), by = "entity_name")
  }
  
  
 
  
  #5.1 Pre-bin Neolithic amended data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
  
  datapath <- "./other_output/rpd/"
  p_neo_sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv") 
  p_neo_window <- 800 # Set the permitted window for analysis
  
  # Set baseperiod - here you will use this to select sites based on the degree of overlap with the baseperiod
  basebeg <- 200
  baseend <- 8200
  basename <- paste0("zt", basebeg, "_", baseend, "k")
  
  p_neo_dat <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>%
    dplyr::left_join(p_entity_neo, by = "ID_ENTITY") %>% #Insert extracted dates of first agriculture
    dplyr::mutate(EST_AGE = (EST_AGE - first_neo)) %>% #Re-centre dates based on first agriculture per site location
    dplyr::select(-c(entity_name, first_neo))
  
  p_neo_min_est_age <- -neo_window # To set ages later
  p_neo_max_est_age <- neo_window # To set ages later
  
  length(unique(p_neo_dat$ID_ENTITY))  ##all entities that were transformed  in long format
  p_neo_ent_vec2 <- unique(p_neo_dat$ID_ENTITY)
  
  p_neo_lsover <- p_neo_ent_vec2;  p_neo_lsover #This line of code enables subsetting of dataset if required (not utilised here)
  
  
  # Binning Z scores
  
  p_neo_Zbinning <- function(dt, bins, binhw){
    # Reshape data
    dat_res <- dt[,c("ID_ENTITY", "zt")] %>% 
      dplyr::group_by(ID_ENTITY) %>% 
      dplyr::mutate(number = row_number()) %>% 
      tidyr::spread(ID_ENTITY, zt) 
    dat_res <- dat_res[,-1] 
    dat_res <- as.matrix(dat_res)
    
    Age_res <- dt[,c("ID_ENTITY", "EST_AGE")] %>% 
      dplyr::group_by(ID_ENTITY) %>% 
      dplyr::mutate(number = row_number()) %>% 
      tidyr::spread(ID_ENTITY, EST_AGE)
    Age_res <- Age_res[,-1]
    Age_res <- as.matrix(Age_res)
    
    # Binning
    result <- matrix(ncol = length(Age_res[1,]), nrow = length(bins))
    colnames(result) <- colnames(dat_res) ##added this in so that column names of the binned matrix are entity IDs
    for (k in 1:length(dat_res[1,])){
      if(length(dat_res[is.na(dat_res[,k]) == F, k]) != 0){
        for (i in 1:length(bins)){
          t <- na.omit(cbind(as.numeric(Age_res[,k]), as.numeric(dat_res[,k])))
          result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] < bins[i] + binhw, 2])
        }
      }
    }
    p_neo_BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
    return(p_neo_BinnedData)
  }
  
  p_neo_bw <- 100 #Standard for composite
  
  
  p_neo_Zbin_dat <- p_neo_Zbinning(dt = p_neo_dat, bins = seq(p_neo_min_est_age, p_neo_max_est_age, p_neo_bw), binhw = p_neo_bw/2) 
  
  # Give ages of each bin as a column
  p_neo_bin_age <- seq(neo_min_est_age,p_neo_max_est_age,p_neo_bw)
  p_neo_Zbin_dat <- cbind(p_neo_bin_age, p_neo_Zbin_dat)
  
  
  # Limit this to the sites in lsover
  class(p_neo_lsover)
  p_neo_lsoverch <- as.character(p_neo_lsover)
  class(p_neo_lsoverch)
  class(colnames(p_neo_Zbin_dat))
  p_neo_Zbin_dat <- p_neo_Zbin_dat[,c('p_neo_bin_age',p_neo_lsoverch)]
  check <- colnames(p_neo_Zbin_dat)
  check[!(check %in% p_neo_lsoverch)]
  p_neo_lsoverch
  
  # Dataframe for curve generation
  p_neo_Zbin_dat <- as.data.frame(p_neo_Zbin_dat)
  
  
  # 5.2 Neolithic amended composite curve (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
  
  ## Set up
  
  # Set names
  p_neo_basename <- paste0("zt", basebeg, "_", baseend, "k")
  p_neo_binname <- paste0("bw", p_neo_bw)
  
  # Locfit (half) window-width parameter
  p_neo_hw <- 500 # bandwidth (smoothing parameter)
  
  
  # Target ages for fitted values
  p_neo_targbeg <- p_neo_min_est_age 
  p_neo_targend <- p_neo_max_est_age
  p_neo_targstep <- p_neo_bw ## this is based on the previous binning procedure
  
  # Array sizes
  maxrecs <- 2000
  maxreps <- 1000
  
  
  
  ## Neo curve plotting
  
  # Curve output path and file
  neo_curvecsvpath <- "./figs/neo_rpd/"
  
  #Read the list of sites to be processed
  p_neo_ns <- as.numeric(colnames(p_neo_Zbin_dat))
  p_neo_ns <- na.omit(p_neo_ns) 
  
  # Arrays for data and fitted values
  p_neo_age <- matrix(NA, ncol=length(p_neo_ns), nrow=maxrecs)
  p_neo_influx <- matrix(NA, ncol=length(p_neo_ns), nrow=maxrecs)
  p_neo_nsamples <- rep(0, maxrecs)
  p_neo_targage <- seq(p_neo_targbeg,p_neo_targend,p_neo_targstep) 
  p_neo_targage.df <- data.frame(p_neo_x=p_neo_targage)
  p_neo_lowage <- p_neo_targage - p_neo_hw; p_neo_highage <- p_neo_targage + p_neo_hw
  p_neo_ntarg <- length(p_neo_targage)
  p_neo_yfit <- matrix(NA, nrow=length(p_neo_targage.df$p_neo_x), ncol=maxreps) 
  
  # Arrays for sample number and effective window span tracking
  p_neo_ndec <- matrix(0, ncol=p_neo_ntarg, nrow=length(p_neo_ns))
  p_neo_ndec_tot <- rep(0, p_neo_ntarg)
  p_neo_xspan <- rep(0, p_neo_ntarg)
  p_neo_ninwin <- matrix(0, ncol=p_neo_ntarg, nrow=length(p_neo_ns))
  p_neo_ninwin_tot <- rep(0, p_neo_ntarg)
  
  
  ## Read data
  ii <- 0
  
  for(i in 2:length(p_neo_Zbin_dat)){ 
    p_neo_indata <- p_neo_Zbin_dat[,c(i,1)] 
    p_neo_indata <- na.omit(p_neo_indata)  
    p_neo_nsamp <-  length(p_neo_indata$p_neo_bin_age) 
    if (p_neo_nsamp > 0) {
      ii <- ii+1
      p_neo_age[1:p_neo_nsamp,ii] <- p_neo_indata$p_neo_bin_age 
      p_neo_influx[1:p_neo_nsamp,ii] <- p_neo_indata[,1]  
      p_neo_nsamples[ii] <- p_neo_nsamp
    }
    print(i)
  }
  
  p_neo_nsites <- ii; p_neo_nsites # number of sites with data
  
  # Trim samples to age range
  p_neo_influx[p_neo_age >= p_neo_targend+p_neo_hw] <- NA
  p_neo_age[p_neo_age >= p_neo_targend+p_neo_hw] <- NA 
  
  # Censor abs(influx) values > 10
  p_neo_influx[abs(p_neo_influx) >= 10] <- NA
  p_neo_age[abs(p_neo_influx) >= 10] <- NA 
  
  ## Find number of sites and samples contributing to fitted values
  
  # Count number of sites that contributed to each fitted value
  ptm <- proc.time()
  for (i in 1:p_neo_ntarg) {
    p_neo_agemax <- -1e32; p_neo_agemin <- 1e32
    for (j in 1:p_neo_nsites) {
      for (k in 1:p_neo_nsamples[j]) {
        if (!is.na(p_neo_age[k,j])) {
          ii <- (p_neo_age[k,j]-p_neo_targage[1])/p_neo_targstep + 1
          if (ii > 0 && ii <= p_neo_ntarg) {p_neo_ndec[j,ii] = 1}
          if (p_neo_age[k,j] >= p_neo_targage[i]-p_neo_hw && p_neo_age[k,j] <= p_neo_targage[i]+p_neo_hw) {
            p_neo_ninwin[j,i] = 1
            if (p_neo_agemax < p_neo_age[k,j]) {p_neo_agemax <- p_neo_age[k,j]}
            if (p_neo_agemin > p_neo_age[k,j]) {p_neo_agemin <- p_neo_age[k,j]}
          }
        }
      }
    }
    p_neo_ndec_tot[i] <- sum(p_neo_ndec[,i])
    p_neo_ninwin_tot[i] <- sum(p_neo_ninwin[,i])
    p_neo_xspan[i] <- p_neo_agemax - p_neo_agemin
  }
  proc.time() - ptm
  

  
  # Curve-fitting
  ptm <- proc.time()
  
  # Reshape matrices into vectors 
  p_neo_x <- as.vector(p_neo_age)
  p_neo_y <- as.vector(p_neo_influx)
  p_neo_lfdata <- data.frame(p_neo_x,p_neo_y)
  p_neo_lfdata <- na.omit(p_neo_lfdata)
  p_neo_x <- p_neo_lfdata$p_neo_x; p_neo_y <- p_neo_lfdata$p_neo_y
  
  # Locfit
  p_neo_loc01 <- locfit::locfit(p_neo_y ~ locfit::lp(p_neo_x, deg=1, h=p_neo_hw), maxk=800, family="qrgauss")
  summary(p_neo_loc01)
 
  # Get fitted values
  p_neo_pred01 <- predict(p_neo_loc01, newdata=p_neo_targage.df, se.fit=TRUE)
  p_neo_loc01_fit <- data.frame(p_neo_targage.df$p_neo_x, p_neo_pred01$fit)
  p_neo_fitname <- paste("p_neo_locfit_",as.character(p_neo_hw), sep="")
  colnames(p_neo_loc01_fit) <- c("p_neo_age", p_neo_fitname)
  head(p_neo_loc01_fit)
  
  #Detrend fitted values
  fit_p_neo_loc01_fit <- lm(p_neo_locfit_500 ~ p_neo_age, data = p_neo_loc01_fit)
  detrend_p_neo_loc01_fit <- dplyr::tibble(p_neo_age = p_neo_loc01_fit$p_neo_age, detrend_p_neo_locfit_500 = resid(fit_p_neo_loc01_fit)) %>% 
    as.data.frame()
  
  proc.time() - ptm
  
  
  # Plot the unresampled (initial) fit
  p_neo_fitname <- paste("p_neo_locfit_",as.character(p_neo_hw), sep="")
  colnames(p_neo_loc01_fit) <- c("p_neo_age", neo_fitname)
  
  
  #Output
  p_neo_curveout <- data.frame(cbind(p_neo_targage.df$p_neo_x, resid(fit_p_neo_loc01_fit), neo_ndec_tot, neo_xspan, neo_ninwin_tot))
  colnames(p_neo_curveout) <- c("age", paste0(q, ".locfit_", q), "nsites", "window", "ninwin")
  p_neo_curveout
  
  p_all_neo_curveout <- cbind(p_all_neo_curveout, p_neo_curveout) #Store output
  p_all_neo_curveout_loc <- cbind(p_all_neo_curveout_loc, dplyr::select(p_neo_curveout, 2)) #Store output
  
  warnings()
  
  
  
  
}


# Store and arrange data for permutation plotting
p_neo_loc_95 <- dplyr::tibble(name = "p_neo_loc_95", value = apply(p_all_neo_curveout_loc, 1, function(x) quantile(x,prob=0.975, na.rm=T))) #Calculate 95% CI
p_neo_loc_05 <- dplyr::tibble(name = "p_neo_loc_05", value = apply(p_all_neo_curveout_loc, 1, function(x) quantile(x,prob=0.025, na.rm=T))) #Calculate 95% CI
p_neo_loc <- dplyr::tibble(name = "p_neo_loc", value = detrend_neo_loc01_fit$detrend_neo_locfit_500)

p_neo_loc_plot <- tidyr::pivot_longer(p_all_neo_curveout_loc, cols = everything()) %>% 
  dplyr::arrange(name) %>% 
  rbind(p_neo_loc) %>% 
  rbind(p_neo_loc_95) %>% 
  rbind(p_neo_loc_05) %>%       
  dplyr::mutate(age = rep(p_all_neo_curveout$age, times = (ncol(p_all_neo_curveout_loc)+3))) %>%  
  dplyr::mutate(size = c(rep(1, times = 17*(ncol(p_all_neo_curveout_loc))), rep(2, times = 17), rep(1, times = 17), rep(1, times = 17))) %>% 
  dplyr::mutate(alpha = c(rep(0.3, times = 17*(ncol(p_all_neo_curveout_loc))), rep(2, times = 17), rep(1, times = 17), rep(1, times = 17)))


# Plot showing ranges of output from small variations in Neolthic tagged calibration dates
ggplot2::ggplot(data = p_neo_loc_plot, mapping = aes(y = value, x = age, color = name, size = size, alpha = alpha))+
  geom_line()+
  scale_color_manual(values = c(rep(c("gray"), times = ncol(p_all_neo_curveout_loc)), "yellow", "blue", "blue"))+
  scale_size(range = c(1,2), guide = F)+
  scale_alpha(range = c(0.3,2), guide = F)+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(750,-750), ylim = c(-1, 1))+
  scale_x_continuous(breaks = c(750, 500, 250, 0, -250, -500, -750), labels = c(-750, -500, -250, 0, 250, 500, 750))+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")+
  annotate("text", x = -300, y = 0.7, label = "Neolithic start date", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = 650, y = -0.7, label = "Before", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = -600, y = -0.7, label = "After", color = "blue", size = 12*0.36, family = "Times")+
  ggplot2::labs(x = "Age before and after the start of the Neolithic", y = "Detrended charcoal z-score composite")+
  ggsave("figs/neo_rpd/Neo_permute.pdf", height = 8, width = 11.28, units = "cm")



# ---------------------------------------------------------



# 6. ARIMA intervention analysis
# ---------------------------------------------------------
# ---------------------------------------------------------


# Set up file for analysis, using Neo locfit output
neo_loc_plot_annual_values <- neo_loc01_fit  %>% 
  dplyr::rename(loc = 2) %>% 
  dplyr::mutate(age = round(neo_age)) %>% 
  dplyr::filter(dplyr::between(age, -200, 750)) %>% 
  dplyr::select(age, loc)

neo_loc_plot_annual <- dplyr::tibble(age = seq(-200, 750, 1)) %>% 
  dplyr::left_join(neo_loc_plot_annual_values, by = "age") %>% 
  dplyr::mutate(loc_annual = zoo::na.spline(loc)) %>% 
  dplyr::mutate(age = -age) %>% 
  dplyr::arrange(rev(age)) 

# Select approriate data to tune model
ts_neo_loc_plot_annual <- tsibble::as_tsibble(neo_loc_plot_annual, key = NULL, index = age) %>%  
  dplyr::filter(age<1)  %>% 
  dplyr::select(-loc)



char_model_arima <- ts_neo_loc_plot_annual %>% #Generate ARIMA model based on data prior to Neolithic
  fabletools::model(fable::ARIMA(loc_annual))

fabletools::report(char_model_arima) #Identify the best ARIMA model

# Plot model forecast and actual output
char_model_arima %>%
  fabletools::forecast(h = 200) %>%  
  autoplot(ts_neo_loc_plot_annual) + 
  geom_line(data = neo_loc_plot_annual, mapping = aes(y = loc_annual, x = age))+
  ggplot2::labs(x = "Age before and after the start of the Neolithic", y = "Z-score composite of charcoal influx")+
  scale_x_continuous(breaks = c(-800, -600, -400, -200, 0, 200))+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")+
  annotate("text", x = -300, y = 0.125, label = "Neolithic start date", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = 100, y = -0.7, label = "After", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = -600, y = -0.7, label = "Before", color = "blue", size = 12*0.36, family = "Times")+
  ggplot2::labs(x = "Age before and after the start of the Neolithic", y = "Z-score composite of charcoal influx")+
  ggsave("figs/neo_rpd/Neo_ARIMA.pdf", height = 8, width = 11.28, units = "cm")

