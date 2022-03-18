#Sweeney, Harrison and Vander Linden 2021. 
#Assessing anthropogenic influence on fire history during the Holocene in the Iberian Peninsula
# ---------------------------------------------------------

#There are seven seperate scripts associated with this research:
# 1. data_import.R 
# 2. charcoal_analysis.R 
# 3. radiocarbon_analysis.R
# 4. correlation_analysis.R
# 5. SEA_analysis.R (this script)
# 6. neolithic_analysis.R
# 7. maps.R


# ---------------------------------------------------------


# SEA analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script uses the data generated within the
# charcoal_analyis.R and radiocarbon_analysis.R scrips to 
# assess the potential influence of particular periods of 
# rapid population change on fire history during the 10-3.5k 
# period within mainland Iberia. 

# As with charcoal_analysis.R, the construction of charcoal
# composite curves are based on the code  published by Patrick 
# J. Bartlein. see "Analysis of the Global Charcoal Database
# – GCDv3". https://pjbartlein.github.io/GCDv3Analysis/
# The code  here allows for reproduction of the results
# generated within the corresponding paper.



# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Incorporate periods of rapid population growth
# 3. Pre-bin data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 4. SEA composite curve (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 5. Permutation work based on variable first and second epoch dates
# 6. ARIMA intervention analysis



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






# 2. Incorporate periods of rapid population growth
# ---------------------------------------------------------
# ---------------------------------------------------------

first_sea <- 7396 # Manually select break points based on SPD pop data
second_sea <- 5369 # Manually select break points based on roc pop data


# ---------------------------------------------------------





# 3. Pre-bin SEA data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

datapath <- "./other_output/rpd/"
sea_sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv") #Load-in influx data
sea_window <- 800 # Set range of data to be considered for each epoch. Note, these should not overlap otherwise relationship becomes blurred

# Set baseperiod - here you will use this to select sites based on the degree of overlap with the baseperiod
basebeg <- 200
baseend <- 8200
basename <- paste0("zt", basebeg, "_", baseend, "k")

sea_dat_first <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>% # Re-centre data based on first SEA
  dplyr::filter(between(EST_AGE, first_sea - sea_window, first_sea + sea_window)) %>% 
  dplyr::mutate(EST_AGE = EST_AGE - first_sea) %>% 
  dplyr::mutate(group = "7500")

sea_dat_second <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>% # Re-centre data based on second SEA
  dplyr::filter(between(EST_AGE, second_sea - sea_window, second_sea + sea_window)) %>% 
  dplyr::mutate(EST_AGE = EST_AGE - second_sea) %>% 
  dplyr::mutate(group = "5500")

sea_dat <- sea_dat_first %>% # Consolidate epoch data for run as per charcoal_analysis.R
  rbind(sea_dat_second) %>% 
  dplyr::arrange(ID_ENTITY, EST_AGE)

sea_min_est_age <- -sea_window # To set ages later
sea_max_est_age <- sea_window # To set ages later

length(unique(sea_dat$ID_ENTITY))  ##all entities that were transformed  in long format
sea_ent_vec2 <- unique(sea_dat$ID_ENTITY)

sea_lsover <-  sea_ent_vec2 #This line of code enables subsetting of dataset if required (not utilised here)



# Binning Z scores
sea_Zbinning <- function(dt, bins, binhw){
  # Reshape data
  dat_res <- dt[,c("ID_ENTITY", "zt")] %>% 
    dplyr::group_by(ID_ENTITY) %>% 
    dplyr::mutate(number = row_number()) %>% 
    tidyr::spread(ID_ENTITY, zt) 
  dat_res <- dat_res[,-1] ##takes away a 'number' column
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
        result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] <= bins[i] + binhw, 2])
      }
    }
  }
  sea_BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
  return(sea_BinnedData)
}

sea_bw <- 100 #Standard for composite

sea_Zbin_dat <- sea_Zbinning(dt = sea_dat, bins = seq(sea_min_est_age, sea_max_est_age, sea_bw), binhw = sea_bw/2) 

# Give ages of each bin as a column
sea_bin_age <- seq(sea_min_est_age, sea_max_est_age, sea_bw)
sea_Zbin_dat <- cbind(sea_bin_age, sea_Zbin_dat)


## Limit this to the sites in lsover
class(sea_lsover)
sea_lsoverch <- as.character(sea_lsover)
class(sea_lsoverch)
class(colnames(sea_Zbin_dat))
sea_Zbin_dat <- sea_Zbin_dat[,c('sea_bin_age',sea_lsoverch)]
check <- colnames(sea_Zbin_dat)
check[!(check %in% sea_lsoverch)]


# Dataframe for curve generation
sea_Zbin_dat <- as.data.frame(sea_Zbin_dat)


# ---------------------------------------------------------




# 4. SEA composite curve (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

## Set up

# Set names
sea_basename <- paste0("zt", basebeg, "_", baseend, "k")
sea_binname <- paste0("bw", sea_bw)

# Locfit (half) window-width parameter
sea_hw <- 500 # bandwidth (smoothing parameter)

# Number of bootstrap samples/replications
nreps <- 1000

# Target ages for fitted values
sea_targbeg <- sea_min_est_age
sea_targend <- sea_max_est_age
sea_targstep <- sea_bw ## this is based on the previous binning procedure

# Array sizes
maxrecs <- 2000
maxreps <- 1000

# Plot output
plotout <- "pdf" 


## SEA curve plotting

# Curve output path and file
sea_curvecsvpath <- "./figs/sea_rpd/"
sea_curvefile <- paste(sea_curvecsvpath,"sea_rpd",".csv", sep="") #adjust to particular data you use now. .
print(sea_curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { 
  sea_pdffile <- paste(sea_curvecsvpath,"sea_rpd", ".pdf", sep="")
  print(sea_pdffile)
}

# Read the list of sites to be processed
sea_ns <- as.numeric(colnames(sea_Zbin_dat))
sea_ns <- na.omit(sea_ns) ##Z score data. note that this excludes the age column and only takes entity ids. so be careful when using the length of this in loop below.

# Arrays for data and fitted values
sea_age <- matrix(NA, ncol=length(sea_ns), nrow=maxrecs)
sea_influx <- matrix(NA, ncol=length(sea_ns), nrow=maxrecs)
sea_nsamples <- rep(0, maxrecs)
sea_targage <- seq(sea_targbeg,sea_targend,sea_targstep) 
sea_targage.df <- data.frame(sea_x=sea_targage)
sea_lowage <- sea_targage - sea_hw; sea_highage <- sea_targage + sea_hw
sea_ntarg <- length(sea_targage)
sea_yfit <- matrix(NA, nrow=length(sea_targage.df$sea_x), ncol=maxreps) 

# Arrays for sample number and effective window span tracking
sea_ndec <- matrix(0, ncol=sea_ntarg, nrow=length(sea_ns))
sea_ndec_tot <- rep(0, sea_ntarg)
sea_xspan <- rep(0, sea_ntarg)
sea_ninwin <- matrix(0, ncol=sea_ntarg, nrow=length(sea_ns))
sea_ninwin_tot <- rep(0, sea_ntarg)


## Read data

# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0

for(i in 2:length(sea_Zbin_dat)){ 
  sea_indata <- sea_Zbin_dat[,c(i,1)] 
  sea_indata <- na.omit(sea_indata) 
  sea_nsamp <-  length(sea_indata$sea_bin_age)
  if (sea_nsamp > 0) {
    ii <- ii+1
    sea_age[1:sea_nsamp,ii] <- sea_indata$sea_bin_age 
    sea_influx[1:sea_nsamp,ii] <- sea_indata[,1] 
    sea_nsamples[ii] <- sea_nsamp
  }
  print(i)
}

sea_nsites <- ii; sea_nsites # Number of sites

# Trim samples to age range
sea_influx[sea_age >= sea_targend+sea_hw] <- NA
sea_age[sea_age >= sea_targend+sea_hw] <- NA 

# Censor abs(influx) values > 10
sea_influx[abs(sea_influx) >= 10] <- NA
sea_age[abs(sea_influx) >= 10] <- NA 


## Find number of sites and samples contributing to fitted values

# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:sea_ntarg) {
  sea_agemax <- -1e32; sea_agemin <- 1e32
  for (j in 1:sea_nsites) {
    for (k in 1:sea_nsamples[j]) {
      if (!is.na(sea_age[k,j])) {
        ii <- (sea_age[k,j]-sea_targage[1])/sea_targstep + 1
        #print (c(i,j,k,ii))
        if (ii > 0 && ii <= sea_ntarg) {sea_ndec[j,ii] = 1}
        if (sea_age[k,j] >= sea_targage[i]-sea_hw && sea_age[k,j] <= sea_targage[i]+sea_hw) {
          sea_ninwin[j,i] = 1
          if (sea_agemax < sea_age[k,j]) {sea_agemax <- sea_age[k,j]}
          if (sea_agemin > sea_age[k,j]) {sea_agemin <- sea_age[k,j]}
        }
      }
    }
  }
  sea_ndec_tot[i] <- sum(sea_ndec[,i])
  sea_ninwin_tot[i] <- sum(sea_ninwin[,i])
  sea_xspan[i] <- sea_agemax - sea_agemin
}
proc.time() - ptm


## Curve-fitting and bootstrapping

# Composite curve
ptm <- proc.time()

# Reshape matrices into vectors
sea_x <- as.vector(sea_age)
sea_y <- as.vector(sea_influx)
sea_lfdata <- data.frame(sea_x,sea_y)
sea_lfdata <- na.omit(sea_lfdata)
sea_x <- sea_lfdata$sea_x; sea_y <- sea_lfdata$sea_y

sea_lfdata_arima <- sea_lfdata

# Locfit
sea_loc01 <- locfit::locfit(sea_y ~ locfit::lp(sea_x, deg=1, h=sea_hw), maxk=800, family="qrgauss")
summary(sea_loc01)

# Get fitted values
sea_pred01 <- predict(sea_loc01, newdata=sea_targage.df, se.fit=TRUE)
sea_loc01_fit <- data.frame(sea_targage.df$sea_x, sea_pred01$fit)
sea_fitname <- paste("sea_locfit_",as.character(sea_hw), sep="")
colnames(sea_loc01_fit) <- c("sea_age", sea_fitname)
head(sea_loc01_fit)

proc.time() - ptm


#Bootstrap-by-site confidence intervals
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications
sea_curvecsvpath
sea_pdffile

if (plotout == "pdf") {pdf(file=sea_pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(sea_x, sea_y, xlab="Age before and after ~7400 and ~5400", ylab="Charcoal z-score composite", xlim=c(750,-750), ylim=c(-1,1), type="n", cex.lab = 1.2, mgp = c(2.4, 1, 0), xaxt = "n") ##make sure xlim is right.
axis(1, at = c(750, 500, 250, 0, -250, -500, -750), labels = c(-750, -500, -250, 0, 250, 500, 750))
grid()

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) { 
  print(i)
  sea_randsitenum <- sample(seq(1:sea_nsites), sea_nsites, replace=TRUE)
  sea_x <- as.vector(sea_age[,sea_randsitenum])
  sea_y <- as.vector(sea_influx[,sea_randsitenum])
  sea_lfdata <- data.frame(sea_x,sea_y)
  sea_lfdata <- na.omit(sea_lfdata)
  sea_x <- sea_lfdata$sea_x; sea_y <- sea_lfdata$sea_y
  sea_locboot <- locfit::locfit(sea_y ~ locfit::lp(sea_x, deg=1, h=sea_hw), maxk=800, maxit=20, family="qrgauss")
  sea_predboot <- predict(sea_locboot, newdata=sea_targage.df, se.fit=TRUE)
  sea_yfit[,i] <- sea_predboot$fit
  lines(sea_targage.df$sea_x, sea_yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}
warnings()

# Step 3 -- Plot the unresampled (initial) fit
sea_fitname <- paste("sea_locfit_",as.character(sea_hw), sep="")
colnames(sea_loc01_fit) <- c("sea_age", sea_fitname)
lines(sea_loc01_fit[,1], sea_loc01_fit[,2], lwd=2, col="red") 

# Step 4 -- Find and add bootstrap CIs
sea_yfit95 <- apply(sea_yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
sea_yfit05 <- apply(sea_yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(sea_targage.df$sea_x, sea_yfit95, lwd=1, col="red")
lines(sea_targage.df$sea_x, sea_yfit05, lwd=1, col="red")

abline(v = 0, col = "blue", lwd = 1, lty = 2)
text.default(x = -25, y = 0.7, "Rapid population growth", col = "blue", pos = 4, cex = 1)
text.default(x = 750, y = -1, "Before", col = "blue", pos = 4, cex = 1)
text.default(x = -550, y = -1, "After", col = "blue", pos = 4, cex = 1)

if (plotout == "pdf") {dev.off()}


## SEA curve plotting: Detrended

# Curve output path and file
sea_curvecsvpath <- "./figs/sea_rpd/"
sea_curvefile <- paste(sea_curvecsvpath,"sea_rpd_detrend",".csv", sep="") #adjust to particular data you use now. .
print(sea_curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { 
  sea_pdffile <- paste(sea_curvecsvpath,"sea_rpd_detrend", ".pdf", sep="")
  print(sea_pdffile)
}

# Read the list of sites to be processed
sea_ns <- as.numeric(colnames(sea_Zbin_dat))
sea_ns <- na.omit(sea_ns) ##Z score data. note that this excludes the age column and only takes entity ids. so be careful when using the length of this in loop below.

# Arrays for data and fitted values
sea_age <- matrix(NA, ncol=length(sea_ns), nrow=maxrecs)
sea_influx <- matrix(NA, ncol=length(sea_ns), nrow=maxrecs)
sea_nsamples <- rep(0, maxrecs)
sea_targage <- seq(sea_targbeg,sea_targend,sea_targstep) 
sea_targage.df <- data.frame(sea_x=sea_targage)
sea_lowage <- sea_targage - sea_hw; sea_highage <- sea_targage + sea_hw
sea_ntarg <- length(sea_targage)
sea_yfit <- matrix(NA, nrow=length(sea_targage.df$sea_x), ncol=maxreps) 

# Arrays for sample number and effective window span tracking
sea_ndec <- matrix(0, ncol=sea_ntarg, nrow=length(sea_ns))
sea_ndec_tot <- rep(0, sea_ntarg)
sea_xspan <- rep(0, sea_ntarg)
sea_ninwin <- matrix(0, ncol=sea_ntarg, nrow=length(sea_ns))
sea_ninwin_tot <- rep(0, sea_ntarg)


## Read data

# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0

for(i in 2:length(sea_Zbin_dat)){ 
  sea_indata <- sea_Zbin_dat[,c(i,1)] 
  sea_indata <- na.omit(sea_indata) 
  sea_nsamp <-  length(sea_indata$sea_bin_age)
  if (sea_nsamp > 0) {
    ii <- ii+1
    sea_age[1:sea_nsamp,ii] <- sea_indata$sea_bin_age 
    sea_influx[1:sea_nsamp,ii] <- sea_indata[,1] 
    sea_nsamples[ii] <- sea_nsamp
  }
  print(i)
}

sea_nsites <- ii; sea_nsites # Number of sites

# Trim samples to age range
sea_influx[sea_age >= sea_targend+sea_hw] <- NA
sea_age[sea_age >= sea_targend+sea_hw] <- NA 

# Censor abs(influx) values > 10
sea_influx[abs(sea_influx) >= 10] <- NA
sea_age[abs(sea_influx) >= 10] <- NA 


## Find number of sites and samples contributing to fitted values

# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:sea_ntarg) {
  sea_agemax <- -1e32; sea_agemin <- 1e32
  for (j in 1:sea_nsites) {
    for (k in 1:sea_nsamples[j]) {
      if (!is.na(sea_age[k,j])) {
        ii <- (sea_age[k,j]-sea_targage[1])/sea_targstep + 1
        #print (c(i,j,k,ii))
        if (ii > 0 && ii <= sea_ntarg) {sea_ndec[j,ii] = 1}
        if (sea_age[k,j] >= sea_targage[i]-sea_hw && sea_age[k,j] <= sea_targage[i]+sea_hw) {
          sea_ninwin[j,i] = 1
          if (sea_agemax < sea_age[k,j]) {sea_agemax <- sea_age[k,j]}
          if (sea_agemin > sea_age[k,j]) {sea_agemin <- sea_age[k,j]}
        }
      }
    }
  }
  sea_ndec_tot[i] <- sum(sea_ndec[,i])
  sea_ninwin_tot[i] <- sum(sea_ninwin[,i])
  sea_xspan[i] <- sea_agemax - sea_agemin
}
proc.time() - ptm


## Curve-fitting and bootstrapping

# Composite curve
ptm <- proc.time()

# Reshape matrices into vectors
sea_x <- as.vector(sea_age)
sea_y <- as.vector(sea_influx)
sea_lfdata <- data.frame(sea_x,sea_y)
sea_lfdata <- na.omit(sea_lfdata)
sea_x <- sea_lfdata$sea_x; sea_y <- sea_lfdata$sea_y

sea_lfdata_arima <- sea_lfdata

# Locfit
sea_loc01 <- locfit::locfit(sea_y ~ locfit::lp(sea_x, deg=1, h=sea_hw), maxk=800, family="qrgauss")
summary(sea_loc01)

# Get fitted values
sea_pred01 <- predict(sea_loc01, newdata=sea_targage.df, se.fit=TRUE)
sea_loc01_fit <- data.frame(sea_targage.df$sea_x, sea_pred01$fit)
sea_fitname <- paste("sea_locfit_",as.character(sea_hw), sep="")
colnames(sea_loc01_fit) <- c("sea_age", sea_fitname)
head(sea_loc01_fit)

#Detrend fitted values
fit_sea_loc01_fit <- lm(sea_locfit_500 ~ sea_age, data = sea_loc01_fit)
detrend_sea_loc01_fit <- dplyr::tibble(sea_age = sea_loc01_fit$sea_age, detrend_sea_locfit_500 = resid(fit_sea_loc01_fit)) %>% 
  as.data.frame()

proc.time() - ptm


#Bootstrap-by-site confidence intervals
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications
sea_curvecsvpath
sea_pdffile

if (plotout == "pdf") {pdf(file=sea_pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(sea_x, sea_y, xlab="Age before and after ~7400 and ~5400", ylab="Detrended charcoal z-score composite", xlim=c(750,-750), ylim=c(-1,1), type="n", cex.lab = 1.2, mgp = c(2.4, 1, 0), xaxt = "n") ##make sure xlim is right.
axis(1, at = c(750, 500, 250, 0, -250, -500, -750), labels = c(-750, -500, -250, 0, 250, 500, 750))
grid()

# Step 2 -- Do the bootstrap iterations, and plot each detrended composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) { 
  print(i)
  sea_randsitenum <- sample(seq(1:sea_nsites), sea_nsites, replace=TRUE)
  sea_x <- as.vector(sea_age[,sea_randsitenum])
  sea_y <- as.vector(sea_influx[,sea_randsitenum])
  sea_lfdata <- data.frame(sea_x,sea_y)
  sea_lfdata <- na.omit(sea_lfdata)
  sea_x <- sea_lfdata$sea_x; sea_y <- sea_lfdata$sea_y
  sea_locboot <- locfit::locfit(sea_y ~ locfit::lp(sea_x, deg=1, h=sea_hw), maxk=800, maxit=20, family="qrgauss")
  sea_predboot <- predict(sea_locboot, newdata=sea_targage.df, se.fit=TRUE)
  sea_boot <- data.frame(sea_predboot$fit, sea_targage.df$sea_x) #Detrend
  fit_sea_boot <- lm(sea_predboot$fit ~ sea_targage.df.sea_x, data = sea_boot) #Detrend
  sea_yfit[,i] <- resid(fit_sea_boot)
  lines(sea_targage.df$sea_x, sea_yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}
warnings()

# Step 3 -- Plot the unresampled (initial) fit
sea_fitname <- paste("sea_detrend_locfit_",as.character(sea_hw), sep="")
colnames(detrend_sea_loc01_fit) <- c("sea_age", sea_fitname)
lines(detrend_sea_loc01_fit[,1], detrend_sea_loc01_fit[,2], lwd=2, col="red") 

# Step 4 -- Find and add bootstrap CIs
sea_yfit95 <- apply(sea_yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
sea_yfit05 <- apply(sea_yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(sea_targage.df$sea_x, sea_yfit95, lwd=1, col="red")
lines(sea_targage.df$sea_x, sea_yfit05, lwd=1, col="red")

abline(v = 0, col = "blue", lwd = 1, lty = 2)
text.default(x = -25, y = 0.7, "Rapid population growth", col = "blue", pos = 4, cex = 1)
text.default(x = 750, y = -1, "Before", col = "blue", pos = 4, cex = 1)
text.default(x = -550, y = -1, "After", col = "blue", pos = 4, cex = 1)

if (plotout == "pdf") {dev.off()}




# ---------------------------------------------------------



# 5. Permutation work based on variable first and second epoch dates (based on non-detrended data)
# ---------------------------------------------------------
# ---------------------------------------------------------

min_range_sea <- 200 #to set the range of values permitted plus or minus epoch date
max_range_sea <- 200 
nreps_sea <- 1000 #how many times a date is selected

p_all_sea_curveout_both <- data.frame(matrix(, nrow=17, ncol=0)) #Set empty dataframe for loop
p_all_sea_curveout_loc_both <- data.frame(matrix(, nrow=17, ncol=0)) #Set empty dataframe for loop

for (z in 1:nreps_sea){ # Loop to re-run binning and composite curve generation (as above)
  
  
  first_sea_rand <- round(runif(1, min = first_sea - min_range_sea, max = first_sea + max_range_sea)) #Randomly select first epoch date +-200 years of actual date
  second_sea_rand <- round(runif(1, min = second_sea - min_range_sea, max = second_sea + max_range_sea)) #Randomly select second epoch date +-200 years of actual date
  
  
  
  # 5.1 Pre-bin data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
  
  datapath <- "./other_output/rpd/"
  p_sea_sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv")
  p_sea_window_both <- 800 # Set the permitted window for analysis
  
  #set baseperiod - here you will use this to select sites based on the degree of overlap with the baseperiod, if you want to
  basebeg <- 200
  baseend <- 8200
  basename <- paste0("zt", basebeg, "_", baseend, "k")
  
  p_sea_dat_first_both <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>% 
    dplyr::filter(between(EST_AGE, first_sea_rand - p_sea_window_both, first_sea_rand + p_sea_window_both)) %>% 
    dplyr::mutate(EST_AGE = EST_AGE - first_sea_rand) %>% 
    dplyr::mutate(group = "7500") %>% 
    dplyr::arrange(ID_ENTITY, EST_AGE)
  
  p_sea_dat_second_both <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = '')) %>% 
    dplyr::filter(between(EST_AGE, second_sea_rand - p_sea_window_both, second_sea_rand + p_sea_window_both)) %>% 
    dplyr::mutate(EST_AGE = EST_AGE - second_sea_rand) %>% 
    dplyr::mutate(group = "5400") %>% 
    dplyr::arrange(ID_ENTITY, EST_AGE)
  
  p_sea_dat_both <- p_sea_dat_first_both %>% 
    rbind(p_sea_dat_second_both) %>% 
    dplyr::arrange(ID_ENTITY, EST_AGE)
  
  p_sea_min_est_age <- -p_sea_window_both
  p_sea_max_est_age <- p_sea_window_both
  p_sea_ent_vec2 <- unique(p_sea_dat_both$ID_ENTITY)

  p_sea_lsover <- p_sea_ent_vec2 #This line of code enables subsetting of dataset if required (not utilised here)
  
  
  # Binning Z scores
  p_sea_Zbinning <- function(dt, bins, binhw){
    # Reshape data
    dat_res <- dt[,c("ID_ENTITY", "zt")] %>% 
      dplyr::group_by(ID_ENTITY) %>% 
      dplyr::mutate(number = row_number()) %>% 
      tidyr::spread(ID_ENTITY, zt) 
    dat_res <- dat_res[,-1] ##takes away a 'number' column
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
          result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] <= bins[i] + binhw, 2])
        }
      }
    }
    p_sea_BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
    return(p_sea_BinnedData)
  }
  
  p_sea_bw <- 100 #Standard for composite
  
  p_sea_Zbin_dat <- p_sea_Zbinning(dt = p_sea_dat_both, bins = seq(p_sea_min_est_age, p_sea_max_est_age, p_sea_bw), binhw = p_sea_bw/2) 
  
  # Give ages of each bin as a column
  p_sea_bin_age <- seq(p_sea_min_est_age, p_sea_max_est_age, p_sea_bw)
  
  p_sea_Zbin_dat <- cbind(p_sea_bin_age, p_sea_Zbin_dat)
  
  ## Limit this to the sites in lsover
  class(p_sea_lsover)
  p_sea_lsoverch <- as.character(p_sea_lsover)
  class(p_sea_lsoverch)
  class(colnames(p_sea_Zbin_dat))
  p_sea_Zbin_dat <- p_sea_Zbin_dat[,c('p_sea_bin_age',p_sea_lsoverch)]
  check <- colnames(p_sea_Zbin_dat)
  check[!(check %in% p_sea_lsoverch)]
  
  
  # Dataframe for curve generation
  p_sea_Zbin_dat <- as.data.frame(p_sea_Zbin_dat)
  
  
  # 5.2 SEA composite data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)

  ## Set up

  # Names
  p_sea_basename <- paste0("zt", basebeg, "_", baseend, "k")
  p_sea_binname <- paste0("bw", p_sea_bw)
  
  # Locfit (half) window-width parameter
  p_sea_hw <- 500 # bandwidth (smoothing parameter)
  
  # Target ages for fitted values
  p_sea_targbeg <- p_sea_min_est_age
  p_sea_targend <- p_sea_max_est_age
  p_sea_targstep <- p_sea_bw ## this is based on the previous binning procedure
  
  # Array sizes
  maxrecs <- 2000
  maxreps <- 1000
  
  # Producing and saving composite
  p_sea_curvecsvpath <- "./figs/sea_rpd/"
  p_sea_curvecsvpath
  
  # Read the list of sites
  p_sea_ns <- as.numeric(colnames(p_sea_Zbin_dat))
  p_sea_ns <- na.omit(p_sea_ns) 
  length(p_sea_ns)
  
  
  # Arrays for data and fitted values
  p_sea_age <- matrix(NA, ncol=length(p_sea_ns), nrow=maxrecs)
  p_sea_influx <- matrix(NA, ncol=length(p_sea_ns), nrow=maxrecs)
  p_sea_nsamples <- rep(0, maxrecs)
  p_sea_targage <- seq(p_sea_targbeg,p_sea_targend,p_sea_targstep) 
  p_sea_targage.df <- data.frame(p_sea_x=p_sea_targage)
  p_sea_lowage <- p_sea_targage - p_sea_hw; p_sea_highage <- p_sea_targage + p_sea_hw
  p_sea_ntarg <- length(p_sea_targage)
  p_sea_yfit <- matrix(NA, nrow=length(p_sea_targage.df$p_sea_x), ncol=maxreps) #the actual curve values (fitted)
  
  # Arrays for sample number and effective window span tracking
  p_sea_ndec <- matrix(0, ncol=p_sea_ntarg, nrow=length(p_sea_ns))
  p_sea_ndec_tot <- rep(0, p_sea_ntarg)
  p_sea_xspan <- rep(0, p_sea_ntarg)
  p_sea_ninwin <- matrix(0, ncol=p_sea_ntarg, nrow=length(p_sea_ns))
  p_sea_ninwin_tot <- rep(0, p_sea_ntarg)
  
  
  ## Read data
  # Read and store the presample (binned) files as matrices of ages and influx values
  ii <- 0
  
  for(i in 2:length(p_sea_Zbin_dat)){ 
    p_sea_indata <- p_sea_Zbin_dat[,c(i,1)] 
    p_sea_indata <- na.omit(p_sea_indata) 
    p_sea_nsamp <-  length(p_sea_indata$p_sea_bin_age)
    if (p_sea_nsamp > 0) {
      ii <- ii+1
      p_sea_age[1:p_sea_nsamp,ii] <- p_sea_indata$p_sea_bin_age 
      p_sea_influx[1:p_sea_nsamp,ii] <- p_sea_indata[,1] 
      p_sea_nsamples[ii] <- p_sea_nsamp
    }
    print(i)
  }
  
  p_sea_nsites <- ii; p_sea_nsites # Number of sites with data
  
  # Trim samples to age range
  p_sea_influx[p_sea_age >= p_sea_targend+p_sea_hw] <- NA
  p_sea_age[p_sea_age >= p_sea_targend+p_sea_hw] <- NA 
  
  # Censor abs(influx) values > 10
  p_sea_influx[abs(p_sea_influx) >= 10] <- NA
  p_sea_age[abs(p_sea_influx) >= 10] <- NA 
  
  ## Find number of sites and samples contributing to fitted values

  # Count number of sites that contributed to each fitted value
  ptm <- proc.time()
  for (i in 1:p_sea_ntarg) {
    p_sea_agemax <- -1e32; p_sea_agemin <- 1e32
    for (j in 1:p_sea_nsites) {
      for (k in 1:p_sea_nsamples[j]) {
        if (!is.na(p_sea_age[k,j])) {
          ii <- (p_sea_age[k,j]-p_sea_targage[1])/p_sea_targstep + 1
          if (ii > 0 && ii <= p_sea_ntarg) {p_sea_ndec[j,ii] = 1}
          if (p_sea_age[k,j] >= p_sea_targage[i]-p_sea_hw && p_sea_age[k,j] <= p_sea_targage[i]+p_sea_hw) {
            p_sea_ninwin[j,i] = 1
            if (p_sea_agemax < p_sea_age[k,j]) {p_sea_agemax <- p_sea_age[k,j]}
            if (p_sea_agemin > p_sea_age[k,j]) {p_sea_agemin <- p_sea_age[k,j]}
          }
        }
      }
    }
    p_sea_ndec_tot[i] <- sum(p_sea_ndec[,i])
    p_sea_ninwin_tot[i] <- sum(p_sea_ninwin[,i])
    p_sea_xspan[i] <- p_sea_agemax - p_sea_agemin
  }
  proc.time() - ptm
  

  
  
  #Curve-fitting and bootstrapping
  
  ptm <- proc.time()
  
  # Reshape matrices into vectors
  p_sea_x <- as.vector(p_sea_age)
  p_sea_y <- as.vector(p_sea_influx)
  p_sea_lfdata <- data.frame(p_sea_x,p_sea_y)
  p_sea_lfdata <- na.omit(p_sea_lfdata)
  p_sea_x <- p_sea_lfdata$p_sea_x; p_sea_y <- p_sea_lfdata$p_sea_y
  
  # Locfit
  p_sea_loc01 <- locfit::locfit(p_sea_y ~ locfit::lp(p_sea_x, deg=1, h=p_sea_hw), maxk=800, family="qrgauss")
  summary(p_sea_loc01)

  # Get  fitted values
  p_sea_pred01 <- predict(p_sea_loc01, newdata=p_sea_targage.df, se.fit=TRUE)
  p_sea_loc01_fit <- data.frame(p_sea_targage.df$p_sea_x, p_sea_pred01$fit)
  p_sea_fitname <- paste("p_sea_locfit_",as.character(p_sea_hw), sep="")
  colnames(p_sea_loc01_fit) <- c("p_sea_age", p_sea_fitname)
  head(p_sea_loc01_fit)
  
  #Detrend fitted values
  # fit_p_sea_loc01_fit <- lm(p_sea_locfit_500 ~ p_sea_age, data = p_sea_loc01_fit)
  # detrend_p_sea_loc01_fit <- dplyr::tibble(p_sea_age = p_sea_loc01_fit$p_sea_age, detrend_p_sea_locfit_500 = resid(fit_p_sea_loc01_fit)) %>% 
  #   as.data.frame()
  
  proc.time() - ptm
  
  p_sea_fitname <- paste("p_sea_locfit_",as.character(p_sea_hw), sep="")
  colnames(p_sea_loc01_fit) <- c("p_sea_age", p_sea_fitname)
 
  #Output
  # p_sea_curveout_both <- data.frame(cbind(p_sea_targage.df$p_sea_x, resid(fit_p_sea_loc01_fit), p_sea_ndec_tot, p_sea_xspan, p_sea_ninwin_tot))
  p_sea_curveout_both <- data.frame(cbind(p_sea_targage.df$p_sea_x, p_sea_loc01_fit$p_sea_locfit_500, p_sea_ndec_tot, p_sea_xspan, p_sea_ninwin_tot))
  # colnames(p_sea_curveout_both) <- c("age", paste0(z, ".detrend_locfit_", first_sea_rand, "_", second_sea_rand), "nsites", "window", "ninwin")
  colnames(p_sea_curveout_both) <- c("age", paste0(z, ".locfit_", first_sea_rand, "_", second_sea_rand), "nsites", "window", "ninwin")
  p_sea_curveout_both
  
  proc.time() - ptm
  
  
  p_all_sea_curveout_both <- cbind(p_all_sea_curveout_both, p_sea_curveout_both) #Store output
  p_all_sea_curveout_loc_both <- cbind(p_all_sea_curveout_loc_both, dplyr::select(p_sea_curveout_both, 2)) #Store output
  
  warnings()
  
  
  
}

# Store and arrange data for permutation plotting
p_loc_95_both <- dplyr::tibble(name = "p_loc_95", value = apply(p_all_sea_curveout_loc_both, 1, function(x) quantile(x,prob=0.975, na.rm=T))) #Calculate 95% CI
p_loc_05_both <- dplyr::tibble(name = "p_loc_05", value = apply(p_all_sea_curveout_loc_both, 1, function(x) quantile(x,prob=0.025, na.rm=T))) #Calculate 95% CI
# detrend_sea_loc <- dplyr::filter(detrend_sea_loc01_fit, dplyr::between(sea_age, -800, 800))
# p_sea_loc <- dplyr::tibble(name = "p_sea_loc", value = detrend_sea_loc$sea_detrend_locfit_500)
p_sea_loc <- dplyr::tibble(name = "p_sea_loc", value = sea_loc01_fit$sea_locfit_500) 

p_sea_plot_both <- tidyr::pivot_longer(p_all_sea_curveout_loc_both, cols = everything()) %>%   
  dplyr::arrange(name) %>% 
  rbind(p_sea_loc) %>% 
  rbind(p_loc_95_both) %>% 
  rbind(p_loc_05_both) %>%     
  dplyr::mutate(age = rep(p_all_sea_curveout_both$age, times = (ncol(p_all_sea_curveout_loc_both)+3))) %>% 
  dplyr::mutate(size = c(rep(1, times = 17*(ncol(p_all_sea_curveout_loc_both))), rep(2, times = 17), rep(1, times = 17), rep(1, times = 17))) %>% 
  dplyr::mutate(alpha = c(rep(0.3, times = 17*(ncol(p_all_sea_curveout_loc_both))), rep(2, times = 17), rep(1, times = 17), rep(1, times = 17))) 

# Plot showing ranges of output from small variations in epoch dates
ggplot2::ggplot(data = p_sea_plot_both, mapping = aes(y = value, x = age, color = name, size = size, alpha = alpha))+
  geom_line()+
  scale_color_manual(values = c(rep(c("gray"), times = ncol(p_all_sea_curveout_loc_both)), "blue", "blue", "yellow"))+
  scale_size(range = c(1,2), guide = F)+
  scale_alpha(range = c(0.1,2), guide = F)+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(750,-750), ylim = c(-1, 1))+
  scale_x_reverse(breaks = c(750, 500, 250, 0, -250, -500, -750), labels = c(-750, -500, -250, 0, 250, 500, 750))+
  labs(x ="Age before and after ~7400 and ~5400", y ="Charcoal z-score composite")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")+
  annotate("text", x = -450, y = 0.7, label = "Rapid population growth", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = 650, y = -0.7, label = "Before", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = -650, y = -0.7, label = "After", color = "blue", size = 12*0.36, family = "Times")
  ggsave("figs/sea_rpd/Sea_permute.pdf", height = 8, width = 11.28, units = "cm")



# ---------------------------------------------------------






# 6. ARIMA intervention analysis
# ---------------------------------------------------------
# ---------------------------------------------------------

# Set up file for analysis, using SEA locfit output
loc01_fit_annual_sea_adj <- dplyr::tibble(sea_age = -200:750) %>% 
  dplyr::left_join(sea_loc01_fit, by = "sea_age") %>% 
  dplyr::mutate(sea_loc01_fit_annual = zoo::na.spline(sea_locfit_500)) %>% #Interpolate annual
  dplyr::mutate(sea_age = -sea_age) %>%  #Because of the BP nature of timing 
  dplyr::arrange(rev(sea_age)) %>% 
  dplyr::rename(loc = sea_loc01_fit_annual)

ggplot2::ggplot(data = loc01_fit_annual_sea_adj, mapping = aes(y = loc, x = sea_age))+
  geom_line()

# Select appropriate data to tune model
ts_sea_locfit_annual <- tsibble::as_tsibble(loc01_fit_annual_sea_adj, key = NULL, index = sea_age) %>%  #Select data prior to epochs
  dplyr::filter(sea_age<1) %>% 
  dplyr::select(-sea_locfit_500)

ts_sea_locfit_annual %>% 
  feasts::gg_tsdisplay(difference(difference(loc)), plot_type='partial')

model_3 <- ts_sea_locfit_annual %>% #Generate ARIMA model based on data prior to epochs
  fabletools::model(fable::ARIMA(loc, trace = T, stepwise = F, approximation = F))

model_3 %>% 
  feasts::gg_tsresiduals()

fabletools::report(model_3) #Identify the best ARIMA model

# Plot model forecast and actual output
model_3 %>% 
  fabletools::forecast(h = 200) %>%   
  autoplot(ts_sea_locfit_annual) +
  geom_line(data = loc01_fit_annual_sea_adj, mapping = aes(y = loc, x = sea_age))+
  scale_x_continuous(limits = c(-500, 200), breaks = c(-500, -400, -300, -200, -100, 0, 100, 200))+
  ggplot2::labs(x = "Age before and after ~7400 and ~5400", y = "Z-score composite of charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue")+
  annotate("text", x = -200, y = 0.75, label = "Rapid population growth", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = 150, y = -0.345, label = "After", color = "blue", size = 12*0.36, family = "Times")+
  annotate("text", x = -450, y = -0.345, label = "Before", color = "blue", size = 12*0.36, family = "Times")
  ggsave("figs/sea_rpd/Sea_arima.pdf", height = 8, width = 13, units = "cm")





