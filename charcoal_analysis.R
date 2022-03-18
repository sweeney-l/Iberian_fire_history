#Sweeney, Harrison and Vander Linden 2021. 
#Assessing anthropogenic influence on fire history during the Holocene in the Iberian Peninsula
# ---------------------------------------------------------

#There are seven seperate scripts associated with this research:
# 1. data_import.R 
# 2. charcoal_analysis.R (the script)
# 3. radiocarbon_analysis.R
# 4. correlation_analysis.R
# 5. SEA_analysis.R
# 6. neolithic_analysis.R
# 7. maps.R


# ---------------------------------------------------------



# RPD charcoal data analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script sets up the charcoal data for subsequent 
# relationship analysis with the radiocarbon data in 
# fire_human_relationship_analysis.R. To generate the data
# used in this script, the data_import.R script should be
# run first.

# It is important to note that this code draws heavily on 
# but slightly adapts scripts already published by Patrick 
# J. Bartlein. see "Analysis of the Global Charcoal Database
# – GCDv3". https://pjbartlein.github.io/GCDv3Analysis/
# The code  here allows for reproduction of the results
# generated within the corresponding paper




# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Calculate influx values and generate individual record csv files (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 3. Transform data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 4. Pre-bin data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# 5. Composite curves (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
#             : A - Single curve
#             : B - Multiple curves with different smoothing windows
#             : C - Multiple curves with different bin widths
# 6. Save data for correlation analysis
# 7. Generate curves based on alternative specifications
# 8. Changepoint analysis
# 9. Trend in charcoal data 
# 10. Site information 

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



# Load in charcoal data
# ---------------------------------------------------------

# Import data generated by data_import.R
rpd_data_input <- rio::import("other_output/rpd/rpd_data_input.csv")

# ---------------------------------------------------------








# 2. Calculate influx values and generate individual record csv files (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

# Read the data file
data_query <- rpd_data_input # Load-in charcoal data
data_query_entity <- rpd_data_input %>% 
  dplyr::group_by(ID_ENTITY) %>% 
  dplyr::summarise(ID_SITE = mean(ID_SITE), 
                   entity_name = unique(entity_name),
                   latitude = mean(latitude),
                   longitude = mean(longitude),
                   elevation = mean(elevation),
                   TYPE = unique(TYPE))

# The main part of the script loops over the individual entities, and does various checks and calculations, including:
# calculation of sedimentation rates and deposition times
# checking for age or depth reversals;
# setting of various indicator variables and flags;
# calculation of alternative quanties (e.g. influx, given concentrations);
# further checking for anomalies
# writing out a .csv file for each site
# consolidating files into a single csv

# debug/log file
csvpath <- "./other_output/rpd/rpd_debug/"
debugname <- "RPD_influx_debug.txt"
# open the debug/log file
debugfile <- file(paste(csvpath, debugname, sep=""), "w")

ls1 <- ls() #To check which entities are successfully saved following transformation

records <- unique(data_query$ID_ENTITY) #Identify the unique records to be analysed

for (j in records) { # Loop over unique records
  nsamp <- 0
  sitedata <- data_query[data_query$ID_ENTITY == j, ]
  nsamp <- length(sitedata$ID_ENTITY)
  
  if (nsamp > 1) { # For those entities with more than one sample
    jchar <- as.character(j) #Record number
    nsampchar <- as.character(nsamp) #Number of samples
    writeLines(paste("Site",jchar,nsampchar,"samples", sep=" "), con = debugfile, sep = "\n") #For debug purposes
    
    # Define variables
    depth <- sitedata$sample_depth
    age <- sitedata$INTCAL2020_median
    quant <- sitedata$quantity
    type <- sitedata$TYPE 
    
    # Recode NA's to missing
    miss <- -9999.0
    depth[is.na(depth)] <- miss
    age[is.na(age)] <- miss
    quant[is.na(quant)] <- miss
    
    # Define new variables for sedimentation and dposition times
    thickness <- rep(miss, nsamp)
    dep_time <- rep(miss, nsamp)
    sed_rate <- rep(miss, nsamp)
    unit_dep_time <- rep(miss, nsamp) 
    
    # Sed rate and deposition time
    # First (top) sample
    if (depth[1] != miss && depth[2] != miss) {
      thickness[1] <- (depth[2] - depth[1])
      dep_time[1] <- age[2] - age[1]
      if (dep_time[1] > 0.0) sed_rate[1] <- thickness[1]/dep_time[1]
      if (sed_rate[1] != miss) unit_dep_time[1] <- 1.0/sed_rate[1] #the 1/sed rate - 1 is a placeholder because it will be *quant
    }
    # Samples 2 to nsamp-1
    for (i in 2:(nsamp-1)) {
      if (depth[1] != miss && depth[2] != miss) {
        thickness[i] <- (depth[i+1] - depth[i]) 
        dep_time[i] <- ((age[i+1] + age[i])/2.0) - ((age[i] + age[i-1])/2.0)
        if (dep_time[i] > 0.0) sed_rate[i] <- thickness[i]/dep_time[i]
        if (sed_rate[i] != miss) unit_dep_time[i] <- 1.0/sed_rate[i] 
      }
    }
    # Last (bottom) sample
    if (depth[nsamp-1] != miss  && depth[nsamp] != miss) {
      thickness[nsamp] <- thickness[nsamp-1] # Replicate thickness
      dep_time[nsamp] <- age[nsamp] - age[nsamp-1]
      sed_rate[nsamp] <- sed_rate[nsamp-1] # Replicate sed_rate
      unit_dep_time[nsamp] <- unit_dep_time[nsamp-1]
    }
    
    # Counts of missing values
    depth_count <- 0
    age_count <- 0
    quant_count <- 0
    sed_rate_count <- 0
    sed_rate_flag <- 1
    depth_count <- sum(depth != miss)
    age_count <- sum(age != miss)
    quant_count <- sum(quant != miss)
    sed_rate_count <- sum(sed_rate != miss)
    if (sed_rate_count != nsamp) sed_rateflag = 0
    
    # Check for age or depth reversals, and zero or negative sedimentation rates (in nonmissing data)
    depth_reversal <- 0
    age_reversal <- 0
    sed_rate_zeroneg <- 0         
    for (i in 2:nsamp) {
      if (age[i] != miss && age[i-1] != miss && age[i] <= age[i-1]) age_reversal=1
      if (depth[i] != miss && depth[i-1] != miss) {
        if (depth[i] <= depth[i-1]) depth_reversal=1
      } 
    }
    for (i in 2:nsamp) {
      if (sed_rate[i] != miss && sed_rate[i] <= 0.0) sed_rate_zeroneg=1
    }
    
    # Set and write out various flags within debug file
    if (depth_count != 0 && depth_count != nsamp) {
      writeLines(paste("**** has a missing depth when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
    }
    if (age_count != 0 && age_count != nsamp) {
      writeLines(paste("**** has a missing age when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
    }
    if (quant_count != 0 && quant_count != nsamp) {
      writeLines(paste("**** has a missing quantity when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
    }
    if (sed_rate_count != 0 && sed_rate_count != nsamp) {
      writeLines(paste("**** has a missing sed rate when some are nonmissing", sep=" "), con = debugfile, sep = "\n")
    }
    if (depth_reversal != 0) {
      writeLines(paste("**** has a depth reversal", sep=" "), con = debugfile, sep = "\n")
    }
    if (age_reversal != 0) {
      writeLines(paste("**** has an age reversal", sep=" "), con = debugfile, sep = "\n")
    }
    if (sed_rate_zeroneg != 0) {
      writeLines(paste("**** has zero or negative sed rates", sep=" "), con = debugfile, sep = "\n")
    }
    
    # Calculate alternative quantities for each record
    
    conc <- rep(miss, nsamp)
    influx <- rep(miss, nsamp)
    
    influx_source <- rep("none", nsamp)
    conc_source <- rep("none", nsamp)
    
    # Select case based on record type
    if (type[1] == "influx") {         # Adopt influx values as they are, calculate concentration
      influx <- quant
      influx_source <- "data"
      if (influx != miss && unit_dep_time != miss && sed_rate != 0.0) {
        conc <- influx * unit_dep_time
        conc_source <- "calculated from influx "
      } else {
        conc <- quant
        conc_source <- "copied from quant "
      }
      writeLines("influx", con = debugfile, sep = "\n")
    }
    
    else if (type[1] == "concentration") {    # Calculate influx, adopt concentration values as they are
      conc <- quant
      conc_source <- "data"
      if (conc != miss && sed_rate != miss && sed_rate != 0.0) {
        influx <- quant * sed_rate
        influx_source <- "calculated from conc "
      } else {
        influx <- quant
        influx_source <- "copied from quant "
      }  
      writeLines("concentration", con = debugfile, sep = "\n")
    }
    
    else if (type[1] == "pollen concentration"){     # Assume quantity is concentration-like
      conc <- quant
      conc_source <- "pollen concentration"
      if (sed_rate != miss && sed_rate != 0.0) {
        influx <- quant * sed_rate
        influx_source <- "calculated from C0P0 (conc) "
      } else {
        influx <- quant
        influx_source <- "copied from quant "
      }    
      writeLines("pollen concentration", con = debugfile, sep = "\n")
    }
    
      else if (type[1] == "per unit weight") {    # Assume quantity is concentration-like
      conc <- quant
      conc_source <- "per unit weight"
      if (sed_rate != miss && sed_rate != 0.0) {
        influx <- quant * sed_rate
        influx_source <- "calculated from per unit weight "
      } else {
        influx <- quant
        influx_source <- "copied from quant "
      }
      writeLines("per unit weight", con = debugfile, sep = "\n")
    }
  
    else {
      conc <- quant
      conc_source <- "copied from quant "
      influx <- quant
      influx_source <- "copied from quant "
      writeLines("Unknown", con = debugfile, sep = "\n") 
    }
  }
  
 
  # Check for records where influx is 0.0
  nzero <- 0
  nzero <- sum(influx != 0.0)
  if (nzero == 0) {
    writeLines(paste("**** has no non-zero influx values", sep=" "), con = debugfile, sep = "\n")
  }
  
  # Produce csv output for each record
  if (nsamp > 1 && nzero > 0) { 
    
    # Generate siteid string
    entname <- unique(sitedata$entity_name)
    entname <- substr(entname, start = 1, stop = 3)
    entidchar <- as.character(j)
    if (j >= 1) entid <- paste("000", entidchar, entname, sep="") ##this needs to be fixed to be suitable for entity names.
    if (j >= 10) entid <- paste("00", entidchar, entname, sep="")
    if (j >= 100) entid <- paste("0", entidchar, entname, sep="")
    if (j >= 1000) entid <- paste(    entidchar, entname, sep="")
    
    
    # Assemble output data and write it out, and save a list of entity IDs that get successfully written out
    ls1 <- append(ls1,entidchar)
    
    siteidvec <- sitedata$ID_SITE
    entidvec <- sitedata$ID_ENTITY
    sampidvec <- sitedata$ID_SAMPLE
    
    outdata <- data.frame(siteidvec, entidvec, sampidvec, depth, age, sed_rate, quant, conc, influx, type, conc_source, influx_source) #removed ID_SAMPLE from this block of code - can add back in if you want
    names(outdata) <- c("ID_SITE", "ID_ENTITY", "ID_SAMPLE", "depth", "est_age", "sed_rate", "quant", "conc",
                        "influx", "type", "conc_source", "influx_source" )
    csvfile <- paste("./other_output/rpd/rpd_influx/",entid,"_RPD_influx.csv", sep="") # Path for csv file
    write.csv(outdata, csvfile, row.names=FALSE)
  } 
  print(j)
}

close(debugfile)


influx_all <- list.files(path = "./other_output/rpd/rpd_influx/", pattern = "influx", full.names = TRUE) %>%    #Combine all influx files
  lapply(readr::read_csv) %>%  
  bind_rows() %>% 
  dplyr::left_join(dplyr::select(data_query_entity, ID_ENTITY, entity_name), by = "ID_ENTITY")

# Check to see whether any ID_entities werent conserved. 

checkdat <- outdata[0,]
for(z in records){
  #compose filename
  tmp <- data_query_entity[data_query_entity$ID_ENTITY == z,]
  entname <- unique(tmp$entity_name)
  entname <- substr(entname, start = 1, stop = 3)
  entidchar <- as.character(z)
  if (z >= 1) entid <- paste("000", entidchar, entname, sep="") ##this needs to be fixed to be suitable for entity names.
  if (z >= 10) entid <- paste("00", entidchar, entname, sep="")
  if (z >= 100) entid <- paste("0", entidchar, entname, sep="")
  if (z >= 1000) entid <- paste(    entidchar, entname, sep="")
  
  csvin <- paste("./other_output/rpd/rpd_influx/",entid,"_RPD_influx.csv", sep="") # the individual files you wrote out earlier
  a <- read.csv(csvin)
  checkdat <- rbind(checkdat, a)
}

ent_infl <- data_query_entity[data_query_entity$ID_ENTITY %in% checkdat$ID_ENTITY,]

write.csv(ent_infl, "./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv",row.names = F)

# ---------------------------------------------------------








# 3. Transform data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

# 1-parameter Box-Cox transformation of charcoal quantities for each entity
# (alpha (shift parameter) is specified, lambda (power transformation parameter) is estimated)

# input .csv files should contain at least the variables "est_age" and "quant",
# identified by labels in a header row  


datapath <- "./other_output/rpd/"
sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv") 

# Set base period ages for z-scores
basebeg <- 200
baseend <- 8200 #Set 8200 as the maximum to ensure data quantity maximised
basename <- paste0("zt", basebeg, "_", baseend, "k")

# Debug/log file
debugpath <- "./other_output/rpd/rpd_debug/"  
debugname <- paste0("trans-and-zscore_debug_", basename, "_.txt")
# Open the debug/log file
debugfile <- file(paste(debugpath, debugname, sep=""), "w")

# Set various path and filenames
sitecsvpath <- paste(datapath,"rpd_influx/",sep="")
transcsvpath <- paste(datapath,"z_trans_rpd/",sep="") # 
statscsvpath <- paste(datapath,"stats_rpd/",sep="")
statsfile <- paste(statscsvpath,"_",basename,"_z_stats.csv", sep="") #JUST 'stats' for the first run for mnx and box cox.


# Read list of sites
ptm <- proc.time()
sites <- sitelist
ent_vec <- unique(sites$ID_ENTITY)
nsites <- length(ent_vec)
print (nsites)


##Main loop
#Loop over the individual sites, doing the following:
  # compose the site .csv file name
  # read the input data
  # discard samples with missing ages
  # discard samples with ages after 2020 CE
  # initial max rescaling
  # set alpha the Box-Cox transformation shift parameter
  # maximum likelihood estimation of lambda
  # Box-Cox transformation of data
  # max rescaling
  # calculate mean and standard deviation of data over base period
  # calculate z-scores
  # write out transformed data for each site


#####################################
#####################################

# Main transforming loop 

# Storage for statistics (vector arrays to store lambda, likelihoods etc)
ent_vec[1000]
sn_save <- rep(0,nsites)
lam_save <- rep(0,nsites)
lik_save <- rep(0,nsites)
tmean_save <- rep(0,nsites)
tstdev_save <- rep(0,nsites)


for (j in 1:nsites) {  #Loop through each .csv file of each record
  
  #1. Set input files
  ent_n<- ent_vec[j]
  ent_nm <- sites[sites$ID_ENTITY == ent_n,"entity_name"]
  entnmsub <- substr(ent_nm, start = 1, stop = 3)
  
  entidchar <- as.character(ent_n)
  if (ent_n >= 1) entid <- paste("000", entidchar, entnmsub, sep="") 
  if (ent_n >= 10) entid <- paste("00", entidchar, entnmsub, sep="")
  if (ent_n >= 100) entid <- paste("0", entidchar, entnmsub, sep="")
  if (ent_n >= 1000) entid <- paste(    entidchar, entnmsub, sep="")
  
  inputfile <- paste(sitecsvpath, entid, "_RPD", "_influx.csv", sep="")
  print(j);print(entid); print(inputfile) #Printing names etc. to check
  
  # 2. Read the input files
  sitedata <- read.csv(inputfile)
  nsamp <- length(sitedata[,"ID_SAMPLE"])
  nsampchar <- as.character(nsamp)
  writeLines(paste("Site",nsampchar,"samples", sep=" "), con = debugfile, sep = "/n") 
  
  # 3. Discard samples with missing (-9999) ages
  sitedata <- sitedata[sitedata$est_age != -9999,]
  
  # 4. Discard samples with ages > -70
  sitedata <- sitedata[sitedata$est_age > -70,]
  
  if(nrow(sitedata) >1){ #To account for those records that dont have any ages
    
    # 5. Initial minimax / max rescaling of data
    # minimax <- (sitedata$influx-min(sitedata$influx))/(max(sitedata$influx)-min(sitedata$influx)) #Choice of minmax vs max rescaling
    max_re <- sitedata$influx/max(sitedata$influx)
    
    # 6. Set `alpha`, the Box-Cox transformation shift parameter
    alpha <- 0.01  # Box-Cox shift parameter

    # 7. Maximum likelihood estimation of lambda
    # derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1
    
    npts <- 201 # number of estimates of lambda
    # y <- minimax + alpha #Choice of minmax vs max rescaling
    y <- max_re + alpha
    n <- length(y)
    logy <- log(y)
    ydot <- exp(mean(logy))
    lasave <- matrix(1:npts)
    liksave <- matrix(1:npts)
    for (i in 1:npts) {
      la <- -2.0+(i-1)*(4/(npts-1))
      if (la != 0.0) yt <- (y^la-1)/la else yt <- logy*(1+(la*logy)/2*(1+(la*logy)/3*(1+(la*logy)/4)))
      zt <- yt/ydot^(la-1)
      loglik <- -n/2*log(sum((zt - mean(zt))^2 ))
      lasave[i] <- la
      liksave[i] <- loglik
    }
    
    # Save the maximum likelihood value and the associated lambda
    maxlh <- liksave[which.max(liksave)]
    lafit <- lasave[which.max(liksave)]
    print (c(entid, maxlh, lafit))
    
    # 8. Box-Cox transformation of data
    if (lafit == 0.0) tall <- log(y) else tall <- (y^lafit - 1)/lafit 
    
    # 9. Minimax / max rescaling
    # tall <- (tall - min(tall))/(max(tall)-min(tall)) #Again, choice of minimax or max transformation
    tall <- tall/max(tall)
    
    # 10. Calculate mean and standard deviation of data over base period
    tmean <- mean(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
    tstdev <- sd(tall[sitedata$est_age >= basebeg & sitedata$est_age <= baseend])
    
    # 11. Calculate z-scores
    ztrans <- (tall-tmean)/tstdev
    
    # 12. Write out transformed data for each record
    siteout <- data.frame(cbind(sitedata$ID_SITE, sitedata$ID_ENTITY, sitedata$ID_SAMPLE, sitedata$est_age, 
                                sitedata$depth, sitedata$influx, max_re, tall, ztrans))
    colnames(siteout) <- c("ID_SITE", "ID_ENTITY", "ID_SAMPLE", "EST_AGE", "DEPTH", "INFLUX", "influxmnx", "tall", "zt")
    
    outputfile <- paste(transcsvpath, entid, "_Ztrans_RPD_", basename, ".csv", sep="") 
    write.table(siteout, outputfile, col.names=TRUE, row.names=FALSE, sep=",")
    transcsvpath
    sn_save[j] <- entid
    lam_save[j] <- lafit
    lik_save[j] <- maxlh
    tmean_save[j] <- tmean
    tstdev_save[j] <- tstdev
    
  }
}


# Write out a file of statistics
stats <- data.frame(cbind(sn_save, lam_save, lik_save, tmean_save, tstdev_save))
colnames(stats) <- c("site", "lambda", "likelihood", "mean", "stdev")
write.table(stats, statsfile, col.names=TRUE, row.names=FALSE, sep=",")

proc.time() - ptm


# Now make a single file with all the transformed data
dat <- data.frame()
for (y in 1:length(ent_vec)) { 
  ent_n<- ent_vec[y]
  ent_nm <- sites[sites$ID_ENTITY == ent_n,"entity_name"]
  entnmsub <- substr(ent_nm, start = 1, stop = 3)
  
  entidchar <- as.character(ent_n)
  if (ent_n >= 1) entid <- paste("000", entidchar, entnmsub, sep="") 
  if (ent_n >= 10) entid <- paste("00", entidchar, entnmsub, sep="")
  if (ent_n >= 100) entid <- paste("0", entidchar, entnmsub, sep="")
  if (ent_n >= 1000) entid <- paste(    entidchar, entnmsub, sep="")
  
  inputfile <- paste(transcsvpath, entid, "_Ztrans_RPD_zt200_8200k.csv", sep="")
  print(y);print(entid); print(inputfile) 
  if (file.exists(inputfile)) sitedata <- read.csv(inputfile)
  dat <- rbind(dat, sitedata)
}

length(unique(dat$ID_ENTITY))
datapath
z_scores <- dat %>% 
  dplyr::left_join(dplyr::select(rpd_data_input, ID_ENTITY, entity_name), by = "ID_ENTITY")
write.csv(dat, paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = ''), row.names = F)

## end. Now you have a dataframe with box cox and z score transformed data.
## Note the following composite curve generation is replicated for the neolithic and SEA analysis 


# ---------------------------------------------------------







# 4. Pre-bin data (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

datapath <- "./other_output/rpd/"
sitelist <- read.csv("./other_output/rpd/rpd_influx/RPD_infl_mtdata.csv") 
debugpath <- "./other_output/rpd/rpd_debug/" 
debugname <- "curve_fit_lat_debug_zt200_8200k.txt"
debugfile <- file(paste(debugpath, debugname, sep=""), "w")

# Set baseperiod - here you will use this to select sites based on the degree of overlap with the baseperiod, if you want to
basebeg <- 200
baseend <- 8200 
basename <- paste0("zt", basebeg, "_", baseend, "k")

dat <- read.csv(paste(datapath, 'RPD_trans_zt200_8200k.csv', sep = ''))
length(unique(dat$ID_ENTITY)) 
ent_vec2 <- unique(dat$ID_ENTITY)

lsover <- ent_vec2 #This line of code enables subsetting of dataset if required (not utilised here)


# Binning Z scores
Zbinning <- function(dt, bins, binhw){
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
  BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix") #originally class = 'matrix'
  return(BinnedData)
}

# Run multiple bin widths to allow for later composite generation
bw <- 100 #Set as standard for composite
bw_50 <- 50 #Alternative bin widths
bw_200 <- 200 #Alternative bin widths
bw_500 <- 500 #Alternative bin widths

Zbin_dat <- Zbinning(dt = dat, bins = seq(0,12000,bw), binhw = bw/2)
Zbin_dat_50 <- Zbinning(dt = dat, bins = seq(0,12000,bw_50), binhw = bw_50/2)
Zbin_dat_200 <- Zbinning(dt = dat, bins = seq(0,12000,bw_200), binhw = bw_200/2)
Zbin_dat_500 <- Zbinning(dt = dat, bins = seq(0,12000,bw_500), binhw = bw_500/2)

#Give ages of each bin as a column
bin_age <- seq(0,12000,bw)
bin_age_50 <- seq(0,12000,bw_50)
bin_age_200 <- seq(0,12000,bw_200)
bin_age_500 <- seq(0,12000,bw_500)

Zbin_dat <- cbind(bin_age, Zbin_dat)
Zbin_dat_50 <- cbind(bin_age_50, Zbin_dat_50)
Zbin_dat_200 <- cbind(bin_age_200, Zbin_dat_200)
Zbin_dat_500 <- cbind(bin_age_500, Zbin_dat_500)

## Limit this to the sites in lsover. 
lsoverch <- as.character(lsover)
Zbin_dat <- Zbin_dat[,c('bin_age',lsoverch)]
Zbin_dat_50 <- Zbin_dat_50[,c('bin_age_50',lsoverch)]
Zbin_dat_200 <- Zbin_dat_200[,c('bin_age_200',lsoverch)]
Zbin_dat_500 <- Zbin_dat_500[,c('bin_age_500',lsoverch)]
check <- colnames(Zbin_dat)
check[!(check %in% lsoverch)]

# Dataframe for curve generation
Zbin_dat <- as.data.frame(Zbin_dat) 
Zbin_dat_50 <- as.data.frame(Zbin_dat_50) 
Zbin_dat_200 <- as.data.frame(Zbin_dat_200) 
Zbin_dat_500 <- as.data.frame(Zbin_dat_500) 

site_bin_info <- Zbin_dat %>% 
  dplyr::as_tibble() %>% 
  dplyr::filter(dplyr::between(bin_age, 3500, 10000)) %>% 
  dplyr::na_if("NaN") 
site_bin_info_n <- colSums(!is.na(site_bin_info)) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(ID_ENTITY = colnames(site_bin_info), .before = 1) %>%
  dplyr::filter(ID_ENTITY != "bin_age") %>% 
  dplyr::mutate(ID_ENTITY = as.numeric(ID_ENTITY)) %>% 
  dplyr::rename(bin_number = value)

# ---------------------------------------------------------









# 5. Composite curves (adapted from Bartlein - see https://pjbartlein.github.io/GCDv3Analysis/)
# ---------------------------------------------------------
# ---------------------------------------------------------

## Set up

# Set names
basename <- paste0("zt", basebeg, "_", baseend, "k")
binname <- paste0("bw", bw)

# Paths for input and output .csv files 
datapath <- "./other_output/rpd/"

# Locfit (half) window-width parameter. Different values set for later composite curve generation
hw <- 500 # bandwidth (smoothing parameter)
hw_250 <- 250 #Alternative bandwidth
hw_1000 <- 1000 #Alternative bandwidth

# Number of bootstrap samples/replications
nreps <- 1000

# Target ages for fitted values
targbeg <- 0
targend <- 12000
targstep <- bw ## this is based on the previous binning procedure

targbeg_bw_50 <- 0
targend_bw_50 <- 12000
targstep_bw_50 <- bw_50 ## this is based on the previous binning procedure

targbeg_bw_200 <- 0
targend_bw_200 <- 12000
targstep_bw_200 <- bw_200 ## this is based on the previous binning procedure

targbeg_bw_500 <- 0
targend_bw_500 <- 12000
targstep_bw_500 <- bw_500 ## this is based on the previous binning procedure


# Array sizes
maxrecs <- 2000
maxreps <- 1000

# Plot output 
plotout <- "pdf"



## A. Individual curve plotting

# Curve output path and file
curvecsvpath <- "./figs/rpd/"
curvefile <- paste(curvecsvpath,"rpd",".csv", sep="") #adjust to particular data you use now. .
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { 
  pdffile <- paste(curvecsvpath, "rpd", ".pdf", sep="")
  print(pdffile)
}

# Read the list of sites to be processed
ns <- as.numeric(colnames(Zbin_dat))
ns <- na.omit(ns)

# Arrays for data and fitted values
age <- matrix(NA, ncol=length(ns), nrow=maxrecs)
influx <- matrix(NA, ncol=length(ns), nrow=maxrecs)
nsamples <- rep(0, maxrecs)
targage <- seq(targbeg,targend,targstep)
targage.df <- data.frame(x=targage)
lowage <- targage - hw; highage <- targage + hw
ntarg <- length(targage)
yfit <- matrix(NA, nrow=length(targage.df$x), ncol=maxreps)

# Arrays for sample number and effective window span tracking
ndec <- matrix(0, ncol=ntarg, nrow=length(ns))
ndec_tot <- rep(0, ntarg)
xspan <- rep(0, ntarg)
ninwin <- matrix(0, ncol=ntarg, nrow=length(ns))
ninwin_tot <- rep(0, ntarg)


## Read data

# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0

for(i in 2:length(Zbin_dat)){ 
  indata <- Zbin_dat[,c(i,1)] 
  indata <- na.omit(indata) 
  nsamp <-  length(indata$bin_age) 
  if (nsamp > 0) {
    ii <- ii+1
    age[1:nsamp,ii] <- indata$bin_age 
    influx[1:nsamp,ii] <- indata[,1]  
    nsamples[ii] <- nsamp
  }
  print(i)
}

nsites <- ii; nsites #Number of sites with data

# Trim samples to age range
influx[age >= targend+hw] <- NA
age[age >= targend+hw] <- NA

# Censor abs(influx) values > 10
influx[abs(influx) >= 10] <- NA
age[abs(influx) >= 10] <- NA  


##Find number of sites and samples contributing to fitted values

# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:ntarg) {
  agemax <- -1e32; agemin <- 1e32
  for (j in 1:nsites) {
    for (k in 1:nsamples[j]) {
      if (!is.na(age[k,j])) {
        ii <- (age[k,j]-targage[1])/targstep + 1
        if (ii > 0 && ii <= ntarg) {ndec[j,ii] = 1}
        if (age[k,j] >= targage[i]-hw && age[k,j] <= targage[i]+hw) {
          ninwin[j,i] = 1
          if (agemax < age[k,j]) {agemax <- age[k,j]}
          if (agemin > age[k,j]) {agemin <- age[k,j]}
        }
      }
    }
  }
  ndec_tot[i] <- sum(ndec[,i])
  ninwin_tot[i] <- sum(ninwin[,i])
  xspan[i] <- agemax - agemin
}
proc.time() - ptm



## Curve-fitting and bootstrapping

# Composite curve
ptm <- proc.time()

# Reshape matrices into vectors 
x <- as.vector(age)
y <- as.vector(influx)
lfdata <- data.frame(x,y)
lfdata <- na.omit(lfdata)
x <- lfdata$x; y <- lfdata$y

lfdata_corr <- lfdata #Create dataframe for later correlation analysis

# Locfit
loc01 <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01)

# Get fitted values
pred01 <- predict(loc01, newdata=targage.df, se.fit=TRUE)
loc01_fit <- data.frame(targage.df$x, pred01$fit)
fitname <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
head(loc01_fit)

proc.time() - ptm


#Bootstrap-by-site confidence intervals
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications
curvecsvpath
pdffile

if (plotout == "pdf") {pdf(file=pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(x, y, xlab="Years cal. BP", ylab="Z-scores of transformed charcoal influx", xlim=c(10000,3500), ylim=c(-3,1), type="n", cex.lab = 1.2, mgp = c(2.4, 1, 0))
grid()

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) {
  print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  locboot <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
  predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
  yfit[,i] <- predboot$fit
  lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}


# Step 3 -- Plot the unresampled (initial) fit
fitname <- paste("Zlocfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="red")

# Step 4 -- Find and add bootstrap CIs
yfit95 <- apply(yfit, 1, function(x) quantile(x,prob=0.975, na.rm=T))
yfit05 <- apply(yfit, 1, function(x) quantile(x,prob=0.025, na.rm=T))
lines(targage.df$x, yfit95, lwd=1, col="red")
lines(targage.df$x, yfit05, lwd=1, col="red")

if (plotout == "pdf") {dev.off()}


# Output - saving the output enables construction of differing curves for comparison
curveout <- data.frame(cbind(targage.df$x, pred01$fit, yfit95, yfit05, ndec_tot, xspan, ninwin_tot))
colnames(curveout) <- c("age", "locfit", "cu95", "cl95", "nsites", "window", "ninwin")
write.csv(curveout, curvefile, row.names = F)

proc.time() - ptm
#End individual curve plotting

#################################



# B. Plot different smoothing windows (as above code, but with multiple smoothing windows on the same plot)

# Curve output path and file
curvecsvpath <- "./figs/rpd/"
curvefile <- paste(curvecsvpath,"rpd_sm",".csv", sep="") #adjust to particular data you use now. .
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") {
  pdffile <- paste(curvecsvpath,"rpd_sm", ".pdf", sep="")
  print(pdffile)
}

#Read the list of sites to be processed
ns_hw_250 <- as.numeric(colnames(Zbin_dat))
ns_hw_250 <- na.omit(ns_hw_250) ##Z score data. note that this excludes the age column and only takes entity ids. so be careful when using the length of this in loop below.

ns_hw_1000 <- as.numeric(colnames(Zbin_dat))
ns_hw_1000 <- na.omit(ns_hw_1000) ##Z score data. note that this excludes the age column and only takes entity ids. so be careful when using the length of this in loop below.


# Arrays for data and fitted values
age_hw_250 <- matrix(NA, ncol=length(ns_hw_250), nrow=maxrecs)
influx_hw_250 <- matrix(NA, ncol=length(ns_hw_250), nrow=maxrecs)
nsamples_hw_250 <- rep(0, maxrecs)
targage_hw_250 <- seq(targbeg,targend,targstep) 
targage.df_hw_250 <- data.frame(x_hw_250=targage_hw_250)
lowage_hw_250 <- targage_hw_250 - hw_250; highage_hw_250 <- targage_hw_250 + hw_250
ntarg_hw_250 <- length(targage_hw_250)
yfit_hw_250 <- matrix(NA, nrow=length(targage.df_hw_250$x_hw_250), ncol=maxreps) #the actual curve values (fitted)

age_hw_1000 <- matrix(NA, ncol=length(ns_hw_1000), nrow=maxrecs)
influx_hw_1000 <- matrix(NA, ncol=length(ns_hw_1000), nrow=maxrecs)
nsamples_hw_1000 <- rep(0, maxrecs)
targage_hw_1000 <- seq(targbeg,targend,targstep) 
targage.df_hw_1000 <- data.frame(x_hw_1000=targage_hw_1000)
lowage_hw_1000 <- targage_hw_1000 - hw_1000; highage_hw_1000 <- targage_hw_1000 + hw_1000
ntarg_hw_1000 <- length(targage_hw_1000)
yfit_hw_1000 <- matrix(NA, nrow=length(targage.df_hw_1000$x_hw_1000), ncol=maxreps) #the actual curve values (fitted)

# Arrays for sample number and effective window span tracking
ndec_hw_250 <- matrix(0, ncol=ntarg_hw_250, nrow=length(ns_hw_250))
ndec_tot_hw_250 <- rep(0, ntarg_hw_250)
xspan_hw_250 <- rep(0, ntarg_hw_250)
ninwin_hw_250 <- matrix(0, ncol=ntarg_hw_250, nrow=length(ns_hw_250))
ninwin_tot_hw_250 <- rep(0, ntarg_hw_250)

ndec_hw_1000 <- matrix(0, ncol=ntarg_hw_1000, nrow=length(ns_hw_1000))
ndec_tot_hw_1000 <- rep(0, ntarg_hw_1000)
xspan_hw_1000 <- rep(0, ntarg_hw_1000)
ninwin_hw_1000 <- matrix(0, ncol=ntarg_hw_1000, nrow=length(ns_hw_1000))
ninwin_tot_hw_1000 <- rep(0, ntarg_hw_1000)


## Read data

# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0
for(i in 2:length(Zbin_dat)){ 
  indata_hw_250 <- Zbin_dat[,c(i,1)] 
  indata_hw_250 <- na.omit(indata_hw_250) 
  nsamp_hw_250 <-  length(indata_hw_250$bin_age) 
  if (nsamp_hw_250 > 0) {
    ii <- ii+1
    age_hw_250[1:nsamp_hw_250,ii] <- indata_hw_250$bin_age 
    influx_hw_250[1:nsamp_hw_250,ii] <- indata_hw_250[,1]  
    nsamples_hw_250[ii] <- nsamp_hw_250
  }
  print(i)
}

nsites_hw_250 <- ii; nsites_hw_250 #Number of sites with data


ii <- 0
for(i in 2:length(Zbin_dat)){ 
  indata_hw_1000 <- Zbin_dat[,c(i,1)] 
  indata_hw_1000 <- na.omit(indata_hw_1000) 
  nsamp_hw_1000 <-  length(indata_hw_1000$bin_age) 
  if (nsamp_hw_1000 > 0) {
    ii <- ii+1
    age_hw_1000[1:nsamp_hw_1000,ii] <- indata_hw_1000$bin_age
    influx_hw_1000[1:nsamp_hw_1000,ii] <- indata_hw_1000[,1]  
    nsamples_hw_1000[ii] <- nsamp_hw_1000
  }
  print(i)
}

nsites_hw_1000 <- ii; nsites_hw_1000 #Number of sites with data


# Trim samples to age range
influx_hw_250[age_hw_250 >= targend+hw_250] <- NA
age_hw_250[age_hw_250 >= targend+hw_250] <- NA 

influx_hw_1000[age_hw_1000 >= targend+hw_1000] <- NA
age_hw_1000[age_hw_1000 >= targend+hw_1000] <- NA 

# Censor abs(influx) values > 10
influx_hw_250[abs(influx_hw_250) >= 10] <- NA
age_hw_250[abs(influx_hw_250) >= 10] <- NA 

influx_hw_1000[abs(influx_hw_1000) >= 10] <- NA
age_hw_1000[abs(influx_hw_1000) >= 10] <- NA 


## Find number of sites and samples contributing to fitted values

# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:ntarg_hw_250) {
  agemax_hw_250 <- -1e32; agemin_hw_250 <- 1e32
  for (j in 1:nsites_hw_250) {
    for (k in 1:nsamples_hw_250[j]) {
      if (!is.na(age_hw_250[k,j])) {
        ii <- (age_hw_250[k,j]-targage_hw_250[1])/targstep + 1
        if (ii > 0 && ii <= ntarg_hw_250) {ndec_hw_250[j,ii] = 1}
        if (age_hw_250[k,j] >= targage_hw_250[i]-hw_250 && age_hw_250[k,j] <= targage_hw_250[i]+hw_250) {
          ninwin_hw_250[j,i] = 1
          if (agemax_hw_250 < age_hw_250[k,j]) {agemax_hw_250 <- age_hw_250[k,j]}
          if (agemin_hw_250 > age_hw_250[k,j]) {agemin_hw_250 <- age_hw_250[k,j]}
        }
      }
    }
  }
  ndec_tot_hw_250[i] <- sum(ndec_hw_250[,i])
  ninwin_tot_hw_250[i] <- sum(ninwin_hw_250[,i])
  xspan_hw_250[i] <- agemax_hw_250 - agemin_hw_250
}
proc.time() - ptm


ptm <- proc.time()
for (i in 1:ntarg_hw_1000) {
  agemax_hw_1000 <- -1e32; agemin_hw_1000 <- 1e32
  for (j in 1:nsites_hw_1000) {
    for (k in 1:nsamples_hw_1000[j]) {
      if (!is.na(age_hw_1000[k,j])) {
        ii <- (age_hw_1000[k,j]-targage_hw_1000[1])/targstep + 1
        if (ii > 0 && ii <= ntarg_hw_1000) {ndec_hw_1000[j,ii] = 1}
        if (age_hw_1000[k,j] >= targage_hw_1000[i]-hw_1000 && age_hw_1000[k,j] <= targage_hw_1000[i]+hw_1000) {
          ninwin_hw_1000[j,i] = 1
          if (agemax_hw_1000 < age_hw_1000[k,j]) {agemax_hw_1000 <- age_hw_1000[k,j]}
          if (agemin_hw_1000 > age_hw_1000[k,j]) {agemin_hw_1000 <- age_hw_1000[k,j]}
        }
      }
    }
  }
  ndec_tot_hw_1000[i] <- sum(ndec_hw_1000[,i])
  ninwin_tot_hw_1000[i] <- sum(ninwin_hw_1000[,i])
  xspan_hw_1000[i] <- agemax_hw_1000 - agemin_hw_1000
}
proc.time() - ptm



## Curve-fitting and bootstrapping

#Composite curve
ptm <- proc.time()

# Reshape matrices into vectors 
x_hw_250 <- as.vector(age_hw_250)
y_hw_250 <- as.vector(influx_hw_250)
lfdata_hw_250 <- data.frame(x_hw_250,y_hw_250)
lfdata_hw_250 <- na.omit(lfdata_hw_250)
x_hw_250 <- lfdata_hw_250$x_hw_250; y_hw_250 <- lfdata_hw_250$y_hw_250

x_hw_1000 <- as.vector(age_hw_1000)
y_hw_1000 <- as.vector(influx_hw_1000)
lfdata_hw_1000 <- data.frame(x_hw_1000,y_hw_1000)
lfdata_hw_1000 <- na.omit(lfdata_hw_1000)
x_hw_1000 <- lfdata_hw_1000$x_hw_1000; y_hw_1000 <- lfdata_hw_1000$y_hw_1000

# Locfit
loc01_hw_250 <- locfit::locfit(y_hw_250 ~ locfit::lp(x_hw_250, deg=1, h=hw_250), maxk=800, family="qrgauss")
summary(loc01_hw_250)

loc01_hw_1000 <- locfit::locfit(y_hw_1000 ~ locfit::lp(x_hw_1000, deg=1, h=hw_1000), maxk=800, family="qrgauss")
summary(loc01_hw_1000)

# Get fitted values
pred01_hw_250 <- predict(loc01_hw_250, newdata=targage.df_hw_250, se.fit=TRUE)
loc01_fit_hw_250 <- data.frame(targage.df_hw_250$x_hw_250, pred01_hw_250$fit)
fitname_hw_250 <- paste("locfit_",as.character(hw_250), sep="")
colnames(loc01_fit_hw_250) <- c("age_hw_250", fitname_hw_250)
head(loc01_fit_hw_250)

proc.time() - ptm

pred01_hw_1000 <- predict(loc01_hw_1000, newdata=targage.df_hw_1000, se.fit=TRUE)
loc01_fit_hw_1000 <- data.frame(targage.df_hw_1000$x_hw_1000, pred01_hw_1000$fit)
fitname_hw_1000 <- paste("locfit_",as.character(hw_1000), sep="")
colnames(loc01_fit_hw_1000) <- c("age_hw_1000", fitname_hw_1000)
head(loc01_fit_hw_1000)

proc.time() - ptm

#Bootstrap-by-site confidence intervals (as before)
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications
curvecsvpath
pdffile

if (plotout == "pdf") {pdf(file=pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(x, y, xlab="Years cal. BP", ylab="Z-scores of transformed charcoal influx", xlim=c(10000,3500), ylim=c(-3,1), type="n", cex.lab = 1.2, mgp = c(2.4, 1, 0)) ##make sure xlim is right.
grid()
legend(8000, -1.5, legend = c("250-year", "500-year", "1000-year"), col = c("blue", "red", "black"), title="Smoothing window (half) width", lty=1, lwd = 2, cex=1, box.lty=0)

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) {
  print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  locboot <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
  predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
  yfit[,i] <- predboot$fit
  lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}

# Step 3 -- Plot the unresampled (initial) fits

#hw_250
fitname_hw_250 <- paste("Zlocfit_",as.character(hw_250), sep="")
colnames(loc01_fit_hw_250) <- c("age_hw_250", fitname_hw_250)
lines(loc01_fit_hw_250[,1], loc01_fit_hw_250[,2], lwd=2, col="blue")

#hw_1000
fitname_hw_1000 <- paste("Zlocfit_",as.character(hw_1000), sep="")
colnames(loc01_fit_hw_1000) <- c("age_hw_1000", fitname_hw_1000)
lines(loc01_fit_hw_1000[,1], loc01_fit_hw_1000[,2], lwd=2, col="black")

#Normal_hw
fitname <- paste("Zlocfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="red") 


# Step 4 -- Find and add bootstrap CIs from the normal curve
lines(targage.df$x, yfit95, lwd=1, col="red")
lines(targage.df$x, yfit05, lwd=1, col="red")



if (plotout == "pdf") {dev.off()}
proc.time() - ptm
#End multiple smoothing bandwith plotting




#################################

# C. Plot different binwidths


# Curve output path and file
curvecsvpath <- "./figs/rpd/"
curvefile <- paste(curvecsvpath,"rpd_bw",".csv", sep="") #adjust to particular data you use now. .
print(curvefile)

# .pdf plot of bootstrap iterations
if (plotout == "pdf") { 
  pdffile <- paste(curvecsvpath,"rpd_bw", ".pdf", sep="")
  print(pdffile)
}

#Read the list of sites to be processed
ns_bw_50 <- as.numeric(colnames(Zbin_dat_50))
ns_bw_50 <- na.omit(ns_bw_50) 
length(ns_bw_50)

ns_bw_200 <- as.numeric(colnames(Zbin_dat_200))
ns_bw_200 <- na.omit(ns_bw_200) 
length(ns_bw_200)

ns_bw_500 <- as.numeric(colnames(Zbin_dat_500))
ns_bw_500 <- na.omit(ns_bw_500) 
length(ns_bw_500)



# Arrays for data and fitted values
age_bw_50 <- matrix(NA, ncol=length(ns_bw_50), nrow=maxrecs)
influx_bw_50 <- matrix(NA, ncol=length(ns_bw_50), nrow=maxrecs)
nsamples_bw_50<- rep(0, maxrecs)
targage_bw_50 <- seq(targbeg_bw_50,targend_bw_50,targstep_bw_50) 
targage.df_bw_50 <- data.frame(x_bw_50=targage_bw_50)
lowage_bw_50 <- targage_bw_50 - bw_50; highage_bw_50 <- targage_bw_50 + bw_50
ntarg_bw_50 <- length(targage_bw_50)
yfit_bw_50 <- matrix(NA, nrow=length(targage.df_bw_50$x_bw_50), ncol=maxreps) 

age_bw_200 <- matrix(NA, ncol=length(ns_bw_200), nrow=maxrecs)
influx_bw_200 <- matrix(NA, ncol=length(ns_bw_200), nrow=maxrecs)
nsamples_bw_200 <- rep(0, maxrecs)
targage_bw_200 <- seq(targbeg_bw_200,targend_bw_200,targstep_bw_200) 
targage.df_bw_200 <- data.frame(x_bw_200=targage_bw_200)
lowage_bw_200 <- targage_bw_200 - bw_200; highage_bw_200 <- targage_bw_200 + bw_200
ntarg_bw_200 <- length(targage_bw_200)
yfit_bw_200 <- matrix(NA, nrow=length(targage.df_bw_200$x_bw_200), ncol=maxreps) 

age_bw_500 <- matrix(NA, ncol=length(ns_bw_500), nrow=maxrecs)
influx_bw_500 <- matrix(NA, ncol=length(ns_bw_500), nrow=maxrecs)
nsamples_bw_500 <- rep(0, maxrecs)
targage_bw_500 <- seq(targbeg_bw_500,targend_bw_500,targstep_bw_500) 
targage.df_bw_500 <- data.frame(x_bw_500=targage_bw_500)
lowage_bw_500 <- targage_bw_500 - bw_500; highage_bw_500 <- targage_bw_500 + bw_500
ntarg_bw_500 <- length(targage_bw_500)
yfit_bw_500 <- matrix(NA, nrow=length(targage.df_bw_500$x_bw_500), ncol=maxreps) 

# Arrays for sample number and effective window span tracking
ndec_bw_50 <- matrix(0, ncol=ntarg_bw_50, nrow=length(ns_bw_50))
ndec_tot_bw_50 <- rep(0, ntarg_bw_50)
xspan_bw_50 <- rep(0, ntarg_bw_50)
ninwin_bw_50 <- matrix(0, ncol=ntarg_bw_50, nrow=length(ns_bw_50))
ninwin_tot_bw_50 <- rep(0, ntarg_bw_50)

ndec_bw_200 <- matrix(0, ncol=ntarg_bw_200, nrow=length(ns_bw_200))
ndec_tot_bw_200 <- rep(0, ntarg_bw_200)
xspan_bw_200 <- rep(0, ntarg_bw_200)
ninwin_bw_200 <- matrix(0, ncol=ntarg_bw_200, nrow=length(ns_bw_200))
ninwin_tot_bw_200 <- rep(0, ntarg_bw_200)

ndec_bw_500 <- matrix(0, ncol=ntarg_bw_500, nrow=length(ns_bw_500))
ndec_tot_bw_500 <- rep(0, ntarg_bw_500)
xspan_bw_500 <- rep(0, ntarg_bw_500)
ninwin_bw_500 <- matrix(0, ncol=ntarg_bw_500, nrow=length(ns_bw_500))
ninwin_tot_bw_500 <- rep(0, ntarg_bw_500)


## Read data

# Read and store the presample (binned) files as matrices of ages and influx values
ii <- 0
for(i in 2:length(Zbin_dat_50)){ 
  indata_bw_50 <- Zbin_dat_50[,c(i,1)] 
  indata_bw_50 <- na.omit(indata_bw_50) 
  nsamp_bw_50 <-  length(indata_bw_50$bin_age) 
  if (nsamp_bw_50 > 0) {
    ii <- ii+1
    age_bw_50[1:nsamp_bw_50,ii] <- indata_bw_50$bin_age 
    influx_bw_50[1:nsamp_bw_50,ii] <- indata_bw_50[,1]  
    nsamples_bw_50[ii] <- nsamp_bw_50
  }
  print(i)
}

nsites_bw_50 <- ii; nsites_bw_50 #Number of sites with data

ii <- 0
for(i in 2:length(Zbin_dat_200)){ 
  indata_bw_200 <- Zbin_dat_200[,c(i,1)] 
  indata_bw_200 <- na.omit(indata_bw_200)  
  nsamp_bw_200 <-  length(indata_bw_200$bin_age) 
  if (nsamp_bw_200 > 0) {
    ii <- ii+1
    age_bw_200[1:nsamp_bw_200,ii] <- indata_bw_200$bin_age 
    influx_bw_200[1:nsamp_bw_200,ii] <- indata_bw_200[,1] 
    nsamples_bw_200[ii] <- nsamp_bw_200
  }
  print(i)
}

nsites_bw_200 <- ii; nsites_bw_200 #Number of sites with data


ii <- 0
for(i in 2:length(Zbin_dat_500)){ 
  indata_bw_500 <- Zbin_dat_500[,c(i,1)] 
  indata_bw_500 <- na.omit(indata_bw_500) 
  nsamp_bw_500 <-  length(indata_bw_500$bin_age) 
  if (nsamp_bw_500 > 0) {
    ii <- ii+1
    age_bw_500[1:nsamp_bw_500,ii] <- indata_bw_500$bin_age 
    influx_bw_500[1:nsamp_bw_500,ii] <- indata_bw_500[,1]  
    nsamples_bw_500[ii] <- nsamp_bw_500
  }
  print(i)
}

nsites_bw_500 <- ii; nsites_bw_500 #Number of sites with data


# Trim samples to age range
influx_bw_50[age_bw_50 >= targend_bw_50+bw_50] <- NA
age_bw_50[age_bw_50 >= targend_bw_50+bw_50] <- NA 

influx_bw_200[age_bw_200 >= targend_bw_200+bw_200] <- NA
age_bw_200[age_bw_200 >= targend_bw_200+bw_200] <- NA 

influx_bw_500[age_bw_500 >= targend_bw_500+bw_500] <- NA
age_bw_50[age_bw_500 >= targend_bw_500+bw_500] <- NA 


# Censor abs(influx) values > 10
influx_bw_50[abs(influx_bw_50) >= 10] <- NA
age_bw_50[abs(influx_bw_50) >= 10] <- NA 

influx_bw_200[abs(influx_bw_200) >= 10] <- NA
age_bw_200[abs(influx_bw_200) >= 10] <- NA 

influx_bw_500[abs(influx_bw_500) >= 10] <- NA
age_bw_500[abs(influx_bw_500) >= 10] <- NA 


## Find number of sites and samples contributing to fitted values


# Count number of sites that contributed to each fitted value
ptm <- proc.time()
for (i in 1:ntarg_bw_50) {
  agemax_bw_50 <- -1e32; agemin_bw_50 <- 1e32
  for (j in 1:nsites_bw_50) {
    for (k in 1:nsamples_bw_50[j]) {
      if (!is.na(age_bw_50[k,j])) {
        ii <- (age_bw_50[k,j]-targage_bw_50[1])/targstep + 1
        if (ii > 0 && ii <= ntarg_bw_50) {ndec_bw_50[j,ii] = 1}
        if (age_bw_50[k,j] >= targage_bw_50[i]-bw_50 && age_bw_50[k,j] <= targage_bw_50[i]+bw_50) {
          ninwin_bw_50[j,i] = 1
          if (agemax_bw_50 < age_bw_50[k,j]) {agemax_bw_50 <- age_bw_50[k,j]}
          if (agemin_bw_50 > age_bw_50[k,j]) {agemin_bw_50 <- age_bw_50[k,j]}
        }
      }
    }
  }
  ndec_tot_bw_50[i] <- sum(ndec_bw_50[,i])
  ninwin_tot_bw_50[i] <- sum(ninwin_bw_50[,i])
  xspan_bw_50[i] <- agemax_bw_50 - agemin_bw_50
}
proc.time() - ptm


ptm <- proc.time()
for (i in 1:ntarg_bw_200) {
  agemax_bw_200 <- -1e32; agemin_bw_200 <- 1e32
  for (j in 1:nsites_bw_200) {
    for (k in 1:nsamples_bw_200[j]) {
      if (!is.na(age_bw_200[k,j])) {
        ii <- (age_bw_200[k,j]-targage_bw_200[1])/targstep + 1
        if (ii > 0 && ii <= ntarg_bw_200) {ndec_bw_200[j,ii] = 1}
        if (age_bw_200[k,j] >= targage_bw_200[i]-bw_200 && age_bw_200[k,j] <= targage_bw_200[i]+bw_200) {
          ninwin_bw_200[j,i] = 1
          if (agemax_bw_200 < age_bw_200[k,j]) {agemax_bw_200 <- age_bw_200[k,j]}
          if (agemin_bw_200 > age_bw_200[k,j]) {agemin_bw_200 <- age_bw_200[k,j]}
        }
      }
    }
  }
  ndec_tot_bw_200[i] <- sum(ndec_bw_200[,i])
  ninwin_tot_bw_200[i] <- sum(ninwin_bw_200[,i])
  xspan_bw_200[i] <- agemax_bw_200 - agemin_bw_200
}
proc.time() - ptm


ptm <- proc.time()
for (i in 1:ntarg_bw_500) {
  agemax_bw_500 <- -1e32; agemin_bw_500 <- 1e32
  for (j in 1:nsites_bw_500) {
    for (k in 1:nsamples_bw_500[j]) {
      if (!is.na(age_bw_500[k,j])) {
        ii <- (age_bw_500[k,j]-targage_bw_500[1])/targstep + 1
        if (ii > 0 && ii <= ntarg_bw_500) {ndec_bw_500[j,ii] = 1}
        if (age_bw_500[k,j] >= targage_bw_500[i]-bw_500 && age_bw_500[k,j] <= targage_bw_500[i]+bw_500) {
          ninwin_bw_500[j,i] = 1
          if (agemax_bw_500 < age_bw_500[k,j]) {agemax_bw_500 <- age_bw_500[k,j]}
          if (agemin_bw_500 > age_bw_500[k,j]) {agemin_bw_500 <- age_bw_500[k,j]}
        }
      }
    }
  }
  ndec_tot_bw_500[i] <- sum(ndec_bw_500[,i])
  ninwin_tot_bw_500[i] <- sum(ninwin_bw_500[,i])
  xspan_bw_500[i] <- agemax_bw_500 - agemin_bw_500
}
proc.time() - ptm



## Curve-fitting and bootstrapping

#Composite curve
ptm <- proc.time()

# Reshape matrices into vectors 
x_bw_50 <- as.vector(age_bw_50)
y_bw_50 <- as.vector(influx_bw_50)
lfdata_bw_50 <- data.frame(x_bw_50,y_bw_50)
lfdata_bw_50 <- na.omit(lfdata_bw_50)
x_bw_50 <- lfdata_bw_50$x_bw_50; y_bw_50 <- lfdata_bw_50$y_bw_50

x_bw_200 <- as.vector(age_bw_200)
y_bw_200 <- as.vector(influx_bw_200)
lfdata_bw_200 <- data.frame(x_bw_200,y_bw_200)
lfdata_bw_200 <- na.omit(lfdata_bw_200)
x_bw_200 <- lfdata_bw_200$x_bw_200; y_bw_200 <- lfdata_bw_200$y_bw_200

x_bw_500 <- as.vector(age_bw_500)
y_bw_500 <- as.vector(influx_bw_500)
lfdata_bw_500 <- data.frame(x_bw_500,y_bw_500)
lfdata_bw_500 <- na.omit(lfdata_bw_500)
x_bw_500 <- lfdata_bw_500$x_bw_500; y_bw_500 <- lfdata_bw_500$y_bw_500

# Locfit
loc01_bw_50 <- locfit::locfit(y_bw_50 ~ locfit::lp(x_bw_50, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01_bw_50)

loc01_bw_200 <- locfit::locfit(y_bw_200 ~ locfit::lp(x_bw_200, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01_bw_200)

loc01_bw_500 <- locfit::locfit(y_bw_500 ~ locfit::lp(x_bw_500, deg=1, h=hw), maxk=800, family="qrgauss")
summary(loc01_bw_500)

# Get  fitted values
pred01_bw_50 <- predict(loc01_bw_50, newdata=targage.df_bw_50, se.fit=TRUE)
loc01_fit_bw_50 <- data.frame(targage.df_bw_50$x_bw_50, pred01_bw_50$fit)
fitname_bw_50 <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit_bw_50) <- c("age_bw_50", fitname_bw_50)
head(loc01_fit_bw_50)

proc.time() - ptm

pred01_bw_200 <- predict(loc01_bw_200, newdata=targage.df_bw_200, se.fit=TRUE)
loc01_fit_bw_200 <- data.frame(targage.df_bw_200$x_bw_200, pred01_bw_200$fit)
fitname_bw_200 <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit_bw_200) <- c("age_bw_200", fitname_bw_200)
head(loc01_fit_bw_200)

proc.time() - ptm

pred01_bw_500 <- predict(loc01_bw_500, newdata=targage.df_bw_500, se.fit=TRUE)
loc01_fit_bw_500 <- data.frame(targage.df_bw_500$x_bw_500, pred01_bw_500$fit)
fitname_bw_500 <- paste("locfit_",as.character(hw), sep="")
colnames(loc01_fit_bw_500) <- c("age_bw_500", fitname_bw_500)
head(loc01_fit_bw_500)

proc.time() - ptm


# Bootstrap-by-site confidence intervals

# Bootstrap samples
ptm <- proc.time()

# Step 1 -- Set up to plot individual replications - grey lines I assume
curvecsvpath
pdffile

if (plotout == "pdf") {pdf(file=pdffile, height = 3.15, width = 3.5, pointsize = 9, family = "Times")}
par(mar = c(4, 4, 1, 1))
plot(x, y, xlab="Years cal. BP", ylab="Z-scores of transformed charcoal influx", xlim=c(10000,3500), ylim=c(-3,1), type="n", cex.lab = 1.2, mgp = c(2.4, 1, 0)) ##make sure xlim is right.
grid()
legend(8000, -1.5, legend = c("50-year", "100-year", "200-year", "500-year"), col = c("blue", "red", "black", "green"), title="Binwidth", lty=1, lwd = 2, cex=1, box.lty=0)

# Step 2 -- Do the bootstrap iterations, and plot each composite curve
set.seed(42) # To get the same sequence of random samples for each run

for (i in 1:nreps) { 
  print(i)
  randsitenum <- sample(seq(1:nsites), nsites, replace=TRUE)
  x <- as.vector(age[,randsitenum])
  y <- as.vector(influx[,randsitenum])
  lfdata <- data.frame(x,y)
  lfdata <- na.omit(lfdata)
  x <- lfdata$x; y <- lfdata$y
  locboot <- locfit::locfit(y ~ locfit::lp(x, deg=1, h=hw), maxk=800, maxit=20, family="qrgauss")
  predboot <- predict(locboot, newdata=targage.df, se.fit=TRUE)
  yfit[,i] <- predboot$fit
  lines(targage.df$x, yfit[,i], lwd=2, col=rgb(0.5,0.5,0.5,0.10))
  if (i %% 10 == 0) {print(i)}
}
warnings()


# Step 3 -- Plot the unresampled (initial) fit

#bw_50
fitname_bw_50 <- paste("Zlocfit_",as.character(bw_50), sep="")
colnames(loc01_fit_bw_50) <- c("age_bw_50", fitname_bw_50)
lines(loc01_fit_bw_50[,1], loc01_fit_bw_50[,2], lwd=2, col="blue")

#bw_200
fitname_bw_200 <- paste("Zlocfit_",as.character(bw_200), sep="")
colnames(loc01_fit_bw_200) <- c("age_bw_200", fitname_bw_200)
lines(loc01_fit_bw_200[,1], loc01_fit_bw_200[,2], lwd=2, col="black")

#bw_500
fitname_bw_500 <- paste("Zlocfit_",as.character(bw_500), sep="")
colnames(loc01_fit_bw_500) <- c("age_bw_500", fitname_bw_500)
lines(loc01_fit_bw_500[,1], loc01_fit_bw_500[,2], lwd=2, col="green")

#Normal_hw
fitname <- paste("Zlocfit_",as.character(hw), sep="")
colnames(loc01_fit) <- c("age", fitname)
lines(loc01_fit[,1], loc01_fit[,2], lwd=2, col="red") 


# Step 4 -- Find and add bootstrap CIs
lines(targage.df$x, yfit95, lwd=1, col="red")
lines(targage.df$x, yfit05, lwd=1, col="red")



if (plotout == "pdf") {dev.off()}




# ---------------------------------------------------------


# 6. Save data for correlation analysis
# ---------------------------------------------------------
# ---------------------------------------------------------

rio::export(lfdata_corr, "./other_output/rpd/lfdata_corr.csv")
rio::export(loc01_fit, "./other_output/rpd/loc01_fit.csv")


# ---------------------------------------------------------


# 7. Generate curves based on alternative specifications
# ---------------------------------------------------------
# ---------------------------------------------------------


# During data import (see data_import.R) and throught the
# code above, there are options to specify changes to the
# data input and hence to the output of the above composite
# curve generation. This necessitates storing the rpd output
# manually within a further csv file. The code below
# indicates how such a approach could be followed. Once the
# approriate data has been saved, with the "type" entered,
# the following code generates curves showing:
# 1.  the number of records prooviding data per data bin;
# 2.  the influence of minimax or max transformation;
# 3.  the influence of higher elevation sites;
# 4.  the impact of a preference for micro rather than macro
#     charcoal records

# Generate graphs of multiple composites
multiple_comp <- rio::import("./other_output/rpd/multiple_composites.csv") %>% #Import manually generated file from saved composite runs
  dplyr::filter(between(age, 3500, 10000))

multiple_comp_records <- multiple_comp %>% 
  dplyr::filter(type == "normal")

multiple_comp_trans <- multiple_comp %>% #Minmax or Max transformation
  dplyr::filter(type %in% c("minmax", "normal")) %>%
  dplyr::mutate(type = if_else(type == "minmax", "minmax", "max"))

multiple_comp_elevation_high <- multiple_comp %>%
  dplyr::filter(type %in% c("only_high_elevation", "normal")) %>%
  dplyr::mutate(type = if_else(type == "only_high_elevation", "higher elevation sites", "excluding higher elevation sites"))

multiple_comp_micro <- multiple_comp %>%
  dplyr::filter(type %in% c("micro_instead", "normal")) %>%
  dplyr::mutate(type = if_else(type == "micro_instead", "micro-charcoal records favoured", "macro- charcoal records favoured"))


ggplot2::ggplot(data = multiple_comp_records, mapping = aes(x = age, y = nsites))+
  ggplot2::geom_bar(stat = "identity")+
  # ylim(-3, 1)+
  scale_x_reverse()+
  labs(x ="Years cal. BP", y ="Number of records per 100 year bin")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/Composites_records.pdf", height = 8, width = 11.28, units = "cm")


ggplot2::ggplot(data = multiple_comp_trans, mapping = aes(x = age, y = locfit, color = type))+
  ggplot2::geom_line()+
  ylim(-0.5, 0.2)+
  scale_x_reverse()+
  labs(color = "Transformation type", x ="Years cal. BP", y ="Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/Composites_trans.pdf", height = 8, width = 13, units = "cm") #Minmax or max


ggplot2::ggplot(data = multiple_comp_elevation_high, mapping = aes(x = age, y = locfit, color = type))+
  ggplot2::geom_line()+
  ylim(-3, 1)+
  scale_x_reverse()+
  labs(color = "Elevation", x ="Years cal. BP", y ="Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/Composites_elevation_high.pdf", height = 8, width = 15, units = "cm") #Minmax or max


ggplot2::ggplot(data = multiple_comp_micro, mapping = aes(x = age, y = locfit, color = type))+
  ggplot2::geom_line()+
  ylim(-3, 1)+
  scale_x_reverse()+
  labs(color = "Macro- or micro or charcoal records", x ="Years cal. BP", y ="Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/Composites_micro.pdf", height = 8, width = 15, units = "cm") #Micro or macro

# ---------------------------------------------------------





# 8. Changepoint analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
change_char <- loc01_fit %>% 
  dplyr::filter(between(age, 3500, 10000)) 
char_cp_binseg <- changepoint::cpt.mean(change_char$Zlocfit_500, method = "BinSeg")
plot(char_cp_binseg, cpt.col = "blue", ylim = c(-0.5,0.5))


# ---------------------------------------------------------



# 9. Trend in charcoal data 
# ---------------------------------------------------------
# ---------------------------------------------------------

filtered_loc01_fit <- loc01_fit %>% 
  dplyr::filter(dplyr::between(age, 3500, 10000)) #Define data of interest
fit_loc01_fit <- lm(Zlocfit_500 ~ age, data = filtered_loc01_fit) #Fit linear model to composite curve
linear_loc01_fit <- dplyr::tibble(age = filtered_loc01_fit$age, Zlocfit_500 = fitted(fit_loc01_fit))
detrend_loc01_fit <- dplyr::tibble(age = filtered_loc01_fit$age, detrend_Zlocfit_500 = resid(fit_loc01_fit)) #Extract residuals


#Plot detrended median (residuals)
ggplot2::ggplot(data = loc01_fit, mapping = aes(y = Zlocfit_500, x = age))+
  geom_line(color = "red")+
  geom_line(data = linear_loc01_fit, mapping = aes(y = Zlocfit_500, x = age))+
  coord_cartesian(xlim = c(10000, 3500), ylim = c(-0.7, 0.3))+
  labs(x ="Years cal. BP", y ="Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/rpd_detrend.pdf", height = 8, width = 11.28, units = "cm")

#Plot detrended median (residuals)
ggplot2::ggplot(data = detrend_loc01_fit, mapping = aes(y = detrend_Zlocfit_500, x = age))+
  geom_line()+
  coord_cartesian(xlim = c(10000, 3500), ylim = c(-0.7, 0.3))+
  labs(x ="Years cal. BP", y ="Detrended charcoal z-scores\n(fitted model residuals)")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggsave("figs/rpd/rpd_detrend_residuals.pdf", height = 8, width = 11.28, units = "cm")

# ---------------------------------------------------------
  

# 10. Site information 
# ---------------------------------------------------------
# ---------------------------------------------------------
charcoal_site_summary <- rpd_data_input %>% 
    dplyr::filter(dplyr::between(INTCAL2020_median, 3500, 10000)) %>% 
    dplyr::group_by(site_name) %>% 
    dplyr::summarise(latitude = round(mean(latitude), 2), 
                     longitude = round(mean(longitude), 2), 
                     elevation = round(mean(elevation)), 
                     macro_micro = (unique(macro_micro)),
                     record_length = round((max(INTCAL2020_median) - min(INTCAL2020_median))/1000,2),
                     resolution = n(),
                     ID_ENTITY = mean(ID_ENTITY)) %>% 
    dplyr::mutate(resolution = round((resolution / record_length), 2)) %>% 
    dplyr::left_join(site_bin_info_n, by = "ID_ENTITY")




