#Sweeney, Harrison and Vander Linden 2021. 
#Assessing anthropogenic influence on fire history during the Holocene in the Iberian Peninsula
# ---------------------------------------------------------

#There are six seperate scripts associated with this research:
# 1. data_import.R 
# 2. charcoal_analysis.R 
# 3. radiocarbon_analysis.R
# 4. correlation_analysis.R (this script)
# 5. SEA_analysis.R
# 6. neolithic_analysis.R


# ---------------------------------------------------------



# Correlation analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
# This script uses the data generated within the
# charcoal_analyis.R and radiocarbon_analysis.R scrips to 
# explore the correlation between fire and population
# during the 10-3.5k period within mainland Iberia. 



# Script elements
# ---------------------------------------------------------
# 1. Packages, paths and data
# 2. Setup the data for analysis
# 3. Quartile plots
# 4. Correlation plots



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




# 2. Setup the data for analysis
# ---------------------------------------------------------

# Load in info from charcoal and radiocarbon code
detrend_loc01_fit <- rio::import("./other_output/rpd/detrend_loc01_fit.csv") #Fitted charcoal composite data
roc_blocks <- rio::import("./other_output/pop/roc_blocks.csv") #ROC of spd

# Fit ROC with same locfit specification as charcoal data
roc_loc <- locfit::locfit(roc_blocks$roc ~ locfit::lp(roc_blocks$age, deg=1, h=500), maxk=800, family="qrgauss") #Apply same locfit smoothig procedure
roc_pred <- predict(roc_loc, newdata=data.frame(roc_blocks$age), se.fit=TRUE)
roc_loc_fit <- data.frame(roc_blocks$age, roc_pred$fit)
fitname <- paste("locfit_",500, sep="")
colnames(roc_loc_fit) <- c("age", fitname)

# Set up joint charcoal and population file
corr_char_rocsm <- detrend_loc01_fit %>% 
  dplyr::rename(char_locfit = 2)  %>% 
  dplyr::left_join(roc_loc_fit, by = "age") %>% 
  dplyr::rename(roc_locfit = locfit_500) %>% 
  dplyr::filter(!is.na(roc_locfit)) #Filter data to same age range

# Establish quartile levels of charcoal and ROC population
q_rocsm <- as.numeric(quantile(corr_char_rocsm$roc_locfit))
q_charsm <- as.numeric(quantile(corr_char_rocsm$char_locfit))


# ---------------------------------------------------------



# 3. Quartile plots
# ---------------------------------------------------------

# Charcoal
ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+ 
  geom_line()+
  geom_hline(yintercept = c(q_charsm[2], q_charsm[3],q_charsm[4]), linetype = "dashed", color = c("red", "blue", "red"))+
  scale_x_reverse()

# ROC population
ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = roc_locfit, x = age))+
  geom_line()+
  geom_hline(yintercept = c(q_rocsm[2], q_rocsm[3],q_rocsm[4]), linetype = "dashed", color = c("red", "blue", "red"))+
  scale_x_reverse()

# Define periods when detrended charcoal and ROC population above or below quartiles
# Charcoal years
yr_char_above_75sm <- corr_char_rocsm %>% 
  dplyr::filter(char_locfit>q_charsm[4]) %>% 
  dplyr::pull(age)
yr_char_below_25sm <- corr_char_rocsm %>% 
  dplyr::filter(char_locfit<q_charsm[2]) %>% 
  dplyr::pull(age)
yr_char_between_25_75sm <- corr_char_rocsm %>%
  dplyr::filter(dplyr::between(char_locfit, q_charsm[2], q_charsm[4])) %>% 
  dplyr::pull(age)

# ROC population years
yr_roc_above_75sm <- corr_char_rocsm %>% 
  dplyr::filter(roc_locfit>q_rocsm[4]) %>% 
  dplyr::pull(age)
yr_roc_below_25sm <- corr_char_rocsm %>% 
  dplyr::filter(roc_locfit<q_rocsm[2]) %>% 
  dplyr::pull(age)
yr_roc_between_25_75sm <- corr_char_rocsm %>%
  dplyr::filter(dplyr::between(roc_locfit, q_rocsm[2], q_rocsm[4])) %>% 
  dplyr::pull(age)


# Figures
# Detrended charcoal locfit with times when ROC population above, below and centre quartiles
ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_roc_above_75sm), linetype = "dashed", color = "red")+
  scale_x_reverse()

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_roc_below_25sm), linetype = "dashed", color = "red")+
  scale_x_reverse()

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_roc_between_25_75sm), linetype = "dashed", color = "red")+
  scale_x_reverse()


# ROC population locfit with times when detrended charcoal above and below quartiles
ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = roc_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_above_75sm[1], xmax = yr_char_above_75sm[8], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_above_75sm[9], xmax = yr_char_above_75sm[11], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_above_75sm[12], xmax = yr_char_above_75sm[16], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Z-score rate of change of SPD")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/rocsm_char_correlation_above_75_roc.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = roc_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_below_25sm[1], xmax = yr_char_below_25sm[5], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_below_25sm[6], xmax = yr_char_below_25sm[12], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_below_25sm[13], xmax = yr_char_below_25sm[16], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Z-score rate of change of SPD")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/rocsm_char_correlation_below_25_roc.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = roc_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[1], xmax = yr_char_between_25_75sm[9], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[10], xmax = yr_char_between_25_75sm[11], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[12], xmax = yr_char_between_25_75sm[19], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[20], xmax = yr_char_between_25_75sm[21], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[22], xmax = yr_char_between_25_75sm[27], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[1,],
            aes(xmin = yr_char_between_25_75sm[28], xmax = yr_char_between_25_75sm[31], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Z-score rate of change of SPD")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/rocsm_char_correlation_between_25_75_roc.pdf", height = 8, width = 11.28, units = "cm")


# Detrended charcoal with times when ROC locfit above and below quartiles
ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_above_75sm[1], xmax = yr_roc_above_75sm[4], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_above_75sm[5], xmax = yr_roc_above_75sm[5], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_above_75sm[6], xmax = yr_roc_above_75sm[10], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_above_75sm[11], xmax = yr_roc_above_75sm[12], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_above_75sm[13], xmax = yr_roc_above_75sm[16], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Detrended charcoal z-scores")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/charsm_roc_correlation_above_75_char.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_below_25sm[1], xmax = yr_roc_below_25sm[8], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_below_25sm[9], xmax = yr_roc_below_25sm[14], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_below_25sm[15], xmax = yr_roc_below_25sm[16], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Detrended charcoal z-scores")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/charsm_roc_correlation_below_25_char.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_rocsm, mapping= aes(y = char_locfit, x = age))+
  geom_line()+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[1], xmax = yr_roc_between_25_75sm[7], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[8], xmax = yr_roc_between_25_75sm[15], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[16], xmax = yr_roc_between_25_75sm[17], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[18], xmax = yr_roc_between_25_75sm[20], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[21], xmax = yr_roc_between_25_75sm[30], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  geom_rect(data = corr_char_rocsm[2,],
            aes(xmin = yr_roc_between_25_75sm[32], xmax = yr_roc_between_25_75sm[32], ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Detrended charcoal z-scores")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))+
  ggplot2::ggsave("./figs/cor/charsm_roc_correlation_between_25_75_char.pdf", height = 8, width = 11.28, units = "cm")

# ---------------------------------------------------------



# 4. Correlation plots
# ---------------------------------------------------------

# Correlation with quartile ranges shown for charcoal and ROC population
ggplot2::ggplot(data = corr_char_rocsm, mapping = aes(roc_locfit, char_locfit))+
  geom_point()+
  geom_hline(yintercept = c(q_charsm[2],q_charsm[4]), linetype = "dashed", color = c("red", "red"))+
  geom_vline(xintercept = c(q_rocsm[2], q_rocsm[4]), linetype = "dashed", color = c("blue", "blue"))

# Correlation with regression line
ggpubr::ggscatter(corr_char_rocsm, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)

# Correlation with loess smoother
ggpubr::ggscatter(corr_char_rocsm, x = "roc_locfit", y = "char_locfit",  add = "loess", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::labs_pubr(base_size = 20)
 
 
# Correlation for quartiles ROC population 
# Below 25th percentile
corr_char_rocsm_below_25_roc <- corr_char_rocsm %>% 
  dplyr::filter(roc_locfit<=q_rocsm[2])

ggpubr::ggscatter(corr_char_rocsm_below_25_roc, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)

# Above 75th percentile
corr_char_rocsm_above_75_roc <- corr_char_rocsm %>% 
  dplyr::filter(roc_locfit>q_rocsm[4])

ggpubr::ggscatter(corr_char_rocsm_above_75_roc, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)

# Between 25th and 75th percentile
corr_char_rocsm_between_25_and_75_roc <- corr_char_rocsm %>% 
  dplyr::filter(dplyr::between(roc_locfit, q_rocsm[2], q_rocsm[4]))

ggpubr::ggscatter(corr_char_rocsm_between_25_and_75_roc, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)


# Correlation for quartiles charcoal 
# Below 25th percentile
corr_char_rocsm_below_25_char <- corr_char_rocsm %>% 
  dplyr::filter(char_locfit<=q_charsm[2])

ggpubr::ggscatter(corr_char_rocsm_below_25_char, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended charcoal accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)

# Above 75th percentile
corr_char_rocsm_above_75_char<- corr_char_rocsm %>% 
  dplyr::filter(char_locfit>q_charsm[4])

ggpubr::ggscatter(corr_char_rocsm_above_75_char, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended z-score composite of charcoal influx")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)

# Between 25th and 75th percentile
corr_char_rocsm_between_25_and_75_char <- corr_char_rocsm %>% 
  dplyr::filter(dplyr::between(char_locfit, q_charsm[2], q_charsm[4]))

ggpubr::ggscatter(corr_char_rocsm_between_25_and_75_char, x = "roc_locfit", y = "char_locfit",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "SPD rate of change", y = "Detrended accumulation composite")+
  ggpubr::stat_cor(size = 8)+
  ggpubr::labs_pubr(base_size = 20)



# ---------------------------------------------------------



