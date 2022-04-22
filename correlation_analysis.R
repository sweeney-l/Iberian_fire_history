#Sweeney, Harrison and Vander Linden 2022. 
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




# 2. Setup the data for analysis
# ---------------------------------------------------------

# Load in info from charcoal and radiocarbon code
spd_data_corr_sm <- rio::import("./other_output/pop/spd_data_sm.csv") %>% #Load in data for joint figure
  dplyr::mutate(group = ggplot2::cut_width(age, width = 100, centre = 3500, labels = F)) %>%  #Adjust into 100 year bins
  dplyr::group_by(group) %>% 
  dplyr::summarise(spd = mean(spd)) %>% 
  dplyr::mutate(age = seq(3500, 10000, 100)) %>% 
  dplyr::select(-group)

spd_data_corr <- rio::import("./other_output/pop/spd_data.csv") %>% #Load in data
  dplyr::mutate(group = ggplot2::cut_width(age, width = 100, centre = 3500, labels = F)) %>%  #Adjust into 100 year bins
  dplyr::group_by(group) %>% 
  dplyr::summarise(spd = mean(spd)) %>% 
  dplyr::mutate(age = seq(3500, 10000, 100)) %>% 
  dplyr::select(-group)

char_smoothed <- rio::import("./other_output/rpd/loc01_fit.csv") %>% #Load in data for joint figure
  dplyr::group_by(age) %>% 
  dplyr::summarise(char = mean(Zlocfit_500)) %>% 
  dplyr::filter(dplyr::between(age, 3500, 10000))

char_not_smoothed <- rio::import("./other_output/rpd/lfdata_corr.csv") %>% #Load in data 
  dplyr::group_by(x) %>% 
  dplyr::summarise(char = mean(y)) %>% 
  dplyr::rename(age = x) %>% 
  dplyr::filter(dplyr::between(age, 3500, 10000))

# Set up joint charcoal and population file
corr_char_spd <- char_not_smoothed %>% 
  dplyr::left_join(spd_data_corr, by = "age")  
  
# Establish quartile levels of charcoal and spd population
q_spd <- as.numeric(quantile(corr_char_spd$spd))
q_char <- as.numeric(quantile(corr_char_spd$char))


# ---------------------------------------------------------



# 3. Quartile and spd changepoint plots
# ---------------------------------------------------------

# Charcoal quartiles
ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+ 
  geom_line()+
  geom_hline(yintercept = c(q_char[2], q_char[3],q_char[4]), linetype = "dashed", color = c("red", "blue", "red"))+
  scale_x_reverse()

# spd population quartiles
ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = spd, x = age))+
  geom_line()+
  geom_hline(yintercept = c(q_spd[2], q_spd[3],q_spd[4]), linetype = "dashed", color = c("red", "blue", "red"))+
  scale_x_reverse()

# spd population changepoint sections
change_spd <- spd_rm200$grid$PrDens #Select probability densities
spd_cp_binseg_mean <- changepoint::cpt.mean(change_spd, method = "BinSeg", "MBIC") #Run changepoint analysis on mean values

first_spd_changepoint <- 10000 - spd_cp_binseg_mean@cpts[1] #First changepoint location. 10k minus because of BP dates
second_spd_changepoint <- 10000 - spd_cp_binseg_mean@cpts[2] #Second changepoint location. 10k minus because of BP dates

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = spd, x = age))+
  geom_line()+
  geom_vline(xintercept = c(first_spd_changepoint, second_spd_changepoint), linetype = "dashed", color = c("red", "blue"))+
  scale_x_reverse()



# Define periods when charcoal and spd population above or below quartiles
# Charcoal years
yr_char_above_75 <- corr_char_spd %>% 
  dplyr::filter(char>q_char[4]) %>% 
  dplyr::pull(age)
yr_char_below_25 <- corr_char_spd %>% 
  dplyr::filter(char<q_char[2]) %>% 
  dplyr::pull(age)
yr_char_between_25_75 <- corr_char_spd %>%
  dplyr::filter(dplyr::between(char, q_char[2], q_char[4])) %>% 
  dplyr::pull(age)

# spd population years
yr_spd_above_75 <- corr_char_spd %>% 
  dplyr::filter(spd > q_spd[4]) %>% 
  dplyr::pull(age)
yr_spd_below_25 <- corr_char_spd %>% 
  dplyr::filter(spd<q_spd[2]) %>% 
  dplyr::pull(age)
yr_spd_between_25_75 <- corr_char_spd %>%
  dplyr::filter(dplyr::between(spd, q_spd[2], q_spd[4])) %>% 
  dplyr::pull(age)



# Figures
# Charcoal with times when spd population above, below and centre quartiles
ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_spd_above_75), linetype = "dashed", color = "red")+
  scale_x_reverse()

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_spd_below_25), linetype = "dashed", color = "red")+
  scale_x_reverse()

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_vline(xintercept = c(yr_spd_between_25_75), linetype = "dashed", color = "red")+
  scale_x_reverse()


# spd population with times when charcoal above and below quartiles
ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = spd, x = age))+
  geom_line()+
  geom_vline(xintercept = yr_char_above_75[1], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[2], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[8], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[9], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[10], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[13], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[14], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[15], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[16], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_above_75[17], color = "blue", alpha = 0.4)+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_above_75[3], xmax = yr_char_above_75[4], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_above_75[5], xmax = yr_char_above_75[7], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_above_75[11], xmax = yr_char_above_75[12], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal. BP", y = "Summed Probability")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/spd_char_correlation_above_75.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = spd, x = age))+
  geom_line()+
  geom_vline(xintercept = yr_char_below_25[1], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_below_25[4], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_below_25[16], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_below_25[17], color = "blue", alpha = 0.4)+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_below_25[2], xmax = yr_char_below_25[3], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_below_25[5], xmax = yr_char_below_25[9], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_below_25[10], xmax = yr_char_below_25[11], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_below_25[12], xmax = yr_char_below_25[13], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_below_25[14], xmax = yr_char_below_25[15], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal. BP", y = "Summed Probability")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/spd_char_correlation_below_25.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = spd, x = age))+
  geom_line()+
  geom_vline(xintercept = yr_char_between_25_75[15], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[16], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[22], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[23], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[27], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[28], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[29], color = "blue", alpha = 0.4)+
  geom_vline(xintercept = yr_char_between_25_75[32], color = "blue", alpha = 0.4)+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[1], xmax = yr_char_between_25_75[6], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[7], xmax = yr_char_between_25_75[9], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[10], xmax = yr_char_between_25_75[11], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[12], xmax = yr_char_between_25_75[14], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[17], xmax = yr_char_between_25_75[18], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[19], xmax = yr_char_between_25_75[21], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[24], xmax = yr_char_between_25_75[26], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[1,],
            aes(xmin = yr_char_between_25_75[30], xmax = yr_char_between_25_75[31], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal. BP", y = "Summed Probability")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/spd_char_correlation_between_25_75.pdf", height = 8, width = 11.28, units = "cm")


# Charcoal with times when spd above and below quartiles
ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_rect(data = corr_char_spd[2,],
            aes(xmin = yr_spd_above_75[1], xmax = yr_spd_above_75[14], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  geom_rect(data = corr_char_spd[2,],
            aes(xmin = yr_spd_above_75[15], xmax = yr_spd_above_75[17], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal. BP", y = "Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/char_spd_correlation_above_75.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_rect(data = corr_char_spd[2,],
            aes(xmin = yr_spd_below_25[1], xmax = yr_spd_below_25[17], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/char_spd_correlation_below_25.pdf", height = 8, width = 11.28, units = "cm")

ggplot2::ggplot(data = corr_char_spd, mapping= aes(y = char, x = age))+
  geom_line()+
  geom_vline(xintercept = yr_spd_between_25_75[1], color = "blue", alpha = 0.4)+
  geom_rect(data = corr_char_spd[2,],
            aes(xmin = yr_spd_between_25_75[2], xmax = yr_spd_between_25_75[32], ymin = -Inf, ymax = Inf), alpha = 0.4, fill = "blue")+
  scale_x_reverse()+
  ggplot2::labs(x = "Years cal BP", y = "Z-scores of transformed charcoal influx")+
  theme(plot.margin = unit(c(4,4,1,1), "mm"), panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "gray", linetype = 3)) +
  theme(text = element_text(family = "Times", size = 12))
  ggplot2::ggsave("./figs/cor/char_spd_correlation_between_25_75.pdf", height = 8, width = 11.28, units = "cm")

# ---------------------------------------------------------



# 4. Correlation plots
# ---------------------------------------------------------

# Visual comparison of charcoal and spd over time
pdf("./figs/cor/spd_char_simple.pdf", height = 3.15, width = 4.5, pointsize = 9, family = "Times")
par(mar = c(4, 4, 1, 4))
plot(fig_char$age, 
       fig_char$locfit, 
       xlab="", 
       ylab="", 
       xlim=c(10000,3500), 
       ylim=c(-4,1), 
       type="n", 
       cex.lab = 1.2, 
       mgp = c(2.4, 1, 0), 
       axes = F)
axis(side = 2, at = c(-3,-2, -1, 0, 1))
mtext("Z-scores of transformed charcoal influx", side = 2, line=2.5, at=-1.5, cex = 1.2)
axis(side = 1, xpd = TRUE)
mtext("Years cal. BP", side = 1, line=2.5, at=7000, cex = 1.2)
polygon(c(fig_char$age, rev(fig_char$age)), c(fig_char$cu95, rev(fig_char$cl95)),  col = "grey90", lty = 0)
lines(fig_char$age, fig_char$locfit, lwd = 2, col = "red")
grid()
  
par(new=T)
plot(fig_spd$calBP,
       fig_spd$PrDens,
       xlim=c(10000,3500), 
       ylim=c(0, 3),
       type="n",
       xlab = "",
       ylab = "",
       axes = F)
axis(side = 4, at = c(0, 0.5, 1, 1.5))
polygon(c(fig_spd$calBP, rev(fig_spd$calBP)), c(fig_spd$lo, rev(fig_spd$hi)),  col = "grey90", lty = 0)
mtext("Summed probability", side = 4, line=2.5, at=0.75, cex = 1.2)
lines(fig_spd$calBP, fig_spd$PrDens, lwd = 2, col = "blue")
dev.off()
  
# Correlation with quartile ranges shown for charcoal and spd population
ggplot2::ggplot(data = corr_char_spd, mapping = aes(spd, char))+
  geom_point()+
  geom_hline(yintercept = c(q_char[2],q_char[4]), linetype = "dashed", color = c("red", "red"))+
  geom_vline(xintercept = c(q_spd[2], q_spd[4]), linetype = "dashed", color = c("blue", "blue"))

# Correlation with regression line
ggpubr::ggscatter(corr_char_spd, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.65, label.y.npc = 0.25)+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/all_correlation.pdf", height = 8, width = 11.28, units = "cm")

# Correlation with loess smoother
ggpubr::ggscatter(corr_char_spd, x = "spd", y = "char",  add = "loess", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/all_correlation_loess.pdf", height = 8, width = 11.28, units = "cm")
 
 
# Correlation for quartiles spd population 
# Below 25th percentile
corr_char_spd_below_25_spd <- corr_char_spd %>% 
  dplyr::filter(spd<=q_spd[2])

ggpubr::ggscatter(corr_char_spd_below_25_spd, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_spd_below_25.pdf", height = 8, width = 11.28, units = "cm")

# Above 75th percentile
corr_char_spd_above_75_spd <- corr_char_spd %>% 
  dplyr::filter(spd>q_spd[4])

ggpubr::ggscatter(corr_char_spd_above_75_spd, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_spd_above_75.pdf", height = 8, width = 11.28, units = "cm")

# Between 25th and 75th percentile
corr_char_spd_between_25_and_75_spd <- corr_char_spd %>% 
  dplyr::filter(dplyr::between(spd, q_spd[2], q_spd[4]))

ggpubr::ggscatter(corr_char_spd_between_25_and_75_spd, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_spd_middle.pdf", height = 8, width = 11.28, units = "cm")


# Correlation for quartiles charcoal 
# Below 25th percentile
corr_char_spd_below_25_char <- corr_char_spd %>% 
  dplyr::filter(char<=q_char[2])

ggpubr::ggscatter(corr_char_spd_below_25_char, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_char_below_25.pdf", height = 8, width = 11.28, units = "cm")

# Above 75th percentile
corr_char_spd_above_75_char<- corr_char_spd %>% 
  dplyr::filter(char>q_char[4])

ggpubr::ggscatter(corr_char_spd_above_75_char, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_char_above_75.pdf", height = 8, width = 11.28, units = "cm")

# Between 25th and 75th percentile
corr_char_spd_between_25_and_75_char <- corr_char_spd %>% 
  dplyr::filter(dplyr::between(char, q_char[2], q_char[4]))

ggpubr::ggscatter(corr_char_spd_between_25_and_75_char, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/correlation_char_middle.pdf", height = 8, width = 11.28, units = "cm")


#Correlation based on spd changepoint sections
corr_char_spd_first_change <- corr_char_spd  %>% 
  dplyr::filter(age>first_spd_changepoint)

ggpubr::ggscatter(corr_char_spd_first_change, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/all_correlation_first_change.pdf", height = 8, width = 11.28, units = "cm")


corr_char_spd_second_change <- corr_char_spd  %>% 
  dplyr::filter(dplyr::between(age, second_spd_changepoint, first_spd_changepoint))

ggpubr::ggscatter(corr_char_spd_second_change, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/all_correlation_second_change.pdf", height = 8, width = 11.28, units = "cm")


corr_char_spd_third_change <- corr_char_spd  %>% 
  dplyr::filter(age<second_spd_changepoint)

ggpubr::ggscatter(corr_char_spd_third_change, x = "spd", y = "char",  add = "reg.line", conf.int = T)+
  ggplot2::labs(x = "Summed Probability", y = "Z-score composite of charcoal influx")+
  ggpubr::stat_cor()+
  ggpubr::labs_pubr(base_size = 10)
  ggplot2::ggsave("./figs/cor/all_correlation_third_change.pdf", height = 8, width = 11.28, units = "cm")


# ---------------------------------------------------------

# Rolling correlation analysis
a_corr_char_spd <- corr_char_spd %>% 
  dplyr::arrange(desc(age))
roll_corr <- zoo::rollapply(a_corr_char_spd, width = 10, function(x) cor(x[,2], x[,3]), by.column = F) 
roll_corr_p <- zoo::rollapply(a_corr_char_spd, width = 10, function(x) cor.test(x[,2], x[,3])$p.value, by.column = F)
  
roll_corr_char_spd <- tail(a_corr_char_spd, -9) %>% 
  cbind(roll_corr, roll_corr_p)

ggplot2::ggplot(data = roll_corr_char_spd, mapping = aes(x = age, y = roll_corr))+
  geom_line()+
  geom_point(data = dplyr::filter(roll_corr_char_spd, roll_corr_p<0.05), mapping = aes(x = age, y = roll_corr), color = "blue") +
  ggplot2::labs(x = "Years cal. BP", y = "Rolling correlation spd and charcoal influx")+
  scale_x_reverse()+
  theme(legend.position="none")




