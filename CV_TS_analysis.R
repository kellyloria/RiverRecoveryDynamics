############
## Time series analysis aggregation of all stream gauge and daily time sieries data
##
## Author: Kelly A. Loria
## Date Created: 2025-01-15
## Email: kellyloria @ gmail.com
##
## ---------------------------
#PC: setwd("R:/Users/kloria/Documents/River_Recovery_Analysis")
#Mac: setwd("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis")

## ---------------------------
## meta data for watershed/landscape differences: Metab_metadata.rds
##    complied in "Metab_watershed_metadata.R
## ---------------------------
## looking at worflow from: https://github.com/nmarzolf91/usgs_metabolism_project
## need to check out https://github.com/nmarzolf91/usgs_metabolism_project/tree/main/code/9_intraannual_variation
## ---------------------------
library(dataRetrieval)
library(nhdplusTools)
library(scales)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(lubridate)

library(PerformanceAnalytics)
library(MuMIn)
library(lme4)
library(lmerTest)
library(car)
## ---------------------------

## ---------------------------
## Load time series data:
metab_ts <- readRDS("./data/Metab_TS.rds")
str(metab_ts)

dat_ts <- readRDS("./data/TS_SPC_pH_FNU.rds")
str(dat_ts)
summary(dat_ts)
unique(dat_ts$site_name)


metab_meta_dat <- readRDS("./data/Metab_metadata.rds")
str(metab_ts)

marz_dat <- read.csv("./data/doi_10_5061_dryad_bcc2fqzj2__v20230622/data_citation/1_site_info.csv")

# 63 different sites 
## ---------------------------
## Step 1. data cleaning 
##  
dat_ts_clean <- dat_ts%>%
  filter(GPP.Rhat < 1.05 &  GPP >0 & ER <0)

unique(dat_ts_clean$long_name)

marz_dat_names <- marz_dat%>%
  select(Site.Number, Latitude, Longitude) %>%
  mutate(source="marzolf")

dat_ts_names <-metab_meta_dat %>%
  select(Site.Number = "site_no", Latitude = "dec_lat_va", Longitude= "dec_long_va") %>%
  mutate(source="appling")

## quick site check
site_check <- rbind(marz_dat_names, dat_ts_names)

unique(site_check)



library(tidyverse)
library(ggrepel)
library(ggmap)
library(ggsci)
library(scales)

library(dplyr)
library(ggplot2)
library(ggrepel)

# devtools::install_github('oswaldosantos/ggsn')
library(ggsn)

#devtools::install_github("briatte/tidykml")
library(tidykml)
library(sf)



range(na.omit(site_check$Longitude))
range(site_check$Latitude)
# 
# Define the bounding box coordinates
bbox <- c(left = -127.300, bottom = 24.900, right = -70.800, top = 52.30)

# Retrieve terrain map layer using get_stamenmap()
map <- get_stadiamap(bbox = bbox, zoom = 5, maptype = 'stamen_terrain')
# need api run for the first time

#library(ggmap)
# Plot the terrain map
stamen_map <- ggmap(map)
stamen_map





coords_map <- ggmap(map) +
  geom_point(data = site_check, aes(x = Longitude, y = Latitude, color = source), size = 2) +
  geom_label_repel(data = site_check, aes(x = Longitude, y = Latitude, label = Site.Number),
                   color = "black", fontface = "bold", size = 3) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() +
  scale_color_viridis_d(option = "viridis",  name = "Water year")  +
  # Add the scale bar
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16, colour = "black", face = "bold"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

