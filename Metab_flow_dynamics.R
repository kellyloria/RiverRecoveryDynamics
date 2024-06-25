############

setwd("R:/Users/kloria/Documents/River_Recovery_Analysis")

metab_ts <- readRDS("./data/Metab_TS.rds")
# Link in other watershed attributes

str(metab_ts)

###########
###########
# Download stream flow from each site 
## ---------------------------
## Data aggregation of all stream gauge data
##
## Author: Kelly A. Loria
## Date Created: 2020-09-26
## Email: kelly.loria@nevada.unr.edu
##
## ---------------------------
## Load packages:
#install.packages("dataRetrieval")
#library(remotes)
#install_github("USGS-R/dataRetrieval", 
#               build_opts = c("--no-resave-data", "--no-manual"),
#               build_vignettes = TRUE)

library(dataRetrieval)
library(nhdplusTools)
library(ggplot2)
library(scales)
library(tidyverse)
#library(zoo)

# File path setup:
if (dir.exists('/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/')){
  inputDir<- '/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/'
  outputDir<- '/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_output/2023_output/' 
}


## ---------------------------
# Source Data from Water Quality Portal
#     List of all values:
#     https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&inline=true&group_cd=%
## ---------------------------
# 10336730 -GB
# IN - 10336700
# TH - 10336698
# WD - 10336676
# BW - 10336660
# GN - 10336645

# TRT - 10336780
# UPTL - 10336610
# UPTU - 103366092

### 1
siteNumber <- "10336730" 
Info <- readNWISsite(siteNumber)
parameterCd <- c("00060")  # discharge

#Raw daily data:
qwData1 <- readNWISdv(siteNumber,parameterCd,
                      "1980-09-30","2023-09-20")

qwData1$dec_lat_va <- Info$dec_lat_va
qwData1$dec_long_va <- Info$dec_long_va
qwData1$huc_cd <- Info$huc_cd
qwData1$drain_area_va <- Info$drain_area_va
qwData1$project_no <- Info$project_no
qwData1$station_nm <- Info$station_nm
qwData1$alt_va <- Info$alt_va

names(Info)



### 2
siteNumber <- "10336660" 
Info <- readNWISsite(siteNumber)
parameterCd <- c("00060")  # discharge

#Raw daily data:
qwData1 <- readNWISdv(siteNumber,parameterCd,
                      "1980-09-30","2023-09-20")

qwData1$dec_lat_va <- Info$dec_lat_va
qwData1$dec_long_va <- Info$dec_long_va
qwData1$huc_cd <- Info$huc_cd
qwData1$drain_area_va <- Info$drain_area_va
qwData1$project_no <- Info$project_no
qwData1$station_nm <- Info$station_nm
qwData1$alt_va <- Info$alt_va

names(Info)

nldi_nwis <- list(featureSource = "nwissite", featureID = "USGS-103366092")
discover_nhdplus_id(nldi_feature = nldi_nwis)

