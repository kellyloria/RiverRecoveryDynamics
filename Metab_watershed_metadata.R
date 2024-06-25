############
## Data aggregation of all stream gauge data
##
## Author: Kelly A. Loria
## Date Created: 2024-06-26
## Email: kellyloria@gmail.com
##
## ---------------------------
#PC: setwd("R:/Users/kloria/Documents/River_Recovery_Analysis")
#Mac: setwd("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis")

metab_ts <- readRDS("./data/Metab_TS.rds")
str(metab_ts)
###########
# Download stream flow + watershed attributes from each site 
## ---------------------------
library(dataRetrieval)
library(nhdplusTools)
library(scales)
library(tidyverse)

## ---------------------------
# Source Data from Water Quality Portal
#     List of all values:
#     https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&inline=true&group_cd=%
## ---------------------------
unique(metab_ts$site_name)
metab_ts$sitenumber <- sub("nwis_", "", metab_ts$site_name)
metab_ts$USGSnumber <- sub("nwis_", "USGS-", metab_ts$site_name)

str(metab_ts)
site_no<- unique(metab_ts$sitenumber)
USGSnumber<- unique(metab_ts$USGSnumber)

## General USGS site info:
Info <- readNWISsite(site_no) %>%
  select(agency_cd, site_no, station_nm, huc_cd, dec_lat_va, dec_long_va, alt_va, drain_area_va)
hist(Info$drain_area_va)

# Create NLDI feature lists for each gage
nldi_features <- lapply(site_no, function(gage) {
  list(featureSource = "nwissite", featureID = paste0("USGS-", gage))
})

# Get the NHDPlus COMID for each USGS gage and attach to Info data frame
Info$COMID <- sapply(nldi_features, function(feature) {
  discover_nhdplus_id(nldi_feature = feature)
})

## ---------------------------
# Bring in the NHD Attribute Data 
# downloaded at https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47
# National Land Cover Database (2016)  https://www.sciencebase.gov/catalog/item/5d66b3b6e4b0c4f70cefb11d
# reach catchments accumulated upstream proportional land cover through the river network

nlcd_acc <- read.csv("/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/NLCD16_ACC_CONUS.TXT")
head(nlcd_acc)
colnames(nlcd_acc) <- c("COMID","Water_Open_pct", "Ice_Snow_Perennial_pct", 
                         "Developed_OpenSpace_pct", "Developed_LowIntensity_pct", 
                         "Developed_MedIntensity_pct", "Developed_HiIntensity_pct", "Barren_Land_pct",
                         "Forest_Deciduous_pct", "Forest_Evergreen_pct", "Forest_mixed_pct",
                         "Shrub_Scrub_pct", "Grass_landHerbaceous_pct", "Grass_PastureHay_pct",
                         "Cultivated_Crops_pct", "Wetlands_Woody_pct", "Wetlands_EmergentHerb_pct",
                         "NoData")
# Subset NLCD by our sites of interest
nlcd_acc <- subset(nlcd_acc, nlcd_acc$COMID %in% Info$COMID) # the number of rows should be = to the number of rows in usgs_comid
# Create one df with all necessary info
usgs_nlcd <- merge(Info, nlcd_acc, by = "COMID", all = TRUE)
head(usgs_nlcd)


### https://www.sciencebase.gov/catalog/item/5703f6b5e4b0328dcb826d06
# reach catchments accumulated upstream proportional Generalized Geology Type attributes through the river network
nlcd_geoGGT <- read.csv("/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/BUSHREED_ACC_CONUS.txt")
head(nlcd_geoGGT)
colnames(nlcd_geoGGT) <- c("COMID", "geol_gneiss", "geol_granitic", "geol_ultramafic","geol_quarternary", "geol_sedimentary", 
                           "geol_volcanic", "geol_water", 
                           "geol_anorthositic", "geol_intermediate_pultonic",
                           "geol_NoData")
# Subset NLCD by our sites of interest
nlcd_geoGGT1 <- subset(nlcd_geoGGT, nlcd_geoGGT$COMID %in% usgs_nlcd$COMID)
nlcd_geoGGT2 <- merge(usgs_nlcd, nlcd_geoGGT1, by = "COMID", all = TRUE)


### https://www.sciencebase.gov/catalog/item/5703f6b5e4b0328dcb826d06
# reach catchments accumulated upstream proportional elevation through the river network
nlcd_basin <- read.csv("/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/BASIN_CHAR_CAT_CONUS.TXT")

# Subset NLCD by our sites of interest
nlcd_basin1 <- subset(nlcd_basin, nlcd_basin$COMID %in% nlcd_geoGGT1$COMID) # the number of rows should be = to the number of rows in usgs_comid
# Names 
# CAT_BASIN_AREA: flowline catchment area in square kilometers.
# CAT_BASIN_SLOPE: flowline catchment's average slope in percent
# CAT_ELEV_MEAN: flowline catchment's mean elevation in meters
# CAT_ELEV_MAX: flowline catchment's maximum elevation in meters
# CAT_STREAM_SLOPE: flowline's average slope in percent.
# CAT_STREAM_LENGTH: flowline's length in kilometers taken directly

# Change column names:
names(nlcd_basin1)
dat_NHDinfo <- merge(nlcd_geoGGT2, nlcd_basin1, by = "COMID", all = TRUE)


## Road density
# https://www.sciencebase.gov/catalog/file/get/57976a0ce4b021cadec97890?f=__disk__52%2F3f%2F23%2F523f238710502f61147307585d9b993fcb71bf48&transform=1&allowOpen=true
nlcd_road <- read.csv("/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/ROAD_DENS_CAT_CONUS.txt")%>%
  select(COMID, CAT_S1400, CAT_S1740, CAT_TOTAL_ROAD_DENS)
head(nlcd_road)
colnames(nlcd_road) <- c("COMID", "road_Rural_den","road_private_den", "road_Total_den")

# Subset NLCD by our sites of interest
nlcd_road1 <- subset(nlcd_road, nlcd_road$COMID %in% Info$COMID) # the number of rows should be = to the number of rows in usgs_comid
# Names 
# CAT_S1200 - Density is defined as the length of road divided by the catchment area,  These roads have one or more lanes of traffic in each direction, may or may not be divided, and usually have at-grade intersections with many other roads and driveways.
# CAT_S1400 - Density of local neighborhood roads divided by the catchment area, rural road- Generally a paved non-arterial street, road, or byway that usually has a single lane of traffic in each direction.
# CAT_S1740 - Density is defined as the length of road divided by the catchment area. A road within private property that is privately maintained for service, extractive, or other purposes
# CAT_TOTAL_ROAD_DENS - Density of all road types per NHDPlusV2 catchment. Density is defined as the length of road divided by the catchment area.
dat_NHDinfo1 <- merge(dat_NHDinfo, nlcd_road1, by = "COMID", all = TRUE)



## Canopy cover
# https://www.sciencebase.gov/catalog/item/570572e2e4b0d4e2b75718bc
nlcd_CNPY11 <- read.csv("/Users/kellyloria/Documents/UNR/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/CNPY11_BUFF100_CONUS.txt")%>%
  select(COMID, CAT_CNPY11_BUFF100, NODATA)
# CAT_CNPY11_BUFF100: Percent NLCD 2011 Tree Canopy in 100 meter riparian buffer per NHDPlus version 2 catchment
# NODATA: Percent of catchment that the source data does not cover. -9999 (flowline) denotes 100% NODATA.
head(nlcd_CNPY11)
colnames(nlcd_CNPY11) <- c("COMID", "Tree_Canopyin100mRip","NODATA")
# Subset NLCD by our sites of interest
nlcd_CNPY11 <- subset(nlcd_CNPY11, nlcd_CNPY11$COMID %in% Info$COMID) # 

dat_NHDinfo <- merge(dat_NHDinfo1, nlcd_CNPY11, by = "COMID", all = TRUE)

## ---------------------------
# Visualize

str(dat_NHDinfo)
names(dat_NHDinfo)

dat_NHDinfo1 <- dat_NHDinfo %>%
  mutate(state = sub(".* ([A-Z]{2})$", "\\1", station_nm))

dat_NHDinfo1 <- within(dat_NHDinfo1, {
  state[state == "BLACK EARTH CREEK NR BREWERY RD AT CROSS PLAINS,WI"] <- "WI"
  state[state == "ARKANSAS RIVER NEAR AVONDALE, CO."] <- "CO"
  state[state == "M FK BEARGRASS CR AT OLD CANNONS LN AT LOUISVILLE,"] <- "KY"
  state[state == "EAST CANYON CREEK AB EAST CYN RES NR MORGAN, UTAH"] <- "UT"
  state[state == "MUDDY FK AT MOCKINGBIRD VALLEY RD AT LOUISVILLE,KY"] <- "KY"
})

unique(dat_NHDinfo1$state)


## ---------------------------
# saveRDS(dat_NHDinfo1, "./data/Metab_metadata.rds")


# Filter relevant columns (Example: columns 8:47 as mentioned)
ws_smry <- dat_NHDinfo1 %>% select(state, COMID, 8:47)

# Principal component analysis (PCA) on watershed summary attributes
# Drop columns with >= 30% missing values, and rows with any missing values
pca_data <- ws_smry %>% 
  select(where(~ sum(is.na(.)) / length(.) < 0.3)) %>%  # Drop columns with >= 30% missing values
  drop_na()                                              # Drop rows with any missing values


domains <- pca_data$state

# Prepare data for PCA by removing non-numeric and non-relevant columns
pca_data <- pca_data %>% 
  select(-state, -COMID) %>% 
  select(where(~ sd(., na.rm = TRUE) != 0)) %>%  # Drop columns with no variance
  as.matrix()

# Categorize summary columns based on their initial patterns
col_patterns <- c("alt_va" = "elevation", "drain_area_va" = "drainage_area", "Water_Open_pct" = "water",
                  "Ice_Snow_Perennial_pct" = "ice_snow", "Developed_OpenSpace_pct" = "open_space",
                  "Developed_LowIntensity_pct" = "open_space", "Developed_MedIntensity_pct" = "developed",
                  "Developed_HiIntensity_pct" = "developed", "Barren_Land_pct" = "open_space",
                  "Forest_Deciduous_pct" = "forest", "Forest_Evergreen_pct" = "forest", 
                  "Forest_mixed_pct" = "forest", "Shrub_Scrub_pct" = "shrub", 
                  "Grass_landHerbaceous_pct" = "shrub", "Grass_PastureHay_pct" = "crops", 
                  "Cultivated_Crops_pct" = "crops", "Wetlands_Woody_pct" = "wetlands", 
                  "Wetlands_EmergentHerb_pct" = "wetlands", "geol_gneiss" = "geology_metamorphic", 
                  "geol_granitic" = "geology_igneous", "geol_ultramafic" = "geology_metamorphic", 
                  "geol_quarternary" = "geology_sedimentary", "geol_sedimentary" = "geology_sedimentary", 
                  "geol_volcanic" = "geology_igneous", "geol_intermediate_pultonic" = "geology_igneous", 
                  "CAT_BASIN_AREA" = "drainage_area", "CAT_BASIN_SLOPE" = "slope", 
                  "CAT_ELEV_MEAN" = "elevation", "CAT_ELEV_MIN" = "elevation", 
                  "CAT_ELEV_MAX" = "elevation", "CAT_STREAM_SLOPE" = "slope", 
                  "CAT_STREAM_LENGTH" = "stream_length", "road_Rural_den" = "road_density", 
                  "road_private_den" = "road_density", "road_Total_den" = "road_density", 
                  "Tree_Canopyin100mRip" = "shading")

# Apply the categorization to column names
smry_categories <- factor(sapply(colnames(pca_data), function(col) {
  col_patterns[col]
}))


pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# Visualize PCA eigenvalues
fviz_eig(pca)

# Visualize PCA biplot with variable categories
watershed_loadplot<- fviz_pca_biplot(pca, geom.var = 'arrow', geom.ind = 'point', title = '',
                col.var = smry_categories)

# Visualize PCA biplot with site domains
watershedplot <- fviz_pca_biplot(pca, geom.var = '', geom.ind = 'point', title = '',
                col.ind = as.factor(domains))

combined_plot <- ggarrange(watershed_loadplot, watershedplot, ncol = 2, nrow = 1)

# ggsave(plot = combined_plot, filename = paste("./figures/Watershed_attributes.png",sep=""),width=12,height=4.5,dpi=300)


