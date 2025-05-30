############
## Meta-Data aggregation of all watershed data
##
## Author: Kelly A. Loria
## Date Created: 2024-06-26
## Email: kellyloria @ gmail.com
##
## ---------------------------
#PC: setwd("R:/Users/kloria/Documents/River_Recovery_Analysis")
#Mac: setwd("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis")

USGS_sites <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/River_synch_cord.csv") %>%
  dplyr::select(ID, Latitude, Longitude)

USGS_sites$sitenumber <- sub("nwis_", "", USGS_sites$ID)
USGS_sites$ID <- sprintf("%08d", as.numeric(USGS_sites$ID))

###########
# Download watershed attributes from each site: 
## ---------------------------
library(dataRetrieval)
library(nhdplusTools)
library(scales)
library(tidyverse)

library(nhdplusTools)


## General USGS site info:
Info <- readNWISsite(USGS_sites$ID) #%>% dplyr::select(sitenumber, Latitude, Longitude)
hist(Info$drain_area_va)

# Create NLDI feature lists for each gage
nldi_features <- lapply(Info$site_no, function(gage) {
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

nlcd_acc <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/NLCD16_ACC_CONUS.TXT")
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
nlcd_geoGGT <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/BUSHREED_ACC_CONUS.txt")
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
nlcd_basin <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/BASIN_CHAR_CAT_CONUS.TXT")

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
nlcd_road <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/ROAD_DENS_CAT_CONUS.txt")%>%
 dplyr::select(COMID, CAT_S1400, CAT_S1740, CAT_TOTAL_ROAD_DENS)
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
nlcd_CNPY11 <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/CNPY11_BUFF100_CONUS.txt") %>%
  dplyr::select(COMID, CAT_CNPY11_BUFF100, NODATA)
# CAT_CNPY11_BUFF100: Percent NLCD 2011 Tree Canopy in 100 meter riparian buffer per NHDPlus version 2 catchment
# NODATA: Percent of catchment that the source data does not cover. -9999 (flowline) denotes 100% NODATA.
head(nlcd_CNPY11)
colnames(nlcd_CNPY11) <- c("COMID", "Tree_Canopyin100mRip","NODATA")
# Subset NLCD by our sites of interest
nlcd_CNPY11 <- subset(nlcd_CNPY11, nlcd_CNPY11$COMID %in% Info$COMID) # 

dat_NHDinfo <- merge(dat_NHDinfo1, nlcd_CNPY11, by = "COMID", all = TRUE)

#### LOOK up hand variable 

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
  state[state == "Little Blue R. at Lees Summit Rd in Independence"] <- "MO"
  state[state == "E. Fk. Black R. bl Lower Taum Sauk Reservoir"] <- "TX"
})

unique(dat_NHDinfo1$state)


####
####

COMID <- Info$COMID

library(nhdplusTools)

### Get stream order:
flowline_attrs <- subset_nhdplus(
  comids = COMID,
  output_file = NULL,
  return_data = TRUE,
  nhdplus_data = "download"
)

df <- flowline_attrs$NHDFlowline_Network

filtered_df <- df %>% 
  filter(comid %in% COMID)

filtered_df1 <- sf::st_drop_geometry(filtered_df)

filtered_df1 <- as.data.frame(filtered_df1)

library(dplyr)

# Left join (all rows from filtered_df1, matching from dat_NHDinfo1)
joined_df <- left_join(filtered_df1, dat_NHDinfo1, by = c("comid" = "COMID"))





## ---------------------------
# saveRDS(joined_df, "/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Metab_metadata_25v2.rds")

dat_NHDinfo1<- readRDS("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Metab_metadata_25v2.rds")

nlcd_BEDPERM <- read.csv("/Users/kellyloria/Documents/UNR/Course\ work/Fall2020Projects/NHD_Tools/NHD_AttributeFiles/BEDPERM_ACC_CONUS.TXT")
head(nlcd_BEDPERM)

# ACC_BEDPERM_1: Percent of all upstream NHDPlus version 2 flowline catchments whose bedrock permeability class is not a principal aquifer
# ACC_BEDPERM_6: Percent of all upstream NHDPlus version 2 flowline catchments whose bedrock permeability class is Unconsolidated sand and gravel, based on divergence routing.
colnames(nlcd_BEDPERM)[2] <- "PerNotAquifer"
colnames(nlcd_BEDPERM)[7] <- "PerSandGravel"

nlcd_BEDPERM1<- nlcd_BEDPERM%>%dplyr::select(COMID, PerNotAquifer, PerSandGravel)

# Subset NLCD by our sites of interest
nlcd_BEDPERM1 <- subset(nlcd_BEDPERM1, nlcd_BEDPERM1$COMID %in% nlcd_BEDPERM1$COMID)
dat_NHDinfo1 <- merge(dat_NHDinfo1, nlcd_BEDPERM1, by = "COMID", all = TRUE)






# Filter relevant columns (Example: columns 8:47 as mentioned)
ws_smry <- dat_NHDinfo1 %>% dplyr::select(state, COMID, 44:81)


ws_smry1 <- dat_NHDinfo1 %>% dplyr::select(state, COMID, 49:58, 65, 71, 72, 79,81)



summary(ws_smry)

domains <- ws_smry$state

# Prepare data for PCA by removing non-numeric and non-relevant columns
pca_data <- ws_smry1 %>% 
  dplyr::select(-state, -COMID) %>% 
  dplyr::select(where(~ sd(., na.rm = TRUE) != 0)) %>%  # Drop columns with no variance
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
screeplot(pca)

library(factoextra)

# Visualize PCA eigenvalues
scree_plt <- fviz_eig(pca)

# ggsave(plot = scree_plt, filename = paste("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/figures/Watershed_attributes_screeplt_25.png",sep=""),width=4,height=3.5,dpi=300)


biplot(prcomp(pca_data, center = TRUE, scale. = TRUE))


# Visualize PCA biplot with variable categories
watershed_loadplot<- fviz_pca_biplot(pca, geom.var = 'arrow', geom.ind = 'point', title = '',
                col.var = smry_categories)

# Visualize PCA biplot with site domains
watershedplot <- fviz_pca_biplot(pca, geom.var = '', geom.ind = 'point', title = '',
                col.ind = as.factor(domains))

combined_plot <- ggarrange(watershed_loadplot, watershedplot, ncol = 2, nrow = 1)

# ggsave(plot = combined_plot, filename = paste("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/figures/Watershed_attributes_25.png",sep=""),width=12,height=4.5,dpi=300)





#####



# Filter relevant columns (Example: columns 8:47 as mentioned)
ws_smry <- dat_NHDinfo1 %>% dplyr::select(state, COMID, 46:53, 56:57, 78:80)

str(ws_smry)


ws_smry1 <- dat_NHDinfo1 %>% dplyr::select(state, COMID, 49:58, 65, 71, 72, 79,81)



summary(ws_smry)

domains <- ws_smry$state

# Prepare data for PCA by removing non-numeric and non-relevant columns
pca_data <- ws_smry %>% 
  dplyr::select(-state, -COMID) %>% 
  dplyr::select(where(~ sd(., na.rm = TRUE) != 0)) %>%  # Drop columns with no variance
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
screeplot(pca)

library(factoextra)

# Visualize PCA eigenvalues
scree_plt <- fviz_eig(pca)

# ggsave(plot = scree_plt, filename = paste("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/figures/Watershed_attributes_screeplt_25.png",sep=""),width=4,height=3.5,dpi=300)


biplot(prcomp(pca_data, center = TRUE, scale. = TRUE))


# Visualize PCA biplot with variable categories
watershed_loadplot<- fviz_pca_biplot(pca, geom.var = 'arrow', geom.ind = 'point', title = '',
                                     col.var = smry_categories)

# Visualize PCA biplot with site domains
watershedplot <- fviz_pca_biplot(pca, geom.var = '', geom.ind = 'point', title = '',
                                 col.ind = as.factor(domains))

combined_plot <- ggarrange(watershed_loadplot, watershedplot, ncol = 2, nrow = 1)

###
library(randomForest)

# Prepare data: assuming ws_smry contains 'synchrony' column
rf_data <- ws_smry[, c("Developed_OpenSpace_pct", "Developed_LowIntensity_pct",
                       "Developed_MedIntensity_pct", "Developed_HiIntensity_pct",
                       "Barren_Land_pct", "Forest_Deciduous_pct",
                       "Forest_Evergreen_pct", "Forest_mixed_pct",
                       "Grass_PastureHay_pct", "Cultivated_Crops_pct",
                       "road_Rural_den", "road_private_den", "road_Total_den")]

rf_model <- randomForest(synchrony ~ ., data = data.frame(synchrony = ws_smry$synchrony, rf_data),
                         importance = TRUE, ntree = 1000)

# Plot importance
varImpPlot(rf_model)


