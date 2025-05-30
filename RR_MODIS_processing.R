library(StreamLight)
library(StreamLightUtils)


working_dir <- c("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Modis_data/")

sites <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/River_synch_cord.csv")


#Make a table for the MODIS request 
request_sites <- sites[, c("ID", "Latitude", "Longitude")] 

#Export your sites as a .csv for the AppEEARS request  
write.table(
  request_sites, 
  paste0(working_dir, "/RR_modis_sites.csv"), 
  sep = ",", 
  row.names = FALSE,
  quote = FALSE, 
  col.names = FALSE
)


MOD_unpack <- AppEEARS_unpack_QC(
  zip_file = "RR_stream_dat.zip", 
  zip_dir = working_dir, 
  request_sites[, "ID"]
)

# Let’s process the LAI data and visualize the results. The black line is the new fitted, interpolated, daily LAI.



MOD_processed <- AppEEARS_proc(
  unpacked_LAI = MOD_unpack,  
  fit_method = "Gu", 
  plot = TRUE,
  AppEEARS_proc =T
)


?AppEEARS_proc()

library(dplyr)
library(purrr)

## save lists as one dataframe: 
MOD_combined <- imap_dfr(MOD_processed, ~ mutate(.x, Site_ID = .y))


TS_plot <- ggplot(MOD_combined, aes(x = Date, y = LAI_proc, color = Site_ID)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_d(option = "viridis", name = "Water year") +
  theme_classic() 


saveRDS(MOD_combined, "/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Modis_data/StreamDat_MOD_combined.rds")


