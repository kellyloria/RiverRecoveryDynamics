############
## Data aggregation of all stream gauge and daily time sieries data
##
## Author: Kelly A. Loria
## Date Created: 2024-07-09
## Email: kellyloria @ gmail.com
##
## ---------------------------
#PC: setwd("R:/Users/kloria/Documents/River_Recovery_Analysis")
#Mac: setwd("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis")

## ---------------------------
## meta data for watershed/landscape differences: Metab_metadata.rds
##    complied in "Metab_watershed_metadata.R
## ---------------------------

## ---------------------------
library(dataRetrieval)
library(nhdplusTools)
library(scales)
library(tidyverse)
library(ggpubr)
## ---------------------------

## ---------------------------
## Load time series data:
metab_ts <- readRDS("./data/Metab_TS.rds")
str(metab_ts)

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

## --------------------------- 
## Explore different WQ products: 

pcode <- readNWISpCode("all")

infoCds <- pcode[grep("nitrate",
                      pcode$parameter_nm,
                      ignore.case = TRUE
), ]

## --------------------------- 
# Identify the minimum start date and maximum end date for each site

site_dates <- metab_ts %>%
  filter(!is.na(sitenumber)) %>%
  group_by(sitenumber) %>%
  summarize(start_date = min(as.Date(date), na.rm = TRUE),
            end_date = max(as.Date(date), na.rm = TRUE))

## --------------------------- 
# pH data pc: 00400
all_data <- list()

# Loop through each row in site_dates
for (i in 1:nrow(site_dates)) {
  # Extract the site number and date range for the current row
  siteNo <- site_dates$sitenumber[i]
  start.date <- site_dates$start_date[i]
  end.date <- site_dates$end_date[i]
  pCode <- "00400"  # pH parameter code
  
  # Print current values for debugging
  print(paste("Processing site:", siteNo, "Start date:", start.date, "End date:", end.date))
  
  # Check for NA values
  if (!is.na(siteNo) && !is.na(start.date) && !is.na(end.date)) {
    # Download the data using readNWISuv
    site_data <- tryCatch({
      readNWISdv(siteNumbers = siteNo, parameterCd = pCode, startDate = start.date, endDate = end.date)
    }, error = function(e) {
      message(paste("download_error", siteNo, ":", e$message))
      return(NULL)
    })
    
    # Append the downloaded data to the list
    if (!is.null(site_data)) {
      all_data[[length(all_data) + 1]] <- site_data
    }
  } else {
    message(paste("Skipping site", siteNo, "due to missing values in site number or date range"))
  }
}


# Define standard column names
standard_columns <- c("agency_cd", "site_no", "dateTime", "X_00400_00000", "tz_cd")

# Initialize an empty list to store the standardized data frames
standardized_data <- list()

# Standardize the columns in each data frame
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  if (!is.null(df) && "X_00400_00000" %in% colnames(df)) { # Ensure the data frame has pH data
    missing_cols <- setdiff(standard_columns, colnames(df))
    df[missing_cols] <- NA  # Add missing columns with NA values
    df <- df[standard_columns]  # Reorder columns to match the standard
    standardized_data[[i]] <- df
  }
}

# Remove NULL entries from the list
standardized_data <- Filter(Negate(is.null), standardized_data)

# Combine all the standardized data frames into one data frame
if (length(standardized_data) > 0) {
  combined_data <- do.call(rbind, standardized_data)
  
  # Remove the tz_cd column if it exists
  if ("tz_cd" %in% colnames(combined_data)) {
    combined_data <- combined_data %>%
      dplyr::select(-tz_cd)
  }
} else {
  combined_data <- NULL
}

# View the combined data
print(combined_data)

hist(combined_data$X_00400_00000)

combined_data_day <- combined_data%>%
  mutate(date= as.Date(dateTime)) %>%
  group_by(site_no, date) %>%
  summarise(pH= mean(X_00400_00000, na.rm=T))
  

pH_plot <- ggplot(data = combined_data_day, aes(x = date, y = pH, color=site_no)) +
  geom_point() + theme_bw() 

## ---------------------------
## rejoin metab and pH data 
 metab_ts_WQ <- metab_ts %>%
  left_join(combined_data_day, by = c("sitenumber" = "site_no", "date"))

## --------------------------- 
## turbidity pc: 63680

# Initialize an empty list to store the data
all_data <- list()

# Loop through each row in site_dates
for (i in 1:nrow(site_dates)) {
  # Extract the site number and date range for the current row
  siteNo <- site_dates$sitenumber[i]
  start.date <- site_dates$start_date[i]
  end.date <- site_dates$end_date[i]
  pCode <- "63680"  # pH parameter code # 80300, 82079, 01350, 61028, 63675, 63676, 72395, 72213, 
  # Print current values for debugging
  print(paste("Processing site:", siteNo, "Start date:", start.date, "End date:", end.date))
  
  # Check for NA values
  if (!is.na(siteNo) && !is.na(start.date) && !is.na(end.date)) {
    # Download the data using readNWISuv
    site_data <- tryCatch({
      readNWISdv(siteNumbers = siteNo, parameterCd = pCode, startDate = start.date, endDate = end.date)
    }, error = function(e) {
      message(paste("download_error", siteNo, ":", e$message))
      return(NULL)
    })
    
    # Append the downloaded data to the list
    if (!is.null(site_data)) {
      all_data[[length(all_data) + 1]] <- site_data
    }
  } else {
    message(paste("Skipping site", siteNo, "due to missing values in site number or date range"))
  }
}

########
# Define standard column names
standard_columns <- c("agency_cd", "site_no", "Date", "X_63680_00003")

# Initialize an empty list to store the standardized data frames
standardized_data <- list()

# Standardize the columns in each data frame
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  if (!is.null(df) && "X_63680_00003" %in% colnames(df)) { # Ensure the data frame has pH data
    missing_cols <- setdiff(standard_columns, colnames(df))
    df[missing_cols] <- NA  # Add missing columns with NA values
    df <- df[standard_columns]  # Reorder columns to match the standard
    standardized_data[[i]] <- df
  }
}

# Remove NULL entries from the list
standardized_data <- Filter(Negate(is.null), standardized_data)

# Combine all the standardized data frames into one data frame
if (length(standardized_data) > 0) {
  combined_data <- do.call(rbind, standardized_data)
  
  # Remove the tz_cd column if it exists
  if ("tz_cd" %in% colnames(combined_data)) {
    combined_data <- combined_data %>%
      dplyr::select(-tz_cd)
  }
} else {
  combined_data <- NULL
}

## ---------------------------
combined_data_FNU <- combined_data
# saveRDS(combined_data_FNU, "./data/FNU_TS_pcode63680.rds")
unique(combined_data_FNU$site_no) # 21 sites with FNU + metabolism 

## ---------------------------
## rejoin metab and FNU data 
metab_ts_WQ1 <- metab_ts_WQ %>%
  left_join(combined_data_FNU, by = c("sitenumber" = "site_no", "date"="Date")) %>%
  rename(turbidity_FNU="X_63680_00003")


## ---------------------------
## Quick visualizations: 
FNU_plot <- ggplot(data = metab_ts_WQ1%>%filter(turbidity_FNU>=0), aes(x = doy, y = log(turbidity_FNU+1), color=sitenumber)) +
  geom_point(size=0.75,  alpha=0.75) + theme_bw() +theme(legend.position = "bottom") + facet_grid(year~.)
#  ggsave(plot = FNU_plot, filename = paste("./figures/FNU_TS.png",sep=""),width=13,height=12,dpi=300)

GPP_FNU_plot <- ggplot(data = metab_ts_WQ1%>%filter(turbidity_FNU>=0), aes(y = GPP, x = log(turbidity_FNU+1), color=sitenumber), shape=as.factor(year)) +
  geom_point(size=0.75, alpha=0.75)+ theme_bw() +
  geom_hline(yintercept = 0, color = "gray50")

ER_FNU_plot <- ggplot(data = metab_ts_WQ1%>%filter(turbidity_FNU>=0), aes(y = c(ER*-1), x = log(turbidity_FNU+1), color=sitenumber), shape=as.factor(year)) +
  geom_point(size=0.75, alpha=0.75) + theme_bw() +
  geom_hline(yintercept = 0, color = "gray50")

FNU_grid <- ggarrange(GPP_FNU_plot, 
                      ER_FNU_plot,
                     ncol = 2, nrow = 1,
                     common.legend = TRUE, 
                     labels=c("A", "B"),
                     legend = "bottom")

# ggsave(plot = FNU_grid, filename = paste("./figures/FNU_metab_TS.png",sep=""),width=9,height=7,dpi=300)


## --------------------------- 
## Chem? nitrogen 00600 NO: 00904

# Initialize an empty list to store the data
all_data <- list()

# Loop through each row in site_dates
for (i in 1:nrow(site_dates)) {
  # Extract the site number and date range for the current row
  siteNo <- site_dates$sitenumber[i]
  start.date <- site_dates$start_date[i]
  end.date <- site_dates$end_date[i]
  pCode <- c("00600", "00605", "00608", "00613")  # pH parameter code # 00600 , 00601
  # Print current values for debugging
  print(paste("Processing site:", siteNo, "Start date:", start.date, "End date:", end.date))
  
  # Check for NA values
  if (!is.na(siteNo) && !is.na(start.date) && !is.na(end.date)) {
    # Download the data using readNWISuv
    site_data <- tryCatch({
      readWQPqw(siteNumbers = siteNo, parameterCd = pCode, 
                  startDate = start.date, endDate = end.date)
    }, error = function(e) {
      message(paste("download_error", siteNo, ":", e$message))
      return(NULL)
    })
    
    # Append the downloaded data to the list
    if (!is.null(site_data)) {
      all_data[[length(all_data) + 1]] <- site_data
    }
  } else {
    message(paste("Skipping site", siteNo, "due to missing values in site number or date range"))
  }
}





########
# Define standard column names
standard_columns <- c("agency_cd", "site_no", "Date", "X_63680_00003")

# Initialize an empty list to store the standardized data frames
standardized_data <- list()

# Standardize the columns in each data frame
for (i in seq_along(all_data)) {
  df <- all_data[[i]]
  if (!is.null(df) && "X_63680_00003" %in% colnames(df)) { # Ensure the data frame has pH data
    missing_cols <- setdiff(standard_columns, colnames(df))
    df[missing_cols] <- NA  # Add missing columns with NA values
    df <- df[standard_columns]  # Reorder columns to match the standard
    standardized_data[[i]] <- df
  }
}

# Remove NULL entries from the list
standardized_data <- Filter(Negate(is.null), standardized_data)

# Combine all the standardized data frames into one data frame
if (length(standardized_data) > 0) {
  combined_data <- do.call(rbind, standardized_data)
  
  # Remove the tz_cd column if it exists
  if ("tz_cd" %in% colnames(combined_data)) {
    combined_data <- combined_data %>%
      dplyr::select(-tz_cd)
  }
} else {
  combined_data <- NULL
}


Phosphorus <- readWQPdata(
  statecode = siteNo, 
  characteristicName = "Phosphorus",
  startDateLo = start.date,
  ignore_attributes = TRUE,
  convertType = FALSE
)

data5 <- readNWISgwl("263819081585801", parameterCd = "72019")

