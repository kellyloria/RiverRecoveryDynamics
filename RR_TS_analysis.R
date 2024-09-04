############
## Time seires analysis aggregation of all stream gauge and daily time sieries data
##
## Author: Kelly A. Loria
## Date Created: 2024-07-29
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

## ---------------------------
## Start evaluating flow characteristics


# scale flow by watershed area:
dat_ts$site_no <- sub("nwis_", "", dat_ts$site_name)
dat_ts$USGSnumber <- sub("nwis_", "USGS-", dat_ts$site_name)
site_no<- unique(dat_ts$sitenumber)
USGSnumber<- unique(dat_ts$USGSnumber)

## General USGS site info:
Info <- readNWISsite(site_no) %>%
  select(agency_cd, site_no, station_nm, huc_cd, dec_lat_va, dec_long_va, alt_va, drain_area_va)
hist(Info$drain_area_va)
hist(Infoq$drain_area_va)

Infoq<- Info%>%
  select(site_no, drain_area_va) %>%
  mutate(drain_area_va = c(2.58999*drain_area_va))


dat_ts <- dat_ts%>%
  left_join(Infoq, by = c("site_no")) %>%
  mutate(scale_Q = c(Q/drain_area_va)) # %>%
#select(-drain_area_va)
hist(dat_ts$Q)
hist(dat_ts$scale_Q)
hist(log(dat_ts$scale_Q+1))

## ---------------------------
### A. Flow Duration Analysis: Estimating flow exceedance probabilities
result <- dat_ts %>%
  dplyr::group_by(site_no) %>%
  dplyr::summarize(
    total_intervals = n(),
    ranked_scale_discharge = list(sort(scale_Q, decreasing = TRUE)),
    ranked_discharge = list(sort(Q, decreasing = TRUE))
  ) %>%
  mutate(
    exceedence_prob = purrr::map(ranked_discharge, ~ seq(1, length(.), length.out = length(.)) / length(.) * 100)
  ) %>%
  unnest(cols = c(ranked_scale_discharge, ranked_discharge, exceedence_prob))

# Plotting for each site
POE_plot <- result %>%
  ggplot(aes(y = ranked_scale_discharge, x = exceedence_prob, color = site_no)) +
  geom_line(linewidth=1.5) + theme_bw() + theme(legend.position = "bottom") +
  #scale_color_manual(values=alpha(c("#3283a8", "#a67d17"),0.5)) +
  labs(y = "Discharge (cms/km)", x = "Exceedence probability",
       title = "Exceedence Probability vs. Discharge by Site") #+ facet_grid(.~Site)
POE_plot

POE10 <- result %>%
  filter(exceedence_prob>9.5 & exceedence_prob <10.5) %>%
  group_by(site_no)%>%
  summarise(
    POE10_Q= mean(ranked_scale_discharge, na.rm =T))

dat_ts1 <- dat_ts%>%
  left_join(POE10, by = c("site_no")) 


POE_plot <- dat_ts1 %>%
  ggplot(aes(y = POE10_Q, x = log(drain_area_va) , color = site_no)) +
  geom_point(size=1) + theme_bw() + theme(legend.position = "bottom") +
  geom_smooth(method="lm", color="grey40", se=F) +
  labs(x = "log(Drainage area) (km)", y =  "Q at 10% POE (m^3 s^-1 km^-1)") #
POE_plot

# ggsave(plot = POE_plot, filename = paste("./figures/10POE_DrainageArea.png",sep=""),width=5,height=7,dpi=300)
### quick stat 
# flow and drainage area
POE_mod <- glm(POE10_Q~log(drain_area_va), data=dat_ts1)
summary(POE_mod)
hist(residuals(POE_mod))

###

POE5 <- result %>%
  filter(exceedence_prob>4.5 & exceedence_prob <5.5) %>%
  group_by(site_no)%>%
  summarise(
    POE5_Q= mean(ranked_scale_discharge, na.rm =T))

dat_ts <- dat_ts%>%
  left_join(POE5, by = c("site_no")) 


POE_plot <- dat_ts %>%
  ggplot(aes(y = POE5_Q, x = log(drain_area_va) , color = site_no)) +
  geom_point(size=1) + theme_bw() + theme(legend.position = "bottom") +
  geom_smooth(method="lm", color="grey40", se=F) +
  labs(x = "log(Drainage area) (km)", y =  "Q at 5% POE (m^3 s^-1 km^-1)") #
POE_plot

# ggsave(plot = POE_plot, filename = paste("./figures/5POE_DrainageArea.png",sep=""),width=7,height=7.5,dpi=300)
### quick stat 
# flow and drainage area
POE_mod <- glm(POE5_Q~log(drain_area_va), data=dat_ts)
summary(POE_mod)
hist(residuals(POE_mod))


## ---------------------------
### B. Create fxn to identify when a site experiences POE of 5% flow
##
# 
# # Assuming 'dat_ts' has columns: site_no, date, and scaled_Q
# dat_ts1 <- dat_ts %>%
#   arrange(site_no, date) %>% # Ensure data is ordered by site and date
#   group_by(site_no) %>%
#   mutate(
#     # Create a lagged column for scaled_Q
#     prev_scaled_Q = lag(scale_Q),
#     # Initialize event column
#     event = NA
#   )
# # Loop through each site and assign events
# for (site in unique(dat_ts1$site_no)) {
#   dat_ts_site <- dat_ts1 %>% filter(site_no == site)
#   
#   for (i in 2:nrow(dat_ts_site)) {
#     if (dat_ts_site$scale_Q[i] >= 3 * dat_ts_site$prev_scaled_Q[i]) {
#       dat_ts_site$event[i] <- "Q3_event"
#     } else {
#       dat_ts_site$event[i] <- "none"
#     }
#   }
#   
#   # Update the main data frame with events for the current site
#   dat_ts1 <- dat_ts1 %>% 
#     filter(site_no != site) %>%
#     bind_rows(dat_ts_site)
# }
# 
# # Ensure no NA values in 'event' column
# dat_ts1 <- dat_ts1 %>%
#   mutate(event = ifelse(is.na(event), "none", event))
# # Print the result to check
# print(dat_ts1)
# 
# 
# 
# POE_plot <- dat_ts1 %>%
#   filter(!is.na(turbidity_FNU))%>%
#   ggplot(aes(x = date, y = GPP, color = site_no, shape=event)) +
#   geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
#   scale_shape_manual(values=c(19,3)) +
#   facet_grid(site_no~.)
#   geom_smooth(method="lm", color="grey40", se=F) +
#   labs(x = "log(Drainage area) (km)", y =  "Q at 5% POE (m^3 s^-1 km^-1)") 
# POE_plot


## POE EVENT?
# Assuming 'dat_ts' has columns: site_no, date, and scaled_Q
dat_ts1 <- dat_ts %>%
  arrange(site_no, date) %>% # Ensure data is ordered by site and date
  group_by(site_no) %>%
  mutate(
    # Initialize event column
    POE_event = NA
  )
# Loop through each site and assign events
for (site in unique(dat_ts1$site_no)) {
  dat_ts_site <- dat_ts1 %>% filter(site_no == site)
  
  for (i in 2:nrow(dat_ts_site)) {
    if (dat_ts_site$scale_Q[i] >= dat_ts_site$POE5_Q[i]) {
      dat_ts_site$POE_event[i] <- "POE_5"
    } else {
      dat_ts_site$POE_event[i] <- "none"
    }
  }
  
  # Update the main data frame with events for the current site
  dat_ts1 <- dat_ts1 %>% 
    filter(site_no != site) %>%
    bind_rows(dat_ts_site)
}

dat_ts1 <- dat_ts1 %>%
  mutate(POE_event = ifelse(is.na(POE_event), "none", POE_event))


POE_plot <- dat_ts1 %>%
  filter(!is.na(turbidity_FNU))%>%
  ggplot(aes(x = date, y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3))    +
  facet_grid(site_no~.)


POE_plot <- dat_ts1 %>%
  filter(!is.na(turbidity_FNU))%>%
  ggplot(aes(x = log(scale_Q+1), y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) 

POE_plot <- dat_ts1 %>%
  filter(!is.na(turbidity_FNU))%>%
  ggplot(aes(x = log(turbidity_FNU+1), y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) 

POE_plot <- dat_ts1 %>%
  filter(!is.na(turbidity_FNU))%>%
  ggplot(aes(x = log(turbidity_FNU+1), y = c(ER*-1), color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) 


#########
###########
dat_tsq<- dat_ts1%>% filter(GPP>0 & ER<0 &!order_group=="unknown")


# start with general glms of GPP and ER to pH FNU and spc 


chart.Correlation(dat_tsq[c(3,33,35,36,39)])
hist(dat_tsq$GPP)
hist(log(dat_tsq$GPP+1))

GPP_mod <-lmer(GPP~ scale(pH) + 
                   scale(turbidity_FNU) + 
                   scale(SPC) + 
                   scale(scale_Q) + order_group +
                  (1|site_name), 
                  #(1|order_group), 
                    data=dat_tsq)
summary(GPP_mod)
hist(residuals(GPP_mod))
r.squaredGLMM(GPP_mod)
vif(GPP_mod)


GPP_ph_plot <- dat_tsq %>%
  ggplot(aes(x = pH, y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

GPP_FNU_plot <- dat_tsq %>%
  ggplot(aes(x = log(turbidity_FNU+1), y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

GPP_SPC_plot <- dat_tsq %>%
  ggplot(aes(x = log(SPC+1), y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

GPP_Q_plot <- dat_tsq %>%
  ggplot(aes(x = scale_Q, y = GPP, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)


GPP_grid <- ggarrange(GPP_ph_plot, 
                      GPP_FNU_plot,
                      GPP_SPC_plot,
                      GPP_Q_plot,
                      ncol = 4, nrow = 1,
                      common.legend = TRUE, 
                     # labels=c("A", "B"),
                      legend = "bottom")

# ggsave(plot = GPP_grid, filename = paste("./figures/GPP_glm.png",sep=""),width=10,height=11,dpi=300)




###
###
# ER
hist(dat_tsq$ER)
dat_tsq$ER_abs <- c(dat_tsq$ER*-1)
chart.Correlation(dat_tsq[c(42,33,35,36,39)])
hist(dat_tsq$ER_abs)
hist(log(dat_tsq$ER_abs+1))

ER_mod <-lmer(ER_abs~ scale(pH) + 
                 scale(turbidity_FNU) + 
                 scale(SPC) + 
                 scale(scale_Q) +
                 (1|site_name), 
               #(1|order_group), 
               data=dat_tsq)
summary(ER_mod)
hist(residuals(ER_mod))
r.squaredGLMM(ER_mod)
vif(ER_mod)


hist(residuals(POE_mod))


ER_ph_plot <- dat_tsq %>%
  ggplot(aes(x = pH, y = ER_abs, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

ER_FNU_plot <- dat_tsq %>%
  ggplot(aes(x = log(turbidity_FNU+1), y = ER_abs, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

ER_SPC_plot <- dat_tsq %>%
  ggplot(aes(x = log(SPC+1), y = ER_abs, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

ER_Q_plot <- dat_tsq %>%
  ggplot(aes(x = scale_Q, y = ER_abs, color = site_no, shape=POE_event )) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)


ER_grid <- ggarrange(ER_ph_plot, 
                      ER_FNU_plot,
                      ER_SPC_plot,
                      ER_Q_plot,
                      ncol = 4, nrow = 1,
                      common.legend = TRUE, 
                      # labels=c("A", "B"),
                      legend = "bottom")

# ggsave(plot = ER_grid, filename = paste("./figures/ER_glm.png",sep=""),width=10,height=11,dpi=300)

## storm even identification :
# Install required packages if not already installed
# install.packages("hydrostats")
# install.packages("zoo")

##############
###############
## maybe get rid of nwis_14182500, nwis_10133980, nwis_03597860*, nwis_03067510 
# Initialize an empty list to store results for each site

library(dplyr)
library(lubridate)

# Function to calculate baseflow separation
BaseflowSeparation <- function(streamflow, filter_parameter = 0.925, passes = 3) {
  suppressWarnings(Ends <- c(1, length(streamflow)) * rep(1, (passes + 1)))
  suppressWarnings(AddToStart <- c(1, -1) * rep(1, passes))
  btP <- streamflow  # Previous pass's baseflow approximation
  qft <- vector(length = length(streamflow))
  bt <- vector(length = length(streamflow))
  bt[1] <- if (streamflow[1] < quantile(streamflow, 0.25)) streamflow[1] else mean(streamflow) / 1.5
  
  for (j in 1:passes) {
    for (i in (Ends[j] + AddToStart[j]):Ends[j + 1]) {
      if ((filter_parameter * bt[i - AddToStart[j]] + ((1 - filter_parameter) / 2) * (btP[i] + btP[i - AddToStart[j]])) > btP[i]) {
        bt[i] <- btP[i]
      } else {
        bt[i] <- filter_parameter * bt[i - AddToStart[j]] + ((1 - filter_parameter) / 2) * (btP[i] + btP[i - AddToStart[j]])
      }
      qft[i] <- streamflow[i] - bt[i]
    }
    if (j < passes) {
      btP <- bt
      bt[Ends[j + 1]] <- if (streamflow[Ends[j + 1]] < mean(btP)) streamflow[Ends[j + 1]] / 1.2 else mean(btP)
    }
  }
  f <- data.frame(bt, qft)
  return(f)
}

# Function to calculate water year
calculate_water_year <- function(date) {
  year <- year(date)
  if (month(date) >= 10) {
    return(year + 1)
  } else {
    return(year)
  }
}

# Initialize an empty list to store results for each site
all_sites_results <- list()

unique_sites <- unique(dat_tsq$site_name)

# Loop through each unique site name
for (site in unique_sites) {
  # Filter data for the current site
  discharge_data <- dat_tsq %>%
    filter(site_name == site) %>%
    dplyr::select(site_name, date, scale_Q)
  
  # Apply the baseflow separation function
  bfs <- BaseflowSeparation(discharge_data$scale_Q, passes = 3)
  bfs[, 1] <- as.numeric(bfs[, 1])
  
  # Combine the baseflow separation results with the original data
  flow_data <- cbind(discharge_data, bfs)
  
  # Add water year column
  flow_data$WaterYear <- sapply(flow_data$date, calculate_water_year)
  flow_data$DOY <- yday(flow_data$date)
  
  # Add water year day of year
  flow_data <- flow_data %>%
    mutate(wy_DOY = ifelse(month(date) >= 10, yday(date) - yday(as.Date(paste0(year(date), "-10-01"))) + 1,
                           yday(date) + (365 - yday(as.Date(paste0(year(date), "-09-30"))))))
  
  # Order data by site_name, WaterYear, and wy_DOY
  flow_data <- flow_data %>%
    arrange(site_name, WaterYear, wy_DOY)
  
  # Identify putative flow events as 10% above baseflow
  flow_data$Event <- ifelse(flow_data$scale_Q > 1.10 * flow_data$bt, 1, 0)
  # Identify putative flow events as 75% above baseflow
  flow_data$Event_LG <- ifelse(flow_data$scale_Q > 1.75 * flow_data$bt, 1, 0)
  
  # Differentiate events within each water year
  flow_data <- flow_data %>%
    group_by(site_name, WaterYear) %>%
    mutate(EventID = cumsum(Event == 1 & (lag(Event, default = 0) == 0))) %>%
    ungroup()
  
  # Calculate event length
  flow_data <- flow_data %>%
    group_by(site_name, WaterYear, EventID) %>%
    mutate(Event_length = sum(Event)) %>%
    ungroup()
  
  # Store the results for the current site
  all_sites_results[[site]] <- flow_data
}

# Combine results for all sites into one data frame
combined_results <- bind_rows(all_sites_results)

# Ensure that the combined_results is ordered correctly
combined_results <- combined_results %>%
  arrange(site_name, WaterYear, wy_DOY)

# Display the structure of the combined results
str(combined_results)

##############
combined_dat1 <- dat_tsq %>%
  full_join(combined_results, by = c("site_no", "site_name","date", "scale_Q"))

# saveRDS(combined_dat1, "./data/Flowevents_metab_TS.rds")


## Quick check 
dat_check <-combined_dat %>%
  filter(site_name== "nwis_09406000")

dat_check2 <-combined_dat %>%
  filter(site_name!= "nwis_14182500" & site_name!= "nwis_10133980" & 
           site_name!= "nwis_03597860" & site_name!= "nwis_03067510")


## make a quick box plots 

Q_plot <- dat_check2 %>%
  ggplot(aes(y = log(turbidity_FNU+1), x = as.factor(Event_LG), color = site_no, fill= site_no)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "right" ) + 
  facet_grid(.~order_group)  


Q_plot <- dat_check2 %>%
  ggplot(aes(y = log(SPC+1), x = as.factor(Event_LG), color = site_no, fill= site_no)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "right" ) + 
  facet_grid(.~order_group) 


Q_plot <- dat_check2 %>%
  ggplot(aes(y = pH, x = as.factor(Event), color = site_no, fill= site_no)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "right" ) + 
  facet_grid(.~order_group) 

Q_plot <- dat_check2 %>%
  ggplot(aes(y = GPP, x = as.factor(Event), color = site_no, fill= site_no)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "right" ) + 
  facet_grid(.~order_group) 

Q_plot <- dat_check2 %>%
  ggplot(aes(y = ER, x = as.factor(Event), color = site_no, fill= site_no)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "right" ) + 
  facet_grid(.~order_group) 


#####

###KINDA Broken here:
str(combined_dat)

df <- as.data.frame(combined_dat1)
str(df)


library(dplyr)

# Function to calculate the average rate before the flood event
calculate_pre_flood_rate <- function(df, rate_col, event_col, bt_col) {
  df %>%
    filter(!is.na(!!sym(rate_col)) & !!sym(bt_col) == 1 & !!sym(event_col) == 0) %>%
    summarise(pre_flood_rate = mean(!!sym(rate_col), na.rm = TRUE)) %>%
    pull(pre_flood_rate)
}

########
### BUG ##

# Function to calculate the average rate before the flood event
calculate_pre_flood_rate <- function(df, rate_col, event_col, bt_col) {
  # Print the initial state of the dataframe for debugging
  print(paste("Initial data for", rate_col))
  print(head(df))
  
  # Filter the data
  filtered_df <- df %>%
    filter(!is.na(!!sym(rate_col)) & !!sym(bt_col) == 1 & !!sym(event_col) == 0)
  
  # Print the filtered data for debugging
  print(paste("Filtered data for", rate_col))
  print(head(filtered_df))
  
  # Check if the filtered dataframe is empty
  if (nrow(filtered_df) == 0) {
    return(NA)
  }
  
  # Calculate the mean pre-flood rate
  pre_flood_rate <- filtered_df %>%
    summarise(pre_flood_rate = mean(!!sym(rate_col), na.rm = TRUE)) %>%
    pull(pre_flood_rate)
  
  # Return the calculated pre-flood rate
  return(pre_flood_rate)
}

# Example usage with debug prints
pre_flood_GPP <- calculate_pre_flood_rate(event_data, "GPP", "Event", "bt")
pre_flood_ER <- calculate_pre_flood_rate(event_data, "ER", "Event", "bt")

# Check the results
print(paste("Pre-flood GPP:", pre_flood_GPP))
print(paste("Pre-flood ER:", pre_flood_ER))




#########
##########











# Function to get the first modelable day after the flood event
get_first_modelable_day <- function(df, date_col, event_col) {
  df %>%
    filter(!!sym(event_col) == 1) %>%
    slice(1) %>%
    pull(!!sym(date_col))
}

# Calculate reduction percentage and response ratios
calculate_metrics <- function(df, rate_col, pre_flood_rate) {
  df %>%
    mutate(
      Reduction_percent = 1 - (!!sym(rate_col) / pre_flood_rate) * 100,
      RateRR = !!sym(rate_col) / pre_flood_rate
    )
}

# Identify unique flood events by site and WaterYear
unique_events <- df %>%
  filter(Event == 1) %>%
  select(site_name, WaterYear, EventID) %>%
  distinct()

# Initialize an empty list to store results for each event
event_results <- list()

# Loop through each unique flood event
for (i in 1:nrow(unique_events)) {
  event_info <- unique_events[i, ]
  
  # Filter data for the current event
  event_data <- df %>%
    filter(site_name == event_info$site_name &
             WaterYear == event_info$WaterYear &
             EventID == event_info$EventID)
  
  # Calculate pre-flood GPP and ER rates
  pre_flood_GPP <- calculate_pre_flood_rate(event_data, "GPP", "Event", "bt")
  pre_flood_ER <- calculate_pre_flood_rate(event_data, "ER", "Event", "bt")
  
  # Get the first modelable day after the flood event
  first_modelable_day <- get_first_modelable_day(event_data, "date", "Event")
  
  # Filter data for the post-flood period
  post_flood_data <- df %>%
    filter(site_name == event_info$site_name &
             WaterYear == event_info$WaterYear &
             date >= first_modelable_day)
  
  # Calculate reduction percentage and response ratios for GPP and ER
  post_flood_metrics <- post_flood_data %>%
    calculate_metrics("GPP", pre_flood_GPP) %>%
    calculate_metrics("ER", pre_flood_ER)
  
  # Store the results for the current event
  event_results[[i]] <- post_flood_metrics
}

# Combine results for all events into one data frame
combined_event_results <- bind_rows(event_results)



########3
library(dplyr)
library(lubridate)

# Function to calculate the average rate before the flood event
calculate_pre_flood_rate <- function(df, rate_col, event_col, bt_col) {
  pre_flood_df <- df %>%
    filter(!!sym(bt_col) == 1 & !!sym(event_col) == 0)
  
  if (nrow(pre_flood_df) > 0) {
    return(mean(pre_flood_df[[rate_col]], na.rm = TRUE))
  } else {
    return(NA)
  }
}

# Function to get the first modelable day after the flood event
get_first_modelable_day <- function(df, date_col, event_col) {
  post_flood_df <- df %>%
    filter(!!sym(event_col) == 1) %>%
    slice(1)
  
  if (nrow(post_flood_df) > 0) {
    return(post_flood_df[[date_col]])
  } else {
    return(NA)
  }
}

# Initialize an empty list to store results for each site
all_sites_results <- list()

# Identify unique flood events by site and WaterYear
unique_events <- df %>%
  filter(Event == 1) %>%
  select(site_name, WaterYear, EventID) %>%
  distinct()

# Loop through each unique flood event
for (i in 1:nrow(unique_events)) {
  event_info <- unique_events[i, ]
  
  # Filter data for the current event
  event_data <- df %>%
    filter(site_name == event_info$site_name &
             WaterYear == event_info$WaterYear &
             EventID == event_info$EventID)
  
  # Calculate pre-flood GPP and ER rates
  pre_flood_GPP <- calculate_pre_flood_rate(event_data, "GPP", "Event", "bt")
  pre_flood_ER <- calculate_pre_flood_rate(event_data, "ER", "Event", "bt")
  
  # Get the first modelable day after the flood event
  first_modelable_day <- get_first_modelable_day(event_data, "date", "Event")
  
  # Filter data for the post-flood period
  post_flood_data <- df %>%
    filter(site_name == event_info$site_name &
             WaterYear == event_info$WaterYear &
             date >= first_modelable_day)
  
  # Calculate reduction percentage and response ratios for GPP and ER
  post_flood_data <- post_flood_data %>%
    mutate(
      Reduction_percent_GPP = ifelse(!is.na(pre_flood_GPP), 1 - (GPP / pre_flood_GPP) * 100, NA),
      Reduction_percent_ER = ifelse(!is.na(pre_flood_ER), 1 - (ER / pre_flood_ER) * 100, NA),
      RateRR_GPP = ifelse(!is.na(pre_flood_GPP), GPP / pre_flood_GPP, NA),
      RateRR_ER = ifelse(!is.na(pre_flood_ER), ER / pre_flood_ER, NA)
    )
  
  # Store the results for the current site
  all_sites_results[[i]] <- post_flood_data
}

# Combine results for all sites into one data frame if needed
combined_results <- bind_rows(all_sites_results)

# Ensure that the combined_results is ordered correctly
combined_results <- combined_results %>%
  arrange(site_name, WaterYear, wy_DOY)

# Display the structure of the combined results
str(combined_results)

