# start with script RR_TS_analysis.R data object: dat_ts
str(dat_ts)

dat_ts$ER_abs <- c(dat_ts$ER*-1)

# List of specific site IDs to loop through
WsID <- c("nwis_01480617", "nwis_01608500", "nwis_0165389205", "nwis_01654500",
          "nwis_03067510", "nwis_03512000", "nwis_04087119", "nwis_04213500",
          "nwis_05406457", "nwis_05406500", "nwis_06893970", "nwis_07075250",
          "nwis_07075270", "nwis_07143672", "nwis_07191222", "nwis_08070200",       
          "nwis_09406000", "nwis_385446094430700", "nwis_385520094420000","nwis_385608094380300")


dat_tsH <- dat_ts %>%
  filter(site_name %in% WsID)%>%
  filter(!is.na(SPC))%>%
  dplyr::select(site_name, date, SPC, temp, Q)%>%
  dplyr::rename(WatershedID=site_name, DateTime=date, SpC=SPC, Temp=temp, Q=Q)

dat_split <- split(dat_tsH, dat_tsH$WatershedID)

str(dat_split)

##################################
## Base vs. Storm Flow separation
#################################

# This function performs baseflow separation for each site.
bs_separation <- function(x){
  site <- x
  Q_ts <- site$Q
  
  # Apply baseflow separation using a filter parameter
  bfs <- BaseflowSeparation(Q_ts, filter_parameter = 0.925, passes = 3)
  
  # Remove very small quickflow values to clean up the data
  bfs$qft[bfs$qft<0.001] <- 0
  bfs$DateTime <- site$DateTime
  
  # Merge the baseflow separation results with the original dataset
  combined <- merge(bfs, site, by="DateTime")
  
  # Add a column to mark if a given time is baseflow ('y') or not ('n')
  combined$base <- ifelse(combined$bt > combined$qft, yes="y", no="n")
  
  return(combined)
}

# Apply the baseflow separation function to each site in H_split
comb_list <- lapply(dat_split, function(x) bs_separation(x))

# Plot an example segment of data to visualize the baseflow separation
ggplot(comb_list$nwis_01480617, aes(DateTime, bt))+geom_line()+
  geom_line(aes(DateTime, qft), color="red")+
  geom_line(aes(DateTime, Q), color="blue")

###################################
##### Split storms and fit models

storm_fit <- function(x) {
  
  ## Split the data into baseflow and non-baseflow segments
  s <- split(x, f = x$base)
  df.base <- s$y  # Select only the baseflow data
  
  # Calculate the difference in days between consecutive measurements
  df.base$lag_diff <- as.numeric(difftime(df.base$DateTime, lag(df.base$DateTime), units = "days"))
  
  # Mark segments where the time gap is greater than 1 day
  df.base$splits <- ifelse(df.base$lag_diff > 1, yes = "y", no = "n")
  
  # Split the baseflow data into separate storms based on the day gaps
  sdf <- split(df.base, cumsum(1:nrow(df.base) %in% which(df.base$splits == "y")))
  
  ## Function to fit models to each storm's data
  fits <- function(df) {
    
    # Perform a linear regression of SpC against itself shifted by one day
    mod <- lm(df$SpC[-1] ~ df$SpC[-length(df$SpC)])
    
    # Extract the model coefficients and calculate additional metrics
    l <- as.data.frame(t(as.matrix(coefficients(mod))))
    colnames(l) <- c("beta", "alpha")
    
    # Calculate mu, which is beta divided by (1 - alpha)
    l$mu <- l$beta / (1 - l$alpha)
    
    # Record the number of data points, adjusted R-squared, and p-value
    l$n <- nrow(df)
    l$R2.adj <- summary(mod)$adj.r.squared
    l$p.val <- summary(mod)$coefficients[2, 4]
    
    return(l)
  }
  
  # Filter out storms with fewer than 12 observations
  sdf2 <- sdf[sapply(sdf, nrow) > 12]
  
  # Apply the fits function to each storm and combine the results
  if (length(sdf2) > 0) {
    c <- ldply(lapply(sdf2, function(x) fits(x)), data.frame)
    
    ## Filter out poor model fits
    cc <- c[-which(c$beta < 0),]        # Remove fits where beta is negative
    cc <- cc[-which(c$R2.adj < 0.95),]  # Remove fits with low adjusted R-squared
    cc <- cc[which(cc$alpha < 1),]      # Remove fits where alpha is 1 or greater
    
    # Re-index the results
    rownames(cc) <- seq(length = nrow(cc))
    
    # Keep only the storms that passed the filtering
    sdf2 <- sdf2[names(sdf2) %in% cc$.id]
    
    if (nrow(cc) > 0) {
      ## Add start and end dates for each storm
      cc$start_date <- NA; cc$end_date <- NA
      for (i in 1:nrow(cc)) {
        cc$start_date[i] <- as.character(sdf2[[i]][1,]$DateTime)
        cc$end_date[i] <- as.character(sdf2[[i]][nrow(sdf2[[i]]),]$DateTime)
      }
      
      # Convert the dates back to POSIXct format
      cc$start_date <- as.POSIXct(as.character(cc$start_date), format = "%Y-%m-%d")
      cc$end_date <- as.POSIXct(as.character(cc$end_date), format = "%Y-%m-%d")
      
      # Rename the ID column
      colnames(cc)[which(colnames(cc) == ".id")] <- "bs_id"
    } else {
      cc <- data.frame()  # Ensure cc is an empty data frame if no fits passed the filter
    }
  } else {
    cc <- data.frame()  # Handle case where sdf2 is empty
    sdf2 <- list()      # Ensure sdf2 is an empty list
  }
  
  # Return the filtered model fits and the corresponding time series
  output <- list(cc, sdf2)
  names(output) <- c("cc", "ts_recov")
  return(output)
}

# Apply the storm_fit function to all datasets in comb_list
storm_list <- lapply(comb_list, function(x) storm_fit(x))

############
############


## visualize
SC.mod_AR <- function(alpha, beta, df) {
  ## Data
  Ndays<-length(df$SpC)
  SC <- df$SpC
  ## Vectors for model output 
  pred_SC<-numeric(Ndays)
  pred_SC[1] <- SC[1]
  ## Process model
  for (j in 2:Ndays) {
    pred_SC[j] = alpha*pred_SC[j-1] + beta
  }
  return(pred_SC)
}



## visual check of outliers
WsID
WsT <- "nwis_01608500"
a <- storm_list[[WsT]]
b <- comb_list[[WsT]]

plot_grid(
  ggplot(a$cc, aes(start_date, beta))+geom_point()+
    coord_cartesian(ylim = c(0,75))+
    geom_ma(ma_fun = SMA, n=15, size=2, linetype = "solid", color="purple")+
    scale_x_datetime(limits = c(as.POSIXct(b$DateTime[1], format="%Y-%m-%d %H:%M:%S"),
                                as.POSIXct(b$DateTime[nrow(b)], format="%Y-%m-%d %H:%M:%S")))+
    theme(axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12))+
    labs(y=expression(beta)),
  ggplot(b, aes(DateTime, SpC))+geom_line(size=1, color="red")+
    coord_cartesian(ylim = c(0,4500))+
    theme(axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12))+
    labs(y="SpC (?S/cm)"),
  ggplot(b, aes(DateTime, Q))+geom_line(size=1, color="blue")+
    coord_cartesian(ylim = c(0,100))+
    theme(axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12))+
    labs(y="Q (cms)"),
  ncol=1, align = "hv"
)


i <- 4 #  i = the index of a particular storm event or baseflow segment

plot(a$ts_recov[[i]]$SpC,
     ylab="Specific Conductance (?S/cm)",
     xlab="Time Points Post QF<BF",
     ylim=c(100,500),cex=1) ## beta is max rate of increase
lines(SC.mod_AR(alpha = a$cc$alpha[i], beta= a$cc$beta[i], a$ts_recov[[i]]),
      col="#588C60",lwd = 3)
a$cc[i,]

#####################################
## All together
#####################################

all.cc <- ldply(lapply(storm_list, function(x) return(x$cc)), data.frame)

combined <- comb_list$nwis_04087119

## Visualize
plot_grid(
  ggplot(combined, aes(DateTime, bt))+geom_line()+
    geom_line(aes(DateTime, qft), color="red")+
    geom_line(aes(DateTime, Q), color="blue")+
    labs(y="Q (cms)"),
  
  ggplot(combined, aes(DateTime, SpC))+geom_line()+
    labs(y="SpC (uS/cm)"),
  
  ggplot(combined, aes(DateTime, base))+geom_point(),
  ncol=1,align="hv")

###########
###########
##########

###### exploring events as flushing or dilluting 
ggplot(combined, aes(Q, SpC, shape=base, color=base))+geom_point()

# Combine all data frames in the list into one data frame
combined_data <- bind_rows(comb_list, .id = "site")

water_year <- function(data) {
  data %>%
    mutate(date = ymd(date)) %>%
    mutate(WaterYear = if_else(month(date) >= 10, year(date) + 1, year(date)))
}

combined_data$date <- as.Date(combined_data$DateTime)
combined_data <-water_year(combined_data)

ggplot(combined_data, aes(x = log(Q+1), y = log(SpC+1), shape = base, color = WaterYear)) +
  geom_point() +
  facet_wrap(~site) +
  labs(x = "log Q (cms)",
       y = "log SPC (uScm)") +
  theme_minimal()

dat_tsH <- dat_ts %>%
  filter(!is.na(SPC))%>%
  dplyr::rename(WatershedID=site_name, DateTime=date, Q_ts=Q)

## Rejoin other data ##
whole_dat <- combined_data %>%
  full_join(dat_tsH, by=join_by(WatershedID, DateTime))

whole_dat_F <- whole_dat%>%filter(!is.na(turbidity_FNU))

ggplot(whole_dat%>%filter(!is.na(turbidity_FNU)), aes(x = log(Q+1), y = log(turbidity_FNU+1), shape = base, color = WaterYear)) +
  geom_point() +
  facet_wrap(~site) +
  labs(x = "log Q (cms)",
       y = "FNU") +
  theme_minimal()

ggplot(whole_dat%>%filter(!is.na(turbidity_FNU)), aes(x = log(Q+1), y = GPP, shape = base, color = WaterYear)) +
  geom_point() +
  facet_wrap(~site) +
  labs(x = "log Q (cms)",
       y = "GPP") +
  theme_minimal()

###

# List of specific site IDs to loop through
WsID <- c("nwis_01480617", "nwis_01608500", "nwis_0165389205", "nwis_01654500",
          "nwis_03067510", "nwis_03512000", "nwis_04087119", "nwis_04213500",
          "nwis_05406457", "nwis_05406500", "nwis_06893970", "nwis_07075250",
          "nwis_07075270", "nwis_07143672", "nwis_07191222", "nwis_08070200",       
          "nwis_09406000", "nwis_385446094430700", "nwis_385520094420000","nwis_385608094380300")

##
##
##

# Initialize an empty list to store all data
combined_data <- list()

# Loop through each site ID
for (site_id in WsID) {
  a <- storm_list[[site_id]]
  b <- comb_list[[site_id]]
  
  # Check if the site has valid data in storm_list
  if (!is.null(a)) {
    # Loop through each segment within the site
    for (i in seq_along(a$ts_recov)) {
      # Create a data frame with the necessary information
      df <- data.frame(
        SpC = a$ts_recov[[i]]$SpC,
        Fitted_SpC = SC.mod_AR(alpha = a$cc$alpha[i], beta = a$cc$beta[i], a$ts_recov[[i]]),
        Segment = i,
        Site = site_id,
        Time_Point = 1:length(a$ts_recov[[i]]$SpC)
      )
      
      # Append to the combined_data list
      combined_data[[length(combined_data) + 1]] <- df
    }
  }
}

# Combine all data into a single data frame
combined_df <- do.call(rbind, combined_data)

# Plot all data on one panel
SPC_storm_seg <- ggplot(combined_df, aes(x = Time_Point, y = SpC, color = Site)) +
  geom_line(aes(y = Fitted_SpC), color="black", size = 1) +  # Fitted line
  geom_point(size = 1.5) +  # Observed data
  labs(x = "Time points post QF<BF", y = "SpC (µS/cm)", 
       title = "Storm segments for different sites") +
  theme_bw() +
  facet_wrap(Segment~ Site, scales = "free_y")  # Facet by site, adjust y-axis individually

# ggsave(plot = SPC_storm_seg, filename = paste("./figures/SPC_storm_seg.png",sep=""),width=10,height=7.5,dpi=300)


###
## Visualize
plot_grid(
  ggplot(whole_dat, aes(DateTime, bt))+geom_line()+
    geom_line(aes(DateTime, qft), color="red")+
    geom_line(aes(DateTime, Q), color="blue")+
    labs(y="Q (cms)"),
  
  ggplot(whole_dat, aes(DateTime, SpC))+geom_line()+
    labs(y="SpC (uS/cm)"),
  
  ggplot(whole_dat, aes(DateTime, base))+geom_point(),
  ncol=1,align="hv")





















############
##
## FNU
##
FNU.mod_AR <- function(alpha, beta, df) {
  ## Data
  Ndays<-length(df$FNU)
  FNU <- df$FNU
  ## Vectors for model output 
  pred_FNU<-numeric(Ndays)
  pred_FNU[1] <- FNU[1]
  ## Process model
  for (j in 2:Ndays) {
    pred_FNU[j] = alpha*pred_FNU[j-1] + beta
  }
  return(pred_FNU)
}


##
##

# Initialize an empty list to store all data
combined_data <- list()

# Loop through each site ID
for (site_id in WsID) {
  a <- storm_list[[site_id]]
  b <- comb_list[[site_id]]
  
  # Check if the site has valid data in storm_list
  if (!is.null(a)) {
    # Loop through each segment within the site
    for (i in seq_along(a$ts_recov)) {
      # Create a data frame with the necessary information
      df <- data.frame(
        SpC = a$ts_recov[[i]]$SpC,
        Fitted_FNU = FNU.mod_AR(alpha = a$cc$alpha[i], beta = a$cc$beta[i], a$ts_recov[[i]]),
        Segment = i,
        Site = site_id,
        Time_Point = 1:length(a$ts_recov[[i]]$SpC)
      )
      
      # Append to the combined_data list
      combined_data[[length(combined_data) + 1]] <- df
    }
  }
}

# Combine all data into a single data frame
combined_df <- do.call(rbind, combined_data)

# Plot all data on one panel
SPC_storm_seg <- ggplot(combined_df, aes(x = Time_Point, y = SpC, color = Site)) +
  geom_line(aes(y = Fitted_SpC), color="black", size = 1) +  # Fitted line
  geom_point(size = 1.5) +  # Observed data
  labs(x = "Time points post QF<BF", y = "SpC (µS/cm)", 
       title = "Storm segments for different sites") +
  theme_minimal() +
  facet_wrap(Segment~ Site, scales = "free_y")  # Facet by site, adjust y-axis individually














# Neat code to look at SPC and road density:
# all.cc$WatershedID <- revalue(all.cc$.id, replace = c("WMT" = "Low Ref",
#                                                       "W22" = "R1",
#                                                       "W26" = "R2",
#                                                       "W79" = "R3",
#                                                       "W130" = "R4",
#                                                       "W207" = "R5",
#                                                       "W13" = "R6",
#                                                       "WGC" = "High Ref"))
# all.cc <- within(all.cc, WatershedID <- factor(WatershedID, levels=c("Low Ref","R1","R2","R3","R4","R5","R6","High Ref")))
# 
# ggplot(all.cc, aes(beta, color=WatershedID, fill=WatershedID))+
#   geom_histogram(alpha=0.3, bins = 50)+
#   scale_x_continuous(limits=c(0.01,75))+
#   facet_wrap(~WatershedID, nrow = 2)
# #theme(legend.position = "none")
# 
# wa_col_ramp <- colorRampPalette(c("#FEF8DE", "#57061C"),)
# wa_col <- wa_col_ramp(8)
# scaleFUN <- function(x) sprintf("%.1f", x)
# 
# ggplot(all.cc, aes(y=beta,x=WatershedID, fill=WatershedID))+
#   geom_boxplot()+
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   scale_y_continuous(trans = "log",breaks=c(0.05,0.5,2,5,20,50,200),labels=comma)+
#   labs(y=expression(beta), x="Watershed ID by Road Density")+
#   theme(axis.title = element_text(size=16), axis.text = element_text(size=14),
#         axis.text.x = element_text(angle=45, hjust=1))+
#   scale_fill_manual(name="Road Density",values = c("Low Ref"=wa_col[1],"R1"=wa_col[2],"R2"=wa_col[3],"R3"=wa_col[4],
#                                                    "R4"=wa_col[5],"R5"=wa_col[6],"R6"=wa_col[7],"High Ref"=wa_col[8]))
# 
# 
# ggplot(all.cc, aes(y=mu,x=WatershedID, fill=WatershedID))+
#   geom_boxplot()+
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   scale_y_continuous(trans = "log",breaks=c(50,200,500),labels=comma,
#                      limits=c(50,500))+
#   labs(y=expression(mu), x="Watershed ID by Road Density")+
#   theme(axis.title = element_text(size=16), axis.text = element_text(size=14),
#         axis.text.x = element_text(angle=45, hjust=1))+
#   scale_fill_manual(name="Road Density",values = c("Low Ref"=wa_col[1],"R1"=wa_col[2],"R2"=wa_col[3],"R3"=wa_col[4],
#                                                    "R4"=wa_col[5],"R5"=wa_col[6],"R6"=wa_col[7],"High Ref"=wa_col[8]))
# 
# 
# plot_grid(
#   
#   ggplot(all.cc, aes(beta, color=.id, fill=.id))+
#     geom_histogram(alpha=0.3, bins = 50)+
#     facet_wrap(~.id,ncol = 1)+
#     theme(legend.position = "none"),
#   
#   ggplot(all.cc, aes(mu, color=.id, fill=.id))+
#     geom_histogram(alpha=0.3, bins = 50)+
#     facet_wrap(~.id,ncol = 1),
#   
#   ncol=2
# )
# 
# ggplot(all.cc, aes(beta, mu, color=.id))+
#   geom_point()

###

## for 