##==============================================================================
## Script for selecting two years of semi-continuous GPP data
## Code author: J.R. Blaszczak
##==============================================================================

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

#  setwd("~/EcoProject_test")

############################
## To create linked file
############################

## Import site data from Appling
# https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982
site <- fread("./data/site_data.tsv")
colnames(site)

## Import StreamLight from Savoy
# https://www.sciencebase.gov/catalog/item/5f974adfd34e198cb77db168
SL <- read.delim("./data/SL_SiteInfo.txt")

colnames(SL)[colnames(SL) == "Site_ID"] <- "site_name"

## Secondary stream order source from hypoxia data set
# https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5
hyp <- read.csv("./data/Distribution, frequency, and global extent of hypoxia in rivers/Data/GRDO_GEE_HA_NHD.csv")

hyp <- hyp[which(hyp$DB_Source == "PC"), c("SiteID","ORD_STRA","NHD_STREAMORDE")]
colnames(hyp)[which(colnames(hyp) == "SiteID")] <- "site_name"

## Merge
df <- left_join(site, SL, "site_name")
df <- left_join(df, hyp, "site_name")
colnames(df)

## subset
sub <- df[,c("site_name","long_name","StreamOrde",
             "site_type","struct.canal_flag","struct.dam_flag","struct.npdes_flag",
             "ORD_STRA","NHD_STREAMORDE")]


##################################################
## Select sites based on data quality
##################################################

# Import and subset model diagnostics
diagnostics <- read.table("./data/diagnostics.tsv",sep = "\t", header=T)

diagnostics <- diagnostics[which(diagnostics$site %in% sub$site_name),]
highq_sites <- diagnostics[which(diagnostics$K600_daily_sigma_Rhat < 1.05 & 
                                   diagnostics$err_obs_iid_sigma_Rhat < 1.05 &
                                   diagnostics$err_proc_iid_sigma_Rhat < 1.05 &
                                   diagnostics$neg_GPP < 20 & diagnostics$pos_ER),] #229
highq_site_names <- unique(highq_sites$site) ## 208


# Subset s based on high sites and site type and flags
s <- sub[which(sub$site_name %in% highq_site_names),] ## 208
s <- s[which(s$site_type == "ST"),] ## 204
s <- s[which(s$struct.dam_flag %in% c(NA,"80","95")),] ## 104
# which have light from Phil
#s_l <- s[!is.na(s$StreamOrde),] # 41

# Import time series
NWIS <- read.table("./data/daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- as.POSIXct(as.character(NWIS$date), format="%Y-%m-%d")
head(NWIS)

NWIS <- NWIS[which(NWIS$GPP.Rhat < 1.05),]
NWIS <- NWIS[which(NWIS$K600.Rhat < 1.05),]
length(levels(factor(NWIS$site_name))) ## 355


## Subset columns and sites
NWIS_sub <- NWIS[,c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                    "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                    "temp.water","discharge","shortwave","velocity")]
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                        "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                        "temp","Q","light","velocity")

## Subset to sites in high_sites (sites with high confidence rating and limited dam interference)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]
# Confirm
length(levels(as.factor(NWIS_sub$site_name))) ## 104

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  count()
## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1]))
maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)

hist(dat_per_year$year)


## subset for sites with a max gap of 14 days
sub_by_gap <- maxgap[which(maxgap$gap <= 14),]

length(levels(as.factor(sub_by_gap$site_name))) ## 92

## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))

## at least 275 days per year
# sub_by_gap <- sub_by_gap[which(sub_by_gap$n <= 275),]

sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum    #sub_by_gap_sum[which(sub_by_gap_sum$gap <= 2),] # 41 when tight dam restriction

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),] ## only sites with two or more years

## Subset to years that meet criteria
sub_by_gap$site_year <- paste(sub_by_gap$site_name,sub_by_gap$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% sub_by_gap$site_year),]
TS_site <- s[which(s$site_name %in% high_q$site_name),]

## Attach the median GPP
TS$GPP_temp <- TS$GPP
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-3):exp(-3), 1)
TS$ER_temp <- TS$ER
TS[which(TS$ER > 0),]$ER_temp <- sample(-exp(-3):-exp(-3), 1)

TS_dat <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = c("GPP_temp", "ER_temp", "K600"), .funs = c(mean, max, min))
colnames(TS_dat) <- c("site_name","GPP_mean", "ER_mean","K600_mean",
                      "GPP_max","ER_max","K600_max",
                      "GPP_min","ER_min","K600_min")
TS_site <- left_join(TS_site, TS_dat, by="site_name")


TS_df <- left_join(TS, TS_site, by="site_name")

## Assign a stream order classification
TS_df$order_group <- "NA"
TS_df[which(TS_df$NHD_STREAMORDE %in% c(1,2)),]$order_group <- "small"
TS_df[which(TS_df$NHD_STREAMORDE %in% c(3,4)),]$order_group <- "mid"
TS_df[which(TS_df$NHD_STREAMORDE >= 5),]$order_group <- "large"
TS_df[which(TS_df$NHD_STREAMORDE %in% c(NA)),]$order_group <- "unknown"

## validate years that meet criteria
## count days per year
TSdat_per_year <- TS %>%
  group_by(site_name, year) %>%
  count()
## identify the max day gap per year
TSgap_per_year <- TS %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1]))
TSmaxgap <- TSgap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)
TS_valid_years <- merge(TSdat_per_year, TSmaxgap, by=c("site_name","year"))

TS_valid_yearsq <-TS_valid_years %>%
  filter(n>270, gap < 12)

# Filter TS_df based on the filtered_valid_yearsq
filtered_TS_df <- TS_df %>%
  semi_join(TS_valid_yearsq, by = c("site_name", "year"))

## choose sites from different groups that meet criteria
View(filtered_TS_df[which(filtered_TS_df$order_group == "small"),])
View(filtered_TS_df[which(filtered_TS_df$order_group == "mid"),])
View(filtered_TS_df[which(filtered_TS_df$order_group == "large"),])

ggplot(filtered_TS_df%>%filter(order_group == "small"), aes(date, GPP_temp, color=site_name)) +
  geom_line() +theme_bw()

ggplot(filtered_TS_df%>%filter(order_group == "mid"), aes(date, GPP_temp, color=site_name)) +
    geom_line() +theme_bw()

ggplot(filtered_TS_df%>%filter(order_group == "large"), aes(date, GPP_temp, color=site_name)) +
  geom_line() +theme_bw()


###########################
## Export
###########################
# saveRDS(filtered_TS_df, "./data/Metab_TS.rds")


##############################
## Plot for talks
###########################
site_subset_list <- split(site_subset, site_subset$site_name)


df <- site_subset_list$nwis_01649190
ratio_QL <- max(df$light)/max(df$Q)
GPP_plot <- ggplot(df, aes(date, GPP_temp))+
  #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
  geom_point(color="chartreuse4", size=2)+geom_line(color="chartreuse4", size=1)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")

data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
  geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
  geom_line(size=1, color="deepskyblue4")+
  scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
  labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
        axis.title.y.right = element_text(size=12, color="deepskyblue4"),
        axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")

GPP_plot + data_plot + plot_layout(ncol = 1)

ggplot(df, aes(light, GPP))+
  geom_point(size=1.5)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily PPFD")+
  theme_bw()+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))

