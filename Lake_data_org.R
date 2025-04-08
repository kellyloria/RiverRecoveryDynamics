# stream data goes from  :2008-01-01 to  :2021-12-31 


corman_dat_metab <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/load_metab_results/load_metab_results.csv")%>%
  mutate(date=as.Date(DateTime),
         year=year(date),
         jday= yday(date))

corman_dat_metadf <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/load_metab_results/corman_lake_cords.csv")
corman_dat_metadf$source <- "corman_lake"


ggplot(corman_dat_metab, aes(x = date, y = GPP, color = factor(year), group = interaction(Lake, year))) +
  geom_line() +  
  geom_point() +  
  facet_wrap(~Lake, scales = 'free_y') +  # scales = 'free_y'
  theme_minimal() +  
  labs(
    title = "Average annual GPP by stream site",
    x = "Day of Year",
    y = "GPP"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    panel.grid.major = element_line(size = 0.5, color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  scale_color_viridis_d(name="Year")  # Use discrete color scale for 'year'


### Alright just one year per site .... 

# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-ntl.402.2
Ladwig_NTL_dat_metab <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/knb-lter-ntl.402.2/model_output.csv")%>%
  mutate(date=as.Date(datetime),
         year=year(date),
         jday= yday(date)) %>%
  filter(date> as.Date("2007-12-31") & date< as.Date("2022-01-01"))

names(Ladwig_NTL_dat_metab)
# NEP in milligramPerMeterCubedPerDay


ggplot(Ladwig_NTL_dat_metab, aes(x = jday, y = NEP, color = factor(year), group = interaction(lake, year))) +
  geom_line() +  
  geom_point() +  
  facet_wrap(~lake, scales = 'free_y') +  # scales = 'free_y'
  theme_minimal() +  
  labs(
    title = "Average annual GPP by stream site",
    x = "Day of Year",
    y = "GPP"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    panel.grid.major = element_line(size = 0.5, color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  scale_color_viridis_d(name="Year")  # Use discrete color scale for 'year'


Ladwig_NTL_dat_metadf <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/knb-lter-ntl.402.2/Ladwig_NTL_lake_cords.csv")
Ladwig_NTL_dat_metadf$source <- "Ladwig_NTL_lake"






Rabaey_dat_metab <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/edi.1531.2/metabolism_daily_summaries.csv")%>%
  mutate(date=as.Date(date),
         year=year(date),
         jday= yday(date)) # %>% filter(date> as.Date("2007-12-31") & date< as.Date("2022-01-01"))

names(Rabaey_dat_metab)

# https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1531.2

ggplot(Rabaey_dat_metab, aes(x = date, y = nep, color = factor(year), group = interaction(name, year))) +
  geom_line() +  
  geom_point() +  
  facet_wrap(~name, scales = 'free_y') +  # scales = 'free_y'
  theme_minimal() +  
  labs(
    title = "Average annual GPP by stream site",
    x = "Day of Year",
    y = "GPP"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    panel.grid.major = element_line(size = 0.5, color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  scale_color_viridis_d(name="Year")  # Use discrete color scale for 'year'

### NOT multi year 



Oleksy_dat_metab <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/edi.1541.2/daily_metabolism.csv") %>%
  mutate(date=as.Date(date),
         year=year(date),
         jday= yday(date)) # %>% filter(date> as.Date("2007-12-31") & date< as.Date("2022-01-01"))

names(Oleksy_dat_metab)

 # https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1541.2

ggplot(Oleksy_dat_metab, aes(x = date, y = GPP_mgO2L, color = factor(year), group = interaction(lakeName, year))) +
  geom_line() +  
  geom_point() +  
  facet_wrap(~lakeName, scales = 'free_y') +  # scales = 'free_y'
  theme_minimal() +  
  labs(
    title = "Average annual GPP by stream site",
    x = "Day of Year",
    y = "GPP"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    panel.grid.major = element_line(size = 0.5, color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  scale_color_viridis_d(name="Year")  # Use discrete color scale for 'year'

### woof some overlap - one multi-year 


# Europe: https://portal.edirepository.org/nis/mapbrowse?packageid=edi.991.1


# /Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/edi.1531.2/metabolism_daily_summaries.csv


unique(corman_dat$Lake)

dat_ts_clean <- dat_ts%>%
  filter(GPP.Rhat < 1.05 &  GPP >0 & ER <0)

# unique(dat_ts_clean$long_name)

marz_dat_names <- marz_dat%>%
  dplyr::select(Site.Number, Latitude, Longitude) %>%
  mutate(source="corman")

dat_ts_names <-metab_meta_dat %>%
  dplyr::select(Site.Number = "site_no", Latitude = "dec_lat_va", Longitude= "dec_long_va") %>%
  mutate(source="appling")

## quick site check
site_check <- rbind(marz_dat_names, dat_ts_names)



lotting_dat_metadf <- read.csv("/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/Lakes_dat/Facundo_Rec/knb-lter-ntl.397.1/lottig_etal_data.csv")
sites <-unique(lotting_dat_metadf$lake)
getwd()

