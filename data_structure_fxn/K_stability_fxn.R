K600_stability <- function(df) {
  
  # ensure date is a Date type
  df$date <- as.Date(df$date)
  
  # get unique sites and years
  sites <- unique(df$site)
  years <- unique(lubridate::year(df$date))
  
  # create output dataframe
  df_summary <- data.frame(site = rep(sites, each = length(years)),
                           year = rep(years, times = length(sites)),
                           slope_K600 = NA_real_,
                           n_obs = NA_integer_)
  
  # loop over site-year combos
  for (i in 1:nrow(df_summary)) {
    
    site <- df_summary$site[i]
    year <- df_summary$year[i]
    
    # filter to this site-year
    df_use <- df %>%
      dplyr::filter(site == !!site,
                    lubridate::year(date) == !!year) %>%
      dplyr::filter(!is.na(K600))
    
    # skip if too few points
    if (nrow(df_use) < 2) next
    
    # convert date to numeric for regression (e.g., days since start of year)
    df_use <- df_use %>%
      dplyr::mutate(day_num = as.numeric(date - as.Date(paste0(year, "-01-01"))))
    
    # linear model: K600 ~ day_num
    mod <- lm(K600 ~ day_num, data = df_use)
    
    # extract slope
    df_summary$slope_K600[i] <- coef(mod)["day_num"]
    df_summary$n_obs[i] <- nrow(df_use)
  }
  
  return(df_summary)
}
