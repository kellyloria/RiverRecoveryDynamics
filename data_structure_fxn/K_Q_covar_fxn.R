K_Q_covar_metab <- function(df) {
  
  # define the number of sites
  sites <- unique(df$site)
  
  # define years
  years <- unique(df %>% 
                    dplyr::mutate(year = lubridate::year(date)) %>% 
                    dplyr::pull(year))
  
  # create a data frame to loop over
  df_summary <- data.frame(site = rep(sites, each = length(years)),
                           year = rep(years, times = length(sites)),
                           cov_K_Q = NA_real_,
                           n_obs = NA_integer_)
  
  # loop over site-year combinations
  for (i in 1:nrow(df_summary)) {
    
    # extract site and year
    site <- df_summary$site[i]
    year <- df_summary$year[i]
    
    # filter data
    df_use <- df %>%
      dplyr::mutate(year = lubridate::year(date)) %>%
      dplyr::filter(site == !!site, year == !!year) %>%
      dplyr::select(K600, discharge) %>%
      dplyr::filter(!is.na(K600), !is.na(discharge))
    
    # if not enough data, skip
    if (nrow(df_use) < 2) next
    
    # compute covariance
    df_summary$cov_K_Q[i] <- cor(df_use$K600, df_use$discharge)
    df_summary$n_obs[i] <- nrow(df_use)
  }
  
  # return summary dataframe
  return(df_summary)
}
