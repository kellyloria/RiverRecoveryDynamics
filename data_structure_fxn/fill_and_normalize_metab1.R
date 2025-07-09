fill_and_normalize_metab1 <- function(df) {
  
  # source in functions
  source('/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/RiverRecoveryDynamics/data_structure_fxn/fillMiss3.R')
  source('/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/RiverRecoveryDynamics/data_structure_fxn/znorm.R')
  
  # define the number of sites and years
  sites <- unique(df$site)
  years <- unique(df %>% dplyr::mutate(year = lubridate::year(date)) %>% dplyr::pull(year))
  
  # create summary tracking dataframe with correctly typed NA values
  df_NA <- data.frame(
    site = rep(sites, each = length(years)),
    year = rep(years, times = length(sites)),
    total_days = as.integer(NA),
    num_NA_before = as.integer(NA),
    num_NA_after = as.integer(NA),
    num_infilled = as.integer(NA),
    start_date = as.Date(NA),
    per_NA = as.numeric(NA),
    complete_year = as.logical(NA)
  )
  
  # storage for filled output
  out <- list()
  
  for (i in 1:nrow(df_NA)) {
    
    # get site-year
    deets <- df_NA[i, ]
    site <- deets$site
    year <- deets$year
    
    # filter to site-year
    df_use <- df %>%
      dplyr::mutate(year = lubridate::year(date)) %>%
      dplyr::filter(site == !!site, year == !!year)
    
    if (nrow(df_use) == 0) {
      df_NA[i, c("start_date", "per_NA", "complete_year")] <- list(NA, NA, FALSE)
      next
    }
    
    # pad to full year
    df_pad <- df_use %>%
      padr::pad(start_val = lubridate::date(glue::glue("{year}-01-01")),
                end_val = lubridate::date(glue::glue("{year}-12-31")),
                interval = "day") %>%
      tidyr::fill(site, resolution, year, .direction = "downup")
    
    # record total days
    df_NA[i, "total_days"] <- nrow(df_pad)
    
    # count NAs before filling
    NAs_before <- sum(is.na(df_pad$GPP))
    df_NA[i, "num_NA_before"] <- NAs_before
    per_NA <- round(NAs_before / nrow(df_pad) * 100, 2)
    df_NA[i, "per_NA"] <- per_NA
    
    # find first/last non-NA dates
    start_date <- df_pad %>% dplyr::filter(!is.na(GPP)) %>% dplyr::pull(date) %>% dplyr::first()
    end_date <- df_pad %>% dplyr::filter(!is.na(GPP)) %>% dplyr::pull(date) %>% dplyr::last()
    df_NA[i, "start_date"] <- start_date
    
    # assess if year is complete
    if (!start_date %in% seq(df_pad$date[1], df_pad$date[90], by = 1) ||
        !end_date %in% seq(df_pad$date[nrow(df_pad) - 90], df_pad$date[nrow(df_pad)], by = 1)) {
      df_NA[i, "complete_year"] <- FALSE
      next
    } else {
      df_NA[i, "complete_year"] <- TRUE
    }
    
    # gap-fill if acceptable NA%
    if (per_NA < 40) {
      df_pad$GPP_filled <- fillMiss3(dataset = data.frame(df_pad), 
                                     var = 'GPP', 
                                     block = 125,
                                     pmiss = 40, 
                                     plot = FALSE)
      
      NAs_after <- sum(is.na(df_pad$GPP_filled))
      df_NA[i, "num_NA_after"] <- NAs_after
      df_NA[i, "num_infilled"] <- NAs_before - NAs_after
      
      out[[i]] <- df_pad
    } else {
      df_NA[i, c("num_NA_after", "num_infilled")] <- list(NA, NA)
    }
    
  } # end for loop
  
  # write the full df_NA to disk
  readr::write_csv(df_NA,
                   '/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/doi_10_5061_dryad_bcc2fqzj2__v20230622/data_citation/summary_site-year_NAs.csv')
  
  # return the compiled object
  out <- do.call(rbind, out)
  return(out)
}
