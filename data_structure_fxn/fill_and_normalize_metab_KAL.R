### MODIFIED By KAL to also infill NAs in ER
fill_and_normalize_metab <- function(df) {
  
  # source in functions
  source('/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/doi_10_5061_dryad_bcc2fqzj2__v20230622/code_fxns/fillMiss3.R')
  source('/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/doi_10_5061_dryad_bcc2fqzj2__v20230622/code_fxns/znorm.R')
  
  # define the number of sites
  sites <- unique(df$site)
  
  # and years
  years <- unique(df %>% 
                    dplyr::mutate(year = lubridate::year(date)) %>% 
                    dplyr::pull(year))
  
  # create a data frame to track missing data
  df_NA <- data.frame(site = rep(sites, each = length(years)),
                      year = rep(years, times = length(sites)))
  
  # where to save output data
  out <- list()
  
  # build for loop around site-year
  for(i in 1:nrow(df_NA)) {
    
    # extract row details
    deets <- df_NA[i,]
    
    site <- deets['site']
    year <- deets['year']
    
    # filter data for site-year
    df_use <- df %>% 
      dplyr::mutate(year = lubridate::year(date)) %>% 
      dplyr::filter(site %in% !!site,
                    year %in% !!year)
    
    # skip if empty
    if(nrow(df_use) == 0){
      df_NA[i, c('start_date_GPP', 'per_NA_GPP', 'start_date_ER', 'per_NA_ER')] <- NA
      df_NA[i, 'complete_year'] <- FALSE
      next
    }
    
    # pad data
    df_pad <- df_use %>% 
      padr::pad(start_val = lubridate::date(glue::glue("{year}-01-01")),
                end_val = lubridate::date(glue::glue("{year}-12-31")),
                'day') %>% 
      tidyr::fill(site, resolution, year, .direction = 'downup')
    
    # Get start and end dates for GPP
    start_date_GPP <- df_pad %>% dplyr::filter(!is.na(GPP)) %>% dplyr::pull(date) %>% dplyr::first()
    end_date_GPP <- df_pad %>% dplyr::filter(!is.na(GPP)) %>% dplyr::pull(date) %>% dplyr::last()
    
    # Get start and end dates for ER
    start_date_ER <- df_pad %>% dplyr::filter(!is.na(ER)) %>% dplyr::pull(date) %>% dplyr::first()
    end_date_ER <- df_pad %>% dplyr::filter(!is.na(ER)) %>% dplyr::pull(date) %>% dplyr::last()
    
    # Compute missing data percentage
    per_NA_GPP <- round(sum(is.na(df_pad$GPP)) / length(df_pad$GPP) * 100, 2)
    per_NA_ER <- round(sum(is.na(df_pad$ER)) / length(df_pad$ER) * 100, 2)
    
    # Save in tracking df
    df_NA[i, c('start_date_GPP', 'per_NA_GPP', 'start_date_ER', 'per_NA_ER')] <- c(start_date_GPP, per_NA_GPP, start_date_ER, per_NA_ER)
    
    # Check if the first date of GPP and ER is within the first 90 days
    if(!start_date_GPP %in% seq(df_pad$date[1], df_pad$date[90], by = 1) ||
       !end_date_GPP %in% seq(df_pad$date[nrow(df_pad) - 90], df_pad$date[nrow(df_pad)], by = 1) ||
       !start_date_ER %in% seq(df_pad$date[1], df_pad$date[90], by = 1) ||
       !end_date_ER %in% seq(df_pad$date[nrow(df_pad) - 90], df_pad$date[nrow(df_pad)], by = 1)){
      df_NA[i, 'complete_year'] <- FALSE
      next
    } else {
      df_NA[i, 'complete_year'] <- TRUE
    }
    
    # Gap filling
    if(per_NA_GPP < 40){
      df_pad$GPP_filled <- fillMiss3(dataset = data.frame(df_pad), 
                                     var = 'GPP', 
                                     block = 125,
                                     pmiss = 40, 
                                     plot = FALSE)
    } else {
      df_pad$GPP_filled <- df_pad$GPP  # Retain original values if not filled
    }
    
    if(per_NA_ER < 40){
      df_pad$ER_filled <- fillMiss3(dataset = data.frame(df_pad), 
                                    var = 'ER', 
                                    block = 125,
                                    pmiss = 40, 
                                    plot = FALSE)
    } else {
      df_pad$ER_filled <- df_pad$ER  # Retain original values if not filled
    }
    
    # Store result
    out[[i]] <- df_pad
  }  # <-- Added missing closing bracket for for-loop
  
  # Save tracker data
  readr::write_csv(df_NA,
                   '/Users/kellyloria/Documents/River_Recovery_Dynamics_Analysis/data/doi_10_5061_dryad_bcc2fqzj2__v20230622/data_citation/summary_site-year_NAs.csv')
  
  # Return compiled object
  if (length(out) > 0) {
    out <- do.call(rbind, out)
  } else {
    out <- NULL  # Prevents error if out is empty
  }
  
  return(out)
}

