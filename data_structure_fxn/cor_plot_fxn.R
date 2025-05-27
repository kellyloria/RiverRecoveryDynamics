
cor_plot_K600_ER_discharge <- function(df) {
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(patchwork)
  
  # Ensure correct date format and add year column
  df <- df %>%
    mutate(date = as.Date(date),
           year = year(date))
  
  sites <- unique(df$site)
  years <- unique(df$year)
  
  # Initialize storage
  cor_results <- data.frame()
  plot_list <- list()
  
  for (s in sites) {
    for (y in years) {
      
      df_sub <- df %>%
        filter(site == s, year == y) %>%
        filter(!is.na(K600), !is.na(ER), !is.na(discharge))
      
      if (nrow(df_sub) < 10) next  # skip sparse data
      
      # Calculate correlations
      cor_ER_K <- cor(df_sub$ER, df_sub$K600, use = "complete.obs")
      cor_K_discharge <- cor(df_sub$K600, df_sub$discharge, use = "complete.obs")
      
      # Record results
      cor_results <- rbind(cor_results, data.frame(
        site = s,
        year = y,
        cor_ER_K600 = round(cor_ER_K, 3),
        cor_K600_discharge = round(cor_K_discharge, 3),
        n_obs = nrow(df_sub)
      ))
      
      # Create plots
      p1 <- ggplot(df_sub, aes(x = K600, y = ER)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(title = paste("ER vs. K600"),
             subtitle = paste(s, y, "- r =", round(cor_ER_K, 2)),
             x = "K600", y = "ER") +
        theme_minimal()
      
      p2 <- ggplot(df_sub, aes(x = discharge, y = K600)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
        labs(title = paste("K600 vs. Discharge"),
             subtitle = paste(s, y, "- r =", round(cor_K_discharge, 2)),
             x = "Discharge", y = "K600") +
        theme_minimal()
      
      p3 <- ggplot(df_sub, aes(x = date, y = K600)) +
        geom_line(color = "black") +
        geom_point(size = 0.8, alpha = 0.5) +
        labs(title = "Daily K600",
             subtitle = paste(s, y),
             x = "Date", y = "K600") +
        theme_minimal()
      
      # Combine plots
      combined_plot <- (p1 | p2) / p3 +
        plot_annotation(title = paste("Site:", s, "Year:", y))
      
      plot_list[[paste0(s, "_", y)]] <- combined_plot
    }
  }
  
  return(list(
    correlation_results = cor_results,
    plots = plot_list
  ))
}
